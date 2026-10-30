"""Vectorized particle advection through a steady 2D velocity field.

This is the numerical heart of the package and the direct replacement
for the legacy ``advect.m``.  Two things make it orders of magnitude
faster than the MATLAB original:

1. **Bilinear grid sampling instead of scattered interpolation.**
   ``advect.m`` evaluated a Delaunay-based ``TriScatteredInterp`` at
   every particle and every time step (its own TODO comment calls this
   "a major bottleneck").  Here the field lives on a regular grid, so
   each velocity lookup is a handful of array operations.
2. **All particles step together.**  The NumPy path advances the whole
   particle cloud per time step with vectorized arithmetic; with Numba
   installed, a parallel compiled kernel integrates each particle's
   full trajectory in one pass over threads.

The integrator is 4th-order Runge–Kutta by default (the legacy code
used forward Euler; pass ``integrator="euler"`` to reproduce it).

Boundary rule: velocity is zero outside the grid, so particles that
leave the domain freeze in place.  NaN velocities in the input field
were replaced by zero when the :class:`~lcs2d.fields.VelocityField`
was built, so masked cells (land, obstacles) behave the same way.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import _kernels
from .fields import VelocityField

__all__ = ["advect", "Trajectories", "resolve_backend"]


def resolve_backend(backend: str = "auto") -> str:
    """Turn a backend request into ``"numba"`` or ``"numpy"``.

    ``"auto"`` picks Numba when it is importable, else NumPy.
    Requesting ``"numba"`` without Numba installed raises with a hint.
    """
    if backend == "auto":
        return "numba" if _kernels.NUMBA_AVAILABLE else "numpy"
    if backend == "numba" and not _kernels.NUMBA_AVAILABLE:
        raise RuntimeError(
            "backend='numba' requested but Numba is not installed. "
            "Install it with: pip install numba  (or pip install lcs2d[fast])")
    if backend not in ("numpy", "numba"):
        raise ValueError(f"unknown backend {backend!r}; "
                         "expected 'auto', 'numpy' or 'numba'")
    return backend


@dataclass
class Trajectories:
    """Recorded particle paths.

    Attributes
    ----------
    t : ndarray, shape (n_saved,)
        Times of the saved snapshots (``t[0] = 0``).
    x, y : ndarray, shape (n_saved, n_particles)
        Particle positions at each saved snapshot.
    """

    t: np.ndarray
    x: np.ndarray
    y: np.ndarray

    def save(self, path: str) -> None:
        np.savez_compressed(path, t=self.t, x=self.x, y=self.y)

    @classmethod
    def load(cls, path: str) -> "Trajectories":
        with np.load(path, allow_pickle=False) as data:
            return cls(data["t"], data["x"], data["y"])


# ---------------------------------------------------------------------------
# NumPy path
# ---------------------------------------------------------------------------

def _sample_uv_numpy(u, v, x0, y0, inv_dx, inv_dy, px, py):
    """Bilinear velocity sampling of a particle cloud; zero outside the grid."""
    ny, nx = u.shape
    gx = (px - x0) * inv_dx
    gy = (py - y0) * inv_dy
    inside = (gx >= 0.0) & (gy >= 0.0) & (gx <= nx - 1.0) & (gy <= ny - 1.0)
    gx = np.clip(gx, 0.0, nx - 1.0)
    gy = np.clip(gy, 0.0, ny - 1.0)
    i = np.minimum(gx.astype(np.intp), nx - 2)
    j = np.minimum(gy.astype(np.intp), ny - 2)
    fx = gx - i
    fy = gy - j
    w00 = (1.0 - fx) * (1.0 - fy)
    w10 = fx * (1.0 - fy)
    w01 = (1.0 - fx) * fy
    w11 = fx * fy
    us = (w00 * u[j, i] + w10 * u[j, i + 1]
          + w01 * u[j + 1, i] + w11 * u[j + 1, i + 1])
    vs = (w00 * v[j, i] + w10 * v[j, i + 1]
          + w01 * v[j + 1, i] + w11 * v[j + 1, i + 1])
    us[~inside] = 0.0
    vs[~inside] = 0.0
    return us, vs


def _step_numpy(u, v, x0, y0, inv_dx, inv_dy, px, py, dt, use_rk4):
    if use_rk4:
        k1u, k1v = _sample_uv_numpy(u, v, x0, y0, inv_dx, inv_dy, px, py)
        k2u, k2v = _sample_uv_numpy(u, v, x0, y0, inv_dx, inv_dy,
                                    px + 0.5 * dt * k1u, py + 0.5 * dt * k1v)
        k3u, k3v = _sample_uv_numpy(u, v, x0, y0, inv_dx, inv_dy,
                                    px + 0.5 * dt * k2u, py + 0.5 * dt * k2v)
        k4u, k4v = _sample_uv_numpy(u, v, x0, y0, inv_dx, inv_dy,
                                    px + dt * k3u, py + dt * k3v)
        return (px + dt * (k1u + 2.0 * (k2u + k3u) + k4u) / 6.0,
                py + dt * (k1v + 2.0 * (k2v + k3v) + k4v) / 6.0)
    su, sv = _sample_uv_numpy(u, v, x0, y0, inv_dx, inv_dy, px, py)
    return px + dt * su, py + dt * sv


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _field_args(field: VelocityField):
    return (field.u, field.v, float(field.x[0]), float(field.y[0]),
            1.0 / field.dx, 1.0 / field.dy)


def _check_dt(dt: float | None, field: VelocityField,
              safety: float = 0.5) -> float:
    """Validate a time-step magnitude, defaulting to ``field.cfl_dt()``."""
    if dt is None:
        return field.cfl_dt(safety)
    dt = float(dt)
    if not np.isfinite(dt) or dt == 0.0:
        raise ValueError(f"time step dt must be finite and nonzero, got {dt}")
    return abs(dt)


def _plan_steps(t: float, dt: float | None, field: VelocityField):
    """Choose (n_steps, signed_dt) so that n_steps * dt == t exactly."""
    t = float(t)
    if not np.isfinite(t) or t == 0.0:
        raise ValueError("integration time t must be finite and nonzero "
                         "(negative t integrates backward in time)")
    dt = _check_dt(dt, field)
    n_steps = max(1, int(round(abs(t) / dt)))
    return n_steps, t / n_steps


def advect(field: VelocityField,
           px: np.ndarray,
           py: np.ndarray,
           t: float,
           dt: float | None = None,
           integrator: str = "rk4",
           backend: str = "auto",
           record: bool = False,
           save_every: int = 1):
    """Advect particles through a steady velocity field.

    Parameters
    ----------
    field : VelocityField
        The flow to advect through.
    px, py : ndarray
        Initial particle positions (any shape; flattened internally).
    t : float
        Total integration time.  Negative ``t`` integrates backward in
        time (used for attracting LCS / backward FTLE).
    dt : float, optional
        Time-step magnitude.  Default: ``field.cfl_dt()``.  The step is
        adjusted slightly so an integer number of steps lands exactly
        on ``t``.
    integrator : {"rk4", "euler"}
        4th-order Runge–Kutta (default) or forward Euler (the legacy
        MATLAB scheme, kept for comparison).
    backend : {"auto", "numpy", "numba"}
        "auto" uses the parallel Numba kernels when available.
    record : bool
        If True, return a :class:`Trajectories` object with saved
        snapshots instead of only the final positions.
    save_every : int
        When recording, keep every ``save_every``-th step (the initial
        and final positions are always kept).

    Returns
    -------
    (x, y) : tuple of ndarray
        Final positions, same shape as the input, when ``record=False``.
    Trajectories
        When ``record=True`` (positions have shape ``(n_saved, n)``).
    """
    if integrator not in ("rk4", "euler"):
        raise ValueError(f"unknown integrator {integrator!r}; "
                         "expected 'rk4' or 'euler'")
    use_rk4 = integrator == "rk4"
    backend = resolve_backend(backend)
    if save_every < 1:
        raise ValueError("save_every must be >= 1")

    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)
    if px.shape != py.shape:
        raise ValueError("px and py must have the same shape")
    if not (np.isfinite(px).all() and np.isfinite(py).all()):
        raise ValueError("particle positions px, py must all be finite")
    shape = px.shape
    fx = np.ascontiguousarray(px.ravel())
    fy = np.ascontiguousarray(py.ravel())

    n_steps, sdt = _plan_steps(t, dt, field)
    args = _field_args(field)

    if not record:
        if backend == "numba":
            ox, oy = _kernels.advect_final_kernel(*args, fx, fy, sdt,
                                                  n_steps, use_rk4)
        else:
            ox, oy = fx.copy(), fy.copy()
            for _ in range(n_steps):
                ox, oy = _step_numpy(*args, ox, oy, sdt, use_rk4)
        return ox.reshape(shape), oy.reshape(shape)

    # record trajectories -----------------------------------------------
    saved_steps = list(range(0, n_steps + 1, save_every))
    if saved_steps[-1] != n_steps:
        saved_steps.append(n_steps)
    times = np.asarray(saved_steps, dtype=np.float64) * sdt

    if backend == "numba":
        tx, ty = _kernels.advect_record_kernel(*args, fx, fy, sdt,
                                               n_steps, save_every, use_rk4)
    else:
        n = fx.size
        tx = np.empty((len(saved_steps), n))
        ty = np.empty((len(saved_steps), n))
        tx[0] = fx
        ty[0] = fy
        ox, oy = fx.copy(), fy.copy()
        k = 1
        for s in range(1, n_steps + 1):
            ox, oy = _step_numpy(*args, ox, oy, sdt, use_rk4)
            if s % save_every == 0 or s == n_steps:
                tx[k] = ox
                ty[k] = oy
                k += 1
    return Trajectories(times, tx, ty)
