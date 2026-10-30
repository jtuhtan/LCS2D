"""Finite-time Lyapunov exponent (FTLE) fields.

Method (auxiliary grid, after Shadden et al. 2005): each seed point is
surrounded by four auxiliary particles offset by a small distance
``delta`` (up / down / left / right — exactly the "secondary seeds" of
the legacy ``advect.m``).  All five clouds are advected for time ``T``,
the flow-map Jacobian is formed by central differences over the
auxiliary particles, and the FTLE is

    sigma = ln( lambda_max( J^T J ) ) / (2 |T|)

where ``J^T J`` is the (right) Cauchy–Green strain tensor.  Ridges of
``sigma`` mark the repelling (forward time) or attracting (backward
time, ``t < 0``) Lagrangian coherent structures.

Two deliberate corrections relative to the legacy ``ftle.m``:

* the central difference divides by ``2 * delta``; the MATLAB
  expression ``./2*delta`` actually *multiplied* by ``delta`` due to
  operator precedence, and
* the exponent uses the standard ``1 / (2|T|)`` normalization, where
  the MATLAB code used ``1 / T`` (twice the standard value).

Neither changes where the ridges are (the legacy script rescaled the
field to [0, 1] anyway), but absolute values now match the literature.
"""

from __future__ import annotations

import numpy as np

from .advect import _check_dt, advect
from .fields import VelocityField
from .result import LCSResult

__all__ = ["ftle_field", "seed_grid"]


def seed_grid(field: VelocityField,
              res: float | None = None,
              margin: float = 0.0):
    """Build a uniform seed grid covering the field's domain.

    Parameters
    ----------
    res : float, optional
        Seed spacing; defaults to the field's own grid spacing.
    margin : float
        Distance to inset the seeds from the domain boundary (useful
        for resampled scattered data whose edges are unreliable — the
        legacy ``seed.m`` used ``2 * resSeed``).
    """
    if res is None:
        res = max(field.dx, field.dy)
    res = float(res)
    if not (np.isfinite(res) and res > 0.0):
        raise ValueError(f"seed spacing res must be positive, got {res}")
    if not margin >= 0.0:
        raise ValueError(f"margin must be >= 0, got {margin}")
    xmin, xmax, ymin, ymax = field.extent
    xmin += margin
    xmax -= margin
    ymin += margin
    ymax -= margin
    if xmin >= xmax or ymin >= ymax:
        raise ValueError("margin leaves no seeding area")

    def axis(lo, hi):
        # whole steps of res that stay inside [lo, hi] (with a small
        # tolerance so an exact multiple keeps its last seed)
        n = int(np.floor((hi - lo) / res + 1e-9)) + 1
        return np.minimum(lo + res * np.arange(n), hi)

    return axis(xmin, xmax), axis(ymin, ymax)


def ftle_field(field: VelocityField,
               t: float,
               dt: float | None = None,
               res: float | None = None,
               delta: float | None = None,
               margin: float = 0.0,
               integrator: str = "rk4",
               backend: str = "auto") -> LCSResult:
    """Compute an FTLE field over a uniform seed grid.

    Parameters
    ----------
    field : VelocityField
        The (steady) velocity field.
    t : float
        Integration time.  Positive: repelling LCS (forward FTLE).
        Negative: attracting LCS (backward FTLE).
    dt : float, optional
        Advection time step; defaults to ``field.cfl_dt(1.0)`` — one
        grid cell per step.  RK4 is accuracy-limited rather than
        stability-limited, and on the double-gyre benchmark this default
        matches a 4x-finer step to ~1e-4 relative.  Halve it once as a
        self-check on new data; if the field doesn't change, it was fine.
    res : float, optional
        Seed-grid spacing; defaults to the velocity grid spacing.
    delta : float, optional
        Auxiliary particle offset; defaults to ``res / 100``.  It should
        be well below ``res`` so the Jacobian is a *local* derivative.
    margin : float
        Inset of the seed grid from the domain edge.
    integrator : {"rk4", "euler"}
    backend : {"auto", "numpy", "numba"}

    Returns
    -------
    LCSResult
        ``values[j, i]`` is the FTLE at seed ``(x[i], y[j])`` in units
        of 1/time.
    """
    if res is None:
        res = max(field.dx, field.dy)
    sx, sy = seed_grid(field, res=res, margin=margin)  # validates res
    if delta is None:
        delta = res / 100.0
    if not delta > 0.0:
        raise ValueError(f"delta must be positive, got {delta}")
    dt = _check_dt(dt, field, safety=1.0)

    mx, my = np.meshgrid(sx, sy)
    cx = mx.ravel()
    cy = my.ravel()
    n = cx.size

    # primary cloud + 4 auxiliary clouds, advected together in one call
    px = np.concatenate([cx, cx, cx, cx - delta, cx + delta])
    py = np.concatenate([cy, cy + delta, cy - delta, cy, cy])
    ex, ey = advect(field, px, py, t, dt=dt,
                    integrator=integrator, backend=backend)

    x_up, y_up = ex[n:2 * n], ey[n:2 * n]
    x_dn, y_dn = ex[2 * n:3 * n], ey[2 * n:3 * n]
    x_lt, y_lt = ex[3 * n:4 * n], ey[3 * n:4 * n]
    x_rt, y_rt = ex[4 * n:5 * n], ey[4 * n:5 * n]

    inv2d = 1.0 / (2.0 * delta)
    j11 = (x_rt - x_lt) * inv2d   # d(x_final) / d(x_initial)
    j12 = (x_up - x_dn) * inv2d   # d(x_final) / d(y_initial)
    j21 = (y_rt - y_lt) * inv2d   # d(y_final) / d(x_initial)
    j22 = (y_up - y_dn) * inv2d   # d(y_final) / d(y_initial)

    # largest eigenvalue of the 2x2 Cauchy-Green tensor C = J^T J,
    # in closed form (C is symmetric positive semi-definite)
    a = j11 * j11 + j21 * j21
    b = j11 * j12 + j21 * j22
    d = j12 * j12 + j22 * j22
    lam_max = 0.5 * ((a + d) + np.sqrt((a - d) ** 2 + 4.0 * b * b))

    sigma = np.zeros(n)
    ok = lam_max > 0.0
    sigma[ok] = np.log(lam_max[ok]) / (2.0 * abs(t))

    params = {"t": t, "dt": dt,
              "res": res, "delta": delta, "margin": margin,
              "integrator": integrator,
              "direction": "forward (repelling LCS)" if t > 0
                           else "backward (attracting LCS)"}
    return LCSResult(kind="ftle", x=sx, y=sy,
                     values=sigma.reshape(sy.size, sx.size), params=params)
