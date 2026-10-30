"""Finite-size Lyapunov exponent (FSLE) fields.

Method (after Aurell et al. 1997; d'Ovidio et al. 2004): each seed is
surrounded by four auxiliary particles at distance ``d0 = delta``.  The
five points are advected together, and ``tau`` is the first time at
which *any* auxiliary particle is separated from the primary by at
least ``df = alpha * d0``.  Then

    lambda = ln(alpha) / tau .

Seeds whose auxiliary particles never separate within ``t_max`` get
``lambda = 0`` (no measurable stretching at this scale — typically gyre
cores and stagnant regions).

This replaces the legacy ``fsle.m`` / ``fsleMAX.m``, which looped over
every seed and offset direction in interpreted MATLAB to find the first
threshold crossing.  Here the search runs inside the advection loop —
vectorized in NumPy, or fully compiled per particle with Numba — and
stops early once every seed has crossed.
"""

from __future__ import annotations

import numpy as np

from . import _kernels
from .advect import _check_dt, _field_args, _step_numpy, resolve_backend
from .fields import VelocityField
from .ftle import seed_grid
from .result import LCSResult

__all__ = ["fsle_field"]


def _fsle_tau_steps_numpy(field, cx, cy, delta, df, dt, max_steps, use_rk4):
    """NumPy twin of the Numba FSLE kernel: first-crossing step per seed."""
    args = _field_args(field)
    n = cx.size
    # rows: primary, up, down, left, right
    px = np.vstack([cx, cx, cx, cx - delta, cx + delta])
    py = np.vstack([cy, cy + delta, cy - delta, cy, cy])
    tau_steps = np.zeros(n, dtype=np.int64)
    active = np.arange(n)
    df2 = df * df
    for s in range(1, max_steps + 1):
        flat_x, flat_y = _step_numpy(*args, px[:, active].ravel(),
                                     py[:, active].ravel(), dt, use_rk4)
        px[:, active] = flat_x.reshape(5, active.size)
        py[:, active] = flat_y.reshape(5, active.size)
        sep2 = ((px[1:, active] - px[0, active]) ** 2
                + (py[1:, active] - py[0, active]) ** 2)
        crossed = (sep2 >= df2).any(axis=0)
        tau_steps[active[crossed]] = s
        active = active[~crossed]
        if active.size == 0:
            break
    return tau_steps


def fsle_field(field: VelocityField,
               t_max: float,
               dt: float | None = None,
               res: float | None = None,
               delta: float | None = None,
               alpha: float = 2.0,
               margin: float = 0.0,
               integrator: str = "rk4",
               backend: str = "auto") -> LCSResult:
    """Compute an FSLE field over a uniform seed grid.

    Parameters
    ----------
    field : VelocityField
        The (steady) velocity field.
    t_max : float
        Longest time to wait for the separation threshold.  Negative
        values integrate backward in time (attracting structures).
    dt : float, optional
        Advection time step; defaults to ``field.cfl_dt()``.
    res : float, optional
        Seed-grid spacing; defaults to the velocity grid spacing.
    delta : float, optional
        Initial separation ``d0``; defaults to ``res / 2``.  Unlike the
        FTLE offset this sets the *physical scale* being probed, so it
        is usually chosen comparable to the seed spacing.
    alpha : float
        Separation ratio ``df / d0`` (must be > 1).  The legacy scripts
        used 2 and 3.
    margin : float
        Inset of the seed grid from the domain edge.
    integrator : {"rk4", "euler"}
    backend : {"auto", "numpy", "numba"}

    Returns
    -------
    LCSResult
        ``values`` holds lambda in 1/time (0 where the threshold was
        never reached); ``extra["tau"]`` holds the separation time in
        time units (NaN where never reached).
    """
    if not alpha > 1.0:
        raise ValueError("alpha (= df/d0) must be greater than 1")
    if integrator not in ("rk4", "euler"):
        raise ValueError(f"unknown integrator {integrator!r}; "
                         "expected 'rk4' or 'euler'")
    use_rk4 = integrator == "rk4"
    backend = resolve_backend(backend)
    if res is None:
        res = max(field.dx, field.dy)
    sx, sy = seed_grid(field, res=res, margin=margin)  # validates res
    if delta is None:
        delta = res / 2.0
    if not delta > 0.0:
        raise ValueError(f"delta must be positive, got {delta}")
    t_max = float(t_max)
    if not np.isfinite(t_max) or t_max == 0.0:
        raise ValueError("t_max must be finite and nonzero "
                         "(negative t_max integrates backward in time)")
    dt = _check_dt(dt, field)
    max_steps = max(1, int(round(abs(t_max) / dt)))
    sdt = np.copysign(dt, t_max)

    mx, my = np.meshgrid(sx, sy)
    cx = np.ascontiguousarray(mx.ravel())
    cy = np.ascontiguousarray(my.ravel())
    df = alpha * delta

    if backend == "numba":
        tau_steps = _kernels.fsle_kernel(*_field_args(field), cx, cy,
                                         delta, df, sdt, max_steps, use_rk4)
    else:
        tau_steps = _fsle_tau_steps_numpy(field, cx, cy, delta, df, sdt,
                                          max_steps, use_rk4)

    tau = np.full(cx.size, np.nan)
    lam = np.zeros(cx.size)
    crossed = tau_steps > 0
    tau[crossed] = tau_steps[crossed] * dt
    lam[crossed] = np.log(alpha) / tau[crossed]

    params = {"t_max": t_max, "dt": dt, "res": res, "delta": delta,
              "alpha": alpha, "margin": margin, "integrator": integrator,
              "direction": "forward (repelling LCS)" if t_max > 0
                           else "backward (attracting LCS)",
              "n_never_separated": int((~crossed).sum())}
    shape = (sy.size, sx.size)
    return LCSResult(kind="fsle", x=sx, y=sy,
                     values=lam.reshape(shape), params=params,
                     extra={"tau": tau.reshape(shape)})
