"""Benchmark the lcs2d engine against the legacy algorithm.

Workload: the FTLE computation from the README quick start — double
gyre on a 401 x 201 grid, seeds every 0.02, t = 10 (5 particle clouds
advected simultaneously).

Three contenders:

* ``legacy-style``: what ``advect.m`` did, reproduced in Python — a
  Delaunay-based scattered interpolant (`scipy` `LinearNDInterpolator`,
  the equivalent of MATLAB's ``TriScatteredInterp``) evaluated at every
  particle every time step, with forward-Euler stepping.  Run on a
  subsample and extrapolated, because the full workload would take
  hours.
* ``numpy``: the vectorized bilinear-sampling RK4 engine.
* ``numba``: the same math compiled per particle (if numba is
  installed).

Run:  python benchmarks/benchmark.py
"""

import time

import numpy as np

import lcs2d
from lcs2d import NUMBA_AVAILABLE


def main() -> None:
    field = lcs2d.double_gyre(nx=401, ny=201)
    seed_res = 0.02
    t_total = 10.0
    dt = field.cfl_dt()
    n_steps = max(1, int(round(t_total / dt)))

    sx, sy = lcs2d.seed_grid(field, res=seed_res)
    n_seeds = sx.size * sy.size
    n_particles = 5 * n_seeds  # primary + 4 auxiliary clouds
    total_ps = n_particles * n_steps

    print("workload: FTLE of the double gyre "
          f"({field.shape[1]}x{field.shape[0]} grid)")
    print(f"  {n_seeds} seeds -> {n_particles} particles, "
          f"{n_steps} time steps of dt={dt:.5f} "
          f"({total_ps / 1e6:.1f} M particle-steps)")
    print()

    rows = []

    # -- legacy-style: scattered interpolant per step + Euler -----------
    from scipy.interpolate import LinearNDInterpolator

    mx, my = np.meshgrid(field.x, field.y)
    pts = np.column_stack([mx.ravel(), my.ravel()])
    tic = time.perf_counter()
    interp_u = LinearNDInterpolator(pts, field.u.ravel())
    interp_v = LinearNDInterpolator(pts, field.v.ravel())
    build_s = time.perf_counter() - tic

    n_sub, steps_sub = 2000, 10
    rng = np.random.default_rng(0)
    px = rng.uniform(0.1, 1.9, n_sub)
    py = rng.uniform(0.1, 0.9, n_sub)
    tic = time.perf_counter()
    for _ in range(steps_sub):
        pu = interp_u(px, py)
        pv = interp_v(px, py)
        px = px + dt * np.nan_to_num(pu)
        py = py + dt * np.nan_to_num(pv)
    sub_s = time.perf_counter() - tic
    legacy_s = build_s + sub_s * total_ps / (n_sub * steps_sub)
    rows.append(("legacy-style (scattered interp + Euler, extrapolated)",
                 legacy_s))

    # -- numpy ----------------------------------------------------------
    # dt is pinned so every contender runs the identical 1257-step
    # workload; the actual FTLE default (CFL 1.0) uses half the steps.
    tic = time.perf_counter()
    lcs2d.ftle_field(field, t=t_total, dt=dt, res=seed_res, backend="numpy")
    rows.append(("lcs2d numpy (bilinear + RK4, vectorized)",
                 time.perf_counter() - tic))

    # -- numba ----------------------------------------------------------
    if NUMBA_AVAILABLE:
        # warm up the JIT on a tiny problem so compilation isn't timed
        lcs2d.ftle_field(field, t=2 * dt, dt=dt, res=0.5, backend="numba")
        tic = time.perf_counter()
        lcs2d.ftle_field(field, t=t_total, dt=dt, res=seed_res,
                         backend="numba")
        rows.append(("lcs2d numba (compiled, parallel)",
                     time.perf_counter() - tic))
    else:
        print("numba not installed -- skipping the compiled backend")

    print(f"{'implementation':58s} {'time':>10s} {'speedup':>9s}")
    base = rows[0][1]
    for name, seconds in rows:
        print(f"{name:58s} {seconds:9.2f}s {base / seconds:8.1f}x")
    print()
    print("note: RK4 does 4 velocity samples per step vs Euler's 1, so "
          "the per-sample speedup is ~4x larger still.")


if __name__ == "__main__":
    main()
