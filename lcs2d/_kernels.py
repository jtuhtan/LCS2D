"""Optional Numba-accelerated compute kernels.

When Numba is installed (``pip install lcs2d[fast]``) these kernels are
compiled to parallel machine code and used automatically.  Every kernel
has a pure-NumPy twin (in :mod:`lcs2d.advect` / :mod:`lcs2d.fsle`); the
test suite enforces agreement between the two backends, so Numba is a
speed-up, never a requirement.

The kernels compile with ``fastmath=True`` (fused multiply-add and
instruction reordering, measured ~12% faster).  Positions may therefore
differ from the NumPy backend at the last-few-ulp level (~1e-14 for the
double-gyre benchmark) - orders of magnitude below the RK4 truncation
error, and bounded by the parity tests.

The bilinear sampler treats any point outside the grid as having zero
velocity: particles that leave the domain freeze in place.  NaN input
velocities were already replaced by zero when the field was built, so
masked regions behave the same way.
"""

from __future__ import annotations

import numpy as np

try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised on numba-free installs
    NUMBA_AVAILABLE = False

__all__ = ["NUMBA_AVAILABLE"]

if NUMBA_AVAILABLE:

    @njit(cache=True, inline="always", fastmath=True)
    def _sample_uv_point(u, v, x0, y0, inv_dx, inv_dy, px, py):
        ny, nx = u.shape
        gx = (px - x0) * inv_dx
        gy = (py - y0) * inv_dy
        if gx < 0.0 or gy < 0.0 or gx > nx - 1.0 or gy > ny - 1.0:
            return 0.0, 0.0
        i = int(gx)
        j = int(gy)
        if i > nx - 2:
            i = nx - 2
        if j > ny - 2:
            j = ny - 2
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
        return us, vs

    @njit(cache=True, inline="always", fastmath=True)
    def _step_point(u, v, x0, y0, inv_dx, inv_dy, px, py, dt, use_rk4):
        if use_rk4:
            k1u, k1v = _sample_uv_point(u, v, x0, y0, inv_dx, inv_dy, px, py)
            k2u, k2v = _sample_uv_point(u, v, x0, y0, inv_dx, inv_dy,
                                        px + 0.5 * dt * k1u, py + 0.5 * dt * k1v)
            k3u, k3v = _sample_uv_point(u, v, x0, y0, inv_dx, inv_dy,
                                        px + 0.5 * dt * k2u, py + 0.5 * dt * k2v)
            k4u, k4v = _sample_uv_point(u, v, x0, y0, inv_dx, inv_dy,
                                        px + dt * k3u, py + dt * k3v)
            return (px + dt * (k1u + 2.0 * (k2u + k3u) + k4u) / 6.0,
                    py + dt * (k1v + 2.0 * (k2v + k3v) + k4v) / 6.0)
        su, sv = _sample_uv_point(u, v, x0, y0, inv_dx, inv_dy, px, py)
        return px + dt * su, py + dt * sv

    @njit(cache=True, parallel=True, fastmath=True)
    def advect_final_kernel(u, v, x0, y0, inv_dx, inv_dy,
                            px, py, dt, n_steps, use_rk4):
        n = px.shape[0]
        out_x = np.empty(n, dtype=np.float64)
        out_y = np.empty(n, dtype=np.float64)
        for p in prange(n):
            xx = px[p]
            yy = py[p]
            for _ in range(n_steps):
                xx, yy = _step_point(u, v, x0, y0, inv_dx, inv_dy,
                                     xx, yy, dt, use_rk4)
            out_x[p] = xx
            out_y[p] = yy
        return out_x, out_y

    @njit(cache=True, parallel=True, fastmath=True)
    def advect_record_kernel(u, v, x0, y0, inv_dx, inv_dy,
                             px, py, dt, n_steps, save_every, use_rk4):
        n = px.shape[0]
        n_saved = n_steps // save_every + 1
        if n_steps % save_every != 0:
            n_saved += 1  # final step is always recorded
        tx = np.empty((n_saved, n), dtype=np.float64)
        ty = np.empty((n_saved, n), dtype=np.float64)
        for p in prange(n):
            xx = px[p]
            yy = py[p]
            tx[0, p] = xx
            ty[0, p] = yy
            k = 1
            for s in range(1, n_steps + 1):
                xx, yy = _step_point(u, v, x0, y0, inv_dx, inv_dy,
                                     xx, yy, dt, use_rk4)
                if s % save_every == 0 or s == n_steps:
                    tx[k, p] = xx
                    ty[k, p] = yy
                    k += 1
        return tx, ty

    @njit(cache=True, parallel=True, fastmath=True)
    def fsle_kernel(u, v, x0, y0, inv_dx, inv_dy,
                    px, py, delta, df, dt, max_steps, use_rk4):
        n = px.shape[0]
        tau_steps = np.zeros(n, dtype=np.int64)  # 0 means "never separated"
        df2 = df * df
        for p in prange(n):
            cx = px[p]
            cy = py[p]
            ox = np.empty(4, dtype=np.float64)
            oy = np.empty(4, dtype=np.float64)
            ox[0] = cx
            oy[0] = cy + delta  # up
            ox[1] = cx
            oy[1] = cy - delta  # down
            ox[2] = cx - delta
            oy[2] = cy          # left
            ox[3] = cx + delta
            oy[3] = cy          # right
            for s in range(1, max_steps + 1):
                cx, cy = _step_point(u, v, x0, y0, inv_dx, inv_dy,
                                     cx, cy, dt, use_rk4)
                crossed = False
                for o in range(4):
                    ox[o], oy[o] = _step_point(u, v, x0, y0, inv_dx, inv_dy,
                                               ox[o], oy[o], dt, use_rk4)
                    ddx = ox[o] - cx
                    ddy = oy[o] - cy
                    if ddx * ddx + ddy * ddy >= df2:
                        crossed = True
                if crossed:
                    tau_steps[p] = s
                    break
        return tau_steps
