"""A ferry crosses the coastal bay: how a passing ship reshapes the LCS.

Builds on ``make_coastal_animation.py`` (run from the repository root):
the first part of the animation is the steady baseline of
:func:`lcs2d.coastal_bay` forming as the integration time grows; once it
has settled, a large ferry enters from the east and crosses the view to
the west along a route 4.6 km offshore, over the bay gyre's saddle, and
the animation slows to half speed. The top panel shows the streamlines of
the instantaneous flow with particles drifting along them, the bottom
panel the forward and backward FTLE.

The ferry adds an unsteady flow to the steady residual circulation:

* **displacement flow**: the hull as a moving 2D dipole (the potential
  flow of a body of radius ``HULL_RADIUS`` pushing water aside), and
* **wake**: a self-propelled ship's wake carries no net momentum (thrust
  balances drag), so it is a propeller jet running aft along the track
  flanked by water dragged forward by the hull,
  ``u ~ (1 - 2 eta^2) exp(-eta^2)`` across the track, which widens
  (``width ~ sqrt(age)``) and decays (e-folding ``WAKE_DECAY``) behind
  the ship.

Both come from one streamfunction, so the added flow is divergence-free;
it is faded to zero onto land. Kelvin waves are left out: their orbital
motion carries almost no depth-averaged transport.

The flow is now time-dependent, so the library's steady engine does not
apply: this script carries its own Numba RK4 kernel with a time-varying
step (fine while the ship is in the domain), and computes FTLE over
sliding windows: forward over [t0, t0 + T], backward over [t0 - T, t0].
Forward FTLE looks ahead, so the ferry's track shows in magenta a little
before the ferry gets there; the backward (attracting) ridges only
appear behind it. Requires Numba (``pip install -e .[fast]``); runs in
about ten minutes.

Run:  python examples/make_ferry_animation.py [output.gif]
"""

import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402
from PIL import Image  # noqa: E402

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover
    sys.exit("this example needs Numba: pip install -e .[fast]")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_coastal_animation as base  # noqa: E402
import lcs2d  # noqa: E402
from lcs2d.fields import _coastline  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else "docs/img/coastal_ferry.gif"

# the ferry (km, h): ~200 m long, 15 knots
SHIP_SPEED = 28.0      # km/h, westward
ROUTE_Y = 4.6          # km offshore
SHIP_X0 = 16.0         # enters the 15 km wide view at clock time 0
SHIP_FROM = 18.0       # the ship (and its wake) start here...
SHIP_TO = -3.0         # ...and the ship is gone past here
HULL_RADIUS = 0.08     # km, dipole radius
WAKE_SPEED = 0.08      # wake jet speed at the stern / ship speed
WAKE_WIDTH = 0.08      # km, wake half-width at the stern
WAKE_GROWTH = 0.1      # h, time for the wake to widen by sqrt(2)
WAKE_DECAY = 0.75      # h, e-folding time of the wake speed

WINDOW_T = 12.0        # FTLE window [h]
RES = base.RES
ICON_SCALE = 2.5       # the ferry icon is drawn this much larger
N_PARTICLES = 1500     # particles drifting along the streamlines
STEADY_STEP = 0.2      # h of drift per frame while the flow is steady
FRAME_MS = 120         # playback; half speed once the ferry is in view


def clock_frames():
    """Clock times t0 [h] of the ferry part: coarse before the ship,
    fine while it crosses (twice the frames, played at half speed),
    medium while the wake decays."""
    return np.concatenate([np.arange(-11.0, -1.0, 2.0),
                           np.arange(-1.0, 0.0, 0.25),
                           np.arange(0.0, 0.6, 0.02),
                           np.arange(0.6, 6.01, 0.3)])


def step_times(t_from, t_to):
    """Step boundaries from t_from to t_to (either direction), small
    while the ship is in the field, moderate while its wake is young."""
    lo, hi = min(t_from, t_to), max(t_from, t_to)
    ship_on = (SHIP_X0 - SHIP_FROM) / SHIP_SPEED - 0.02
    ship_off = (SHIP_X0 - SHIP_TO) / SHIP_SPEED + 0.02
    edges = [lo]
    t = lo
    while t < hi:
        if ship_on <= t < ship_off:
            dt = 0.0015
        elif ship_off <= t < ship_off + 4.0 * WAKE_DECAY:
            dt = 0.01
        else:
            dt = 0.02
        t = min(t + dt, hi)
        edges.append(t)
    edges = np.array(edges)
    return edges if t_to >= t_from else edges[::-1].copy()


@njit(fastmath=True, cache=True)
def _bilinear(g, x0, y0, inv_dx, inv_dy, px, py):
    ny, nx = g.shape
    fx = (px - x0) * inv_dx
    fy = (py - y0) * inv_dy
    if fx < 0.0 or fy < 0.0 or fx > nx - 1 or fy > ny - 1:
        return 0.0
    i = min(int(fx), nx - 2)
    j = min(int(fy), ny - 2)
    ax = fx - i
    ay = fy - j
    return ((1 - ax) * (1 - ay) * g[j, i] + ax * (1 - ay) * g[j, i + 1]
            + (1 - ax) * ay * g[j + 1, i] + ax * ay * g[j + 1, i + 1])


@njit(fastmath=True, cache=True)
def _ship_velocity(x, y, t):
    """Velocity of the ferry's flow, from its streamfunction psi
    (u = -dpsi/dy, v = dpsi/dx) differentiated analytically."""
    xs = SHIP_X0 - SHIP_SPEED * t
    if xs > SHIP_FROM:
        return 0.0, 0.0
    dx = x - xs
    dy = y - ROUTE_Y
    u = 0.0
    v = 0.0
    # moving dipole, psi = U a^2 dy / (r^2 + c): a body of radius a
    # moving with velocity (-U, 0); cut off beyond 3 km, where it is
    # under 6 % of the residual current
    r2 = dx * dx + dy * dy
    if xs > SHIP_TO and r2 < 9.0:
        k = SHIP_SPEED * HULL_RADIUS * HULL_RADIUS
        q = r2 + 0.25 * HULL_RADIUS * HULL_RADIUS
        u -= k * (q - 2.0 * dy * dy) / (q * q)
        v -= 2.0 * k * dx * dy / (q * q)
    # momentumless wake behind the ship, from where it started:
    # psi = -A w eta exp(-eta^2), eta = dy / w, so
    # u = A (1 - 2 eta^2) exp(-eta^2) runs aft on the track
    if dx > 0.0 and x < SHIP_FROM:
        age = dx / SHIP_SPEED
        if age < 6.0 * WAKE_DECAY:
            grow = np.sqrt(1.0 + age / WAKE_GROWTH)
            w = WAKE_WIDTH * grow
            eta = dy / w
            if abs(eta) < 4.0:
                g = np.exp(-eta * eta)
                rise = np.exp(-dx / HULL_RADIUS)
                fade = WAKE_SPEED * SHIP_SPEED * np.exp(-age / WAKE_DECAY)
                amp = fade * (1.0 - rise) / grow
                # d/dx of w and of amp
                dw = WAKE_WIDTH / (2.0 * grow * WAKE_GROWTH * SHIP_SPEED)
                damp = (amp * (-1.0 / (WAKE_DECAY * SHIP_SPEED) - dw / w)
                        + fade * rise / (HULL_RADIUS * grow))
                profile = (1.0 - 2.0 * eta * eta) * g
                u += amp * profile
                # psi_x = -(amp w)' eta g - amp w profile d(eta)/dx,
                # with d(eta)/dx = -eta w' / w
                v -= ((damp * w + amp * dw) * eta * g
                      - amp * profile * eta * dw)
    return u, v


@njit(fastmath=True, cache=True)
def _velocity(u, v, wet, x0, y0, inv_dx, inv_dy, px, py, t):
    ub = _bilinear(u, x0, y0, inv_dx, inv_dy, px, py)
    vb = _bilinear(v, x0, y0, inv_dx, inv_dy, px, py)
    s = _bilinear(wet, x0, y0, inv_dx, inv_dy, px, py)
    if s > 0.0:
        us, vs = _ship_velocity(px, py, t)
        ub += s * us
        vb += s * vs
    return ub, vb


@njit(parallel=True, fastmath=True, cache=True)
def _advance(u, v, wet, x0, y0, inv_dx, inv_dy, px, py, times):
    for p in prange(px.size):
        x = px[p]
        y = py[p]
        for k in range(times.size - 1):
            t = times[k]
            h = times[k + 1] - t
            k1u, k1v = _velocity(u, v, wet, x0, y0, inv_dx, inv_dy, x, y, t)
            k2u, k2v = _velocity(u, v, wet, x0, y0, inv_dx, inv_dy,
                                 x + 0.5 * h * k1u, y + 0.5 * h * k1v,
                                 t + 0.5 * h)
            k3u, k3v = _velocity(u, v, wet, x0, y0, inv_dx, inv_dy,
                                 x + 0.5 * h * k2u, y + 0.5 * h * k2v,
                                 t + 0.5 * h)
            k4u, k4v = _velocity(u, v, wet, x0, y0, inv_dx, inv_dy,
                                 x + h * k3u, y + h * k3v, t + h)
            x += h * (k1u + 2.0 * k2u + 2.0 * k3u + k4u) / 6.0
            y += h * (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0
        px[p] = x
        py[p] = y


@njit(parallel=True, fastmath=True, cache=True)
def _grid_velocity(u, v, wet, x0, y0, inv_dx, inv_dy, gx, gy, t):
    gu = np.empty((gy.size, gx.size))
    gv = np.empty((gy.size, gx.size))
    for j in prange(gy.size):
        for i in range(gx.size):
            gu[j, i], gv[j, i] = _velocity(u, v, wet, x0, y0, inv_dx,
                                           inv_dy, gx[i], gy[j], t)
    return gu, gv


class Flow:
    """The steady coastal field plus the ferry, advected by _advance."""

    def __init__(self, field):
        self.field = field
        mx, my = np.meshgrid(field.x, field.y)
        # 1 at sea, fading to 0 within 300 m of the shore and on land
        wet = np.clip((my - _coastline(mx)) / 0.3, 0.0, 1.0)
        self.args = (field.u, field.v, wet, field.x[0], field.y[0],
                     1.0 / field.dx, 1.0 / field.dy)

    def advect(self, px, py, t_from, t_to):
        px = np.array(px, dtype=np.float64)
        py = np.array(py, dtype=np.float64)
        _advance(*self.args, px, py, step_times(t_from, t_to))
        return px, py

    def velocity_grid(self, gx, gy, t):
        """Instantaneous velocity (u, v) at clock time t on a grid."""
        return _grid_velocity(*self.args, gx, gy, t)


def ftle(flow, sx, sy, t0, t1):
    """FTLE of the seed grid over the flow from t0 to t1."""
    delta = RES / 100.0
    mx, my = np.meshgrid(sx, sy)
    cx, cy = mx.ravel(), my.ravel()
    n = cx.size
    px = np.concatenate([cx, cx, cx - delta, cx + delta])
    py = np.concatenate([cy + delta, cy - delta, cy, cy])
    ex, ey = flow.advect(px, py, t0, t1)
    inv2d = 1.0 / (2.0 * delta)
    j11 = (ex[3 * n:] - ex[2 * n:3 * n]) * inv2d
    j12 = (ex[:n] - ex[n:2 * n]) * inv2d
    j21 = (ey[3 * n:] - ey[2 * n:3 * n]) * inv2d
    j22 = (ey[:n] - ey[n:2 * n]) * inv2d
    a = j11 * j11 + j21 * j21
    b = j11 * j12 + j21 * j22
    d = j12 * j12 + j22 * j22
    lam = 0.5 * ((a + d) + np.sqrt((a - d) ** 2 + 4.0 * b * b))
    sigma = np.log(np.maximum(lam, 1e-300)) / (2.0 * abs(t1 - t0))
    return sigma.reshape(sy.size, sx.size)


def ferry_outline(xs, scale):
    """A ferry hull in plan view, bow to the west, ~200 m x 30 m."""
    length, beam = 0.2 * scale, 0.03 * scale
    shape = np.array([(-0.5, 0.0), (-0.3, 0.5), (0.5, 0.5),
                      (0.5, -0.5), (-0.3, -0.5)])
    return np.column_stack([xs + shape[:, 0] * length,
                            ROUTE_Y + shape[:, 1] * beam])


def main() -> None:
    plt.rcParams["font.family"] = "Arial"
    field = lcs2d.coastal_bay()
    flow = Flow(field)
    x0, x1, y0, y1 = base.WINDOW
    sx = np.arange(x0, x1 + RES / 2, RES)
    sy = np.arange(y0, y1 + RES / 2, RES)
    mx, my = np.meshgrid(sx, sy)
    sea = my > _coastline(mx)
    extent = (sx[0] - RES / 2, sx[-1] + RES / 2,
              sy[0] - RES / 2, sy[-1] + RES / 2)

    t_ferry = clock_frames()
    rng = np.random.default_rng(1)

    def sea_points(n):
        """n random points at sea in the view."""
        px = rng.uniform(x0, x1, 4 * n + 16)
        py = rng.uniform(y0, y1, 4 * n + 16)
        keep = py > _coastline(px) + 0.05
        return px[keep][:n], py[keep][:n]

    def drift(px, py, t_from, t_to):
        """Advance the particles; re-seed those that left the view (they
        get no tail in that frame)."""
        qx, qy = flow.advect(px, py, t_from, t_to)
        gone = ((qx < x0) | (qx > x1) | (qy > y1)
                | (qy < _coastline(qx) + 0.01))
        tail_x, tail_y = px.copy(), py.copy()
        tail_x[gone] = np.nan
        tail_y[gone] = np.nan
        qx[gone], qy[gone] = sea_points(int(gone.sum()))
        return qx, qy, tail_x, tail_y

    # each frame: (title, LCS rgb, particle x, y, tail x, y,
    #              clock time of the streamlines or None, ship x or None)
    frames = []
    px, py = sea_points(N_PARTICLES)
    clock = t_ferry[0] - WINDOW_T  # steady here: no ship for hours

    # part 1: the steady baseline forming, before the ferry is anywhere
    # near: the steady engine's recorded runs from make_coastal_animation
    print("part 1: baseline LCS forming...", flush=True)
    fwd = base.ftle_snapshots(field, sx, sy, WINDOW_T)[1::2]
    bwd = base.ftle_snapshots(field, sx, sy, -WINDOW_T)[1::2]
    for (t, f), (_, b) in zip(fwd, bwd):
        px, py, tx, ty = drift(px, py, clock, clock + STEADY_STEP)
        clock += STEADY_STEP
        frames.append((f"BASELINE  ·  INTEGRATION TIME {t:4.1f} h",
                       base.blend(base.normalized(f, sea) * sea,
                                  base.normalized(b, sea) * sea),
                       px, py, tx, ty, None, None))

    # part 2: sliding windows while the ferry comes and goes; until it is
    # near, the flow is still steady and the particles keep drifting
    # STEADY_STEP per frame, from then on they keep the clock
    print(f"part 2: {t_ferry.size} frames with the ferry...", flush=True)
    t_prev = t_ferry[0]
    for k, t0 in enumerate(t_ferry):
        f = ftle(flow, sx, sy, t0, t0 + WINDOW_T)
        b = ftle(flow, sx, sy, t0, t0 - WINDOW_T)
        if t0 < -0.1:
            px, py, tx, ty = drift(px, py, t0 - STEADY_STEP, t0)
        else:
            px, py, tx, ty = drift(px, py, t_prev, t0)
        t_prev = t0
        xs = SHIP_X0 - SHIP_SPEED * t0
        text = (f"FERRY DUE IN {-t0:4.2f} h" if t0 < 0.0
                else f"FERRY ENTERED {t0:4.2f} h AGO")
        frames.append((text,
                       base.blend(base.normalized(f, sea) * sea,
                                  base.normalized(b, sea) * sea),
                       px, py, tx, ty, t0 if t0 >= 0.0 else None,
                       xs if x0 - 0.5 < xs < x1 + 0.5 else None))
        print(f"  frame {k + 1:2d}/{t_ferry.size}  t0 = {t0:+5.2f} h",
              flush=True)

    # render: streamlines and particles (top), LCS (bottom), the ferry
    # on both
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8.0, 8.0))
    for ax in (top, bottom):
        base.styled_axes(fig, ax)
    gx = np.arange(x0 - 0.5, x1 + 0.5 + 1e-9, 0.025)
    gy = np.arange(y0, y1 + 0.5 + 1e-9, 0.025)
    mgx, mgy = np.meshgrid(gx, gy)
    sea_grid = mgy > _coastline(mgx)

    def streamlines(t):
        # t = None: the steady flow, long before the ferry
        gu, gv = flow.velocity_grid(gx, gy, -50.0 if t is None else t)
        gu = np.ma.masked_where(~sea_grid, gu)
        gv = np.ma.masked_where(~sea_grid, gv)
        return top.streamplot(gx, gy, gu, gv, color=base.SHORE,
                              density=1.8, linewidth=0.45, arrowsize=0.6,
                              zorder=1)

    stream = streamlines(None)
    tails = LineCollection([], colors=base.DYE, linewidths=1.0,
                           alpha=0.55, zorder=2)
    top.add_collection(tails)
    dots = top.scatter(px, py, s=5, c=base.DYE, linewidths=0, zorder=3)
    im = bottom.imshow(frames[0][1], origin="lower", extent=extent,
                       interpolation="bilinear", zorder=1)
    head_top = base.title(top, "a")
    head_bottom = base.title(bottom, "b")
    base.legend(bottom)
    ships = [Polygon(ferry_outline(-10.0, ICON_SCALE), closed=True,
                     facecolor=base.INK, edgecolor=base.PAGE, lw=0.6,
                     zorder=6) for _ in (top, bottom)]
    for ax, ship in zip((top, bottom), ships):
        ax.add_patch(ship)
    fig.tight_layout()

    images, durations = [], []
    for text, rgb, fx, fy, tx, ty, t_stream, xs in frames:
        if t_stream is not None:  # the ferry's flow changes them
            stream.lines.remove()
            for arrow in [a for a in top.patches if a not in ships]:
                arrow.remove()
            stream = streamlines(t_stream)
        dots.set_offsets(np.column_stack([fx, fy]))
        tails.set_segments(np.stack([np.column_stack([tx, ty]),
                                     np.column_stack([fx, fy])], axis=1))
        im.set_data(rgb)
        head_top.set_text(f"a  STREAMLINES  ·  {text}")
        head_bottom.set_text("b  FORWARD AND BACKWARD FTLE, "
                             f"{WINDOW_T:g} h WINDOWS")
        for ship in ships:
            ship.set_xy(ferry_outline(-10.0 if xs is None else xs,
                                      ICON_SCALE))
        images.append(base.grab(fig, 90))
        durations.append(FRAME_MS if t_stream is None else 2 * FRAME_MS)
    plt.close(fig)

    # quantize every frame with the final frame's palette to avoid
    # palette flicker between frames
    ref = images[-1].quantize(colors=256, method=Image.MEDIANCUT)
    images = [frame.quantize(palette=ref) for frame in images]
    durations[-1] = base.FINAL_HOLD_MS
    images[0].save(OUT, save_all=True, append_images=images[1:],
                   duration=durations, loop=0, optimize=True)
    print(f"written: {OUT}  ({len(images)} frames, "
          f"{os.path.getsize(OUT) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
