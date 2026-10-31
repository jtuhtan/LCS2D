"""Render the LCS of a coastal residual circulation forming, as GIFs.

The flow is :func:`lcs2d.coastal_bay`: a weak alongshore current past a
cape and a bay, with a pair of headland eddies and a bay gyre. Each
frame shows the forward FTLE (repelling LCS, magenta) and the backward
FTLE (attracting LCS, light blue) for a longer integration time, so the
fences between the eddies and the coastal current draw themselves.

Three outputs, in ``docs/img`` by default:

* ``coastal_lcs.gif``: the LCS forming (the README hero);
* ``coastal_bay.gif``: the same under a panel of dye stripes carried by
  the flow for the same time, so the dye can be seen folding along the
  ridges and wrapping up inside the eddies;
* ``coastal_bay.png``: a still of the flow and the finished LCS.

Rather than one ``ftle_field`` call per frame, the four particle clouds
of the FTLE stencil are advected ONCE per direction with their positions
recorded, and every frame's FTLE is taken from a recorded snapshot. The
seeds cover only the 15 x 7 km window shown; the field around it keeps
the escaping particles moving. With Numba it runs in about a minute.

Run:  python examples/make_coastal_animation.py [output_dir]
"""

import io
import os
import sys

import numpy as np
from PIL import Image

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import lcs2d  # noqa: E402
from lcs2d.fields import _coastline  # noqa: E402

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "docs/img"
WINDOW = (0.0, 15.0, 0.0, 7.0)  # x0, x1, y0, y1 of the view [km]
RES = 0.04             # seed spacing [km] (376 x 176 seeds)
T_MAX = 30.0           # longest integration time [h]
N_FRAMES = 48
FRAME_MS = 100         # 10 fps
FINAL_HOLD_MS = 2500   # linger on the finished field

# TalTech night palette
PAGE = np.array([0x12, 0x0F, 0x26]) / 255
MAGENTA = np.array([0xE4, 0x06, 0x7E]) / 255
LIGHT_BLUE = np.array([0x4D, 0xBE, 0xD2]) / 255
LAND = "#252049"
SHORE = "#9396B0"
INK = "#f2f1f8"
INK_2 = "#b6b4cd"
DYE = "#F2C14E"        # tracer stripes: night yellow


def record(field, px, py, t_max):
    """Particle positions at the N_FRAMES times k * t_max / N_FRAMES."""
    per_frame = int(np.ceil(abs(t_max) / field.cfl_dt(1.0) / N_FRAMES))
    n_steps = per_frame * N_FRAMES
    traj = lcs2d.advect(field, px, py, t_max, dt=abs(t_max) / n_steps,
                        record=True, save_every=per_frame)
    return np.abs(traj.t[1:]), traj.x[1:], traj.y[1:]


def ftle_snapshots(field, sx, sy, t_max):
    """FTLE of every seed at the N_FRAMES times, from one run."""
    delta = RES / 100.0
    mx, my = np.meshgrid(sx, sy)
    cx, cy = mx.ravel(), my.ravel()
    n = cx.size
    # the same stencil as lcs2d.ftle_field: up, down, left, right
    px = np.concatenate([cx, cx, cx - delta, cx + delta])
    py = np.concatenate([cy + delta, cy - delta, cy, cy])
    times, tx, ty = record(field, px, py, t_max)
    inv2d = 1.0 / (2.0 * delta)
    frames = []
    for t, ex, ey in zip(times, tx, ty):
        j11 = (ex[3 * n:] - ex[2 * n:3 * n]) * inv2d
        j12 = (ex[:n] - ex[n:2 * n]) * inv2d
        j21 = (ey[3 * n:] - ey[2 * n:3 * n]) * inv2d
        j22 = (ey[:n] - ey[n:2 * n]) * inv2d
        a = j11 * j11 + j21 * j21
        b = j11 * j12 + j21 * j22
        d = j12 * j12 + j22 * j22
        lam = 0.5 * ((a + d) + np.sqrt((a - d) ** 2 + 4.0 * b * b))
        sigma = np.log(np.maximum(lam, 1e-300)) / (2.0 * t)
        frames.append((float(t), sigma.reshape(sy.size, sx.size)))
    return frames


def normalized(values, sea):
    # FTLE amplitude decays ~1/t while ridges sharpen: stretch each
    # frame between its own median and 99.5th percentile over the sea
    # (the double-gyre animation does the same with the 98th) to show
    # the structure, not the amplitude
    lo, hi = np.percentile(values[sea], [50.0, 99.5])
    if hi <= lo:
        return values * 0.0
    return np.clip((values - lo) / (hi - lo), 0.0, 1.0)


def blend(fwd, bwd):
    """Magenta (repelling) and light blue (attracting) on the night page."""
    rgb = (PAGE + fwd[..., None] * (MAGENTA - PAGE)
           + bwd[..., None] * (LIGHT_BLUE - PAGE))
    return np.clip(rgb, 0.0, 1.0)


def styled_axes(fig, ax):
    fig.patch.set_facecolor(PAGE)
    ax.set_facecolor(PAGE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.set_xlabel("alongshore distance (km)", color=INK_2, fontsize=9)
    ax.set_ylabel("offshore (km)", color=INK_2, fontsize=9)
    ax.set_xlim(WINDOW[0], WINDOW[1])
    ax.set_ylim(WINDOW[2], WINDOW[3])
    ax.set_aspect("equal")
    xs = np.linspace(WINDOW[0], WINDOW[1], 600)
    ax.fill_between(xs, WINDOW[2], _coastline(xs), color=LAND, zorder=3)
    ax.plot(xs, _coastline(xs), color=SHORE, lw=1.0, zorder=4)


def title(ax, text=""):
    return ax.set_title(text, loc="left", color=INK, fontsize=10,
                        fontweight="bold")


def legend(ax, at=0.86):
    for x, text, color in ((1.0, "repelling", MAGENTA),
                           (at, "attracting", LIGHT_BLUE)):
        ax.text(x, 1.02, text, color=color, fontsize=9, fontweight="bold",
                ha="right", va="bottom", transform=ax.transAxes)


def grab(fig, dpi):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor=PAGE)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def save_gif(frames, path):
    # quantize every frame with the final frame's palette to avoid
    # palette flicker between frames
    ref = frames[-1].quantize(colors=256, method=Image.MEDIANCUT)
    frames = [f.quantize(palette=ref) for f in frames]
    durations = [FRAME_MS] * (len(frames) - 1) + [FINAL_HOLD_MS]
    frames[0].save(path, save_all=True, append_images=frames[1:],
                   duration=durations, loop=0, optimize=True)
    print(f"written: {path}  ({len(frames)} frames, "
          f"{os.path.getsize(path) / 1e6:.1f} MB)")


def main() -> None:
    plt.rcParams["font.family"] = "Arial"
    os.makedirs(OUT_DIR, exist_ok=True)
    field = lcs2d.coastal_bay()
    x0, x1, y0, y1 = WINDOW
    sx = np.arange(x0, x1 + RES / 2, RES)
    sy = np.arange(y0, y1 + RES / 2, RES)
    mx, my = np.meshgrid(sx, sy)
    sea = my > _coastline(mx)

    print(f"advecting {4 * sx.size * sy.size:,} particles forward "
          f"and backward for {T_MAX:g} h...", flush=True)
    fwd = ftle_snapshots(field, sx, sy, T_MAX)
    bwd = ftle_snapshots(field, sx, sy, -T_MAX)
    times = [t for t, _ in fwd]
    lcs = [blend(normalized(f, sea) * sea, normalized(b, sea) * sea)
           for (_, f), (_, b) in zip(fwd, bwd)]
    extent = (sx[0] - RES / 2, sx[-1] + RES / 2,
              sy[0] - RES / 2, sy[-1] + RES / 2)

    # dye: 0.5 km wide stripes every km, seeded from the upstream edge
    # of the field so dyed water keeps flowing into the view
    gx, gy = np.meshgrid(np.arange(field.x[0], x1, 0.04),
                         np.arange(y0, y1, 0.04))
    dyed = (gy > _coastline(gx) + 0.02) & (np.floor(2.0 * gx) % 2 == 0)
    gx, gy = gx[dyed], gy[dyed]
    _, tracer_x, tracer_y = record(field, gx, gy, T_MAX)

    # 1) the LCS alone
    fig, ax = plt.subplots(figsize=(8.0, 4.3))
    styled_axes(fig, ax)
    im = ax.imshow(lcs[0], origin="lower", extent=extent,
                   interpolation="bilinear", zorder=1)
    head = title(ax, "COASTAL")
    legend(ax)
    fig.tight_layout()
    frames = []
    for t, rgb in zip(times, lcs):
        im.set_data(rgb)
        head.set_text(f"COASTAL LCS FORMING  ·  INTEGRATION TIME {t:4.1f} h")
        frames.append(grab(fig, 100))
    plt.close(fig)
    save_gif(frames, os.path.join(OUT_DIR, "coastal_lcs.gif"))

    # 2) tracers carried by the flow (top) over the LCS forming (bottom)
    cols = (field.x >= x0 - 0.5) & (field.x <= x1 + 0.5)
    rows = field.y <= y1 + 0.5
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8.0, 8.0))
    for ax in (top, bottom):
        styled_axes(fig, ax)
    dots = top.scatter(gx, gy, s=0.5, c=DYE, linewidths=0, zorder=2)
    head_top = title(top, "a")
    im = bottom.imshow(lcs[0], origin="lower", extent=extent,
                       interpolation="bilinear", zorder=1)
    head_bottom = title(bottom, "b")
    legend(bottom)
    fig.tight_layout()
    frames = []
    for k, t in enumerate(times):
        dots.set_offsets(np.column_stack([tracer_x[k], tracer_y[k]]))
        im.set_data(lcs[k])
        head_top.set_text(f"a  DYE STRIPES CARRIED BY THE FLOW  ·  {t:4.1f} h")
        head_bottom.set_text(f"b  FORWARD AND BACKWARD FTLE  ·  "
                             f"t = {t:4.1f} h")
        frames.append(grab(fig, 90))
    plt.close(fig)
    save_gif(frames, os.path.join(OUT_DIR, "coastal_bay.gif"))

    # 3) still: the flow (top) and its finished LCS (bottom)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8.0, 8.2))
    for ax in (top, bottom):
        styled_axes(fig, ax)
    speed = field.speed()[np.ix_(rows, cols)] * (1 / 3.6)  # m/s
    stream = top.streamplot(field.x[cols], field.y[rows],
                            field.u[np.ix_(rows, cols)],
                            field.v[np.ix_(rows, cols)],
                            color=speed, cmap="viridis", density=2.4,
                            linewidth=0.7, arrowsize=0.7, zorder=2)
    bar = fig.colorbar(stream.lines, ax=top, shrink=0.8, pad=0.02)
    bar.set_label("speed (m/s)", color=INK_2, fontsize=9)
    bar.ax.tick_params(colors=INK_2, labelsize=8)
    bar.outline.set_edgecolor(INK_2)
    title(top, "a  RESIDUAL FLOW: CAPE, BAY AND THREE EDDIES")
    bottom.imshow(lcs[-1], origin="lower", extent=extent,
                  interpolation="bilinear", zorder=1)
    title(bottom, f"b  FORWARD AND BACKWARD FTLE  ·  t = {times[-1]:g} h")
    legend(bottom, at=0.84)
    fig.tight_layout()
    png = os.path.join(OUT_DIR, "coastal_bay.png")
    fig.savefig(png, dpi=110, facecolor=PAGE)
    plt.close(fig)
    print(f"written: {png}")


if __name__ == "__main__":
    main()
