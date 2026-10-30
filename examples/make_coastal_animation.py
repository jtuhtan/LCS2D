"""Render the LCS of a coastal residual circulation forming, as a GIF.

The flow is :func:`lcs2d.coastal_bay`: a weak alongshore current past a
cape and a bay, with a pair of headland eddies and a bay gyre. Each
frame shows the forward FTLE (repelling LCS, magenta) and the backward
FTLE (attracting LCS, light blue) for a longer integration time, so the
fences between the eddies and the coastal current draw themselves.

Rather than one ``ftle_field`` call per frame, the four particle clouds
of the FTLE stencil are advected ONCE per direction with their positions
recorded, and every frame's FTLE is taken from a recorded snapshot. The
seeds cover only the 15 x 7 km window shown; the field around it keeps
the escaping particles moving. With Numba it runs in under two minutes.

Run:  python examples/make_coastal_animation.py [output.gif] [figure.png]
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

OUT_GIF = sys.argv[1] if len(sys.argv) > 1 else "docs/img/coastal_lcs.gif"
OUT_PNG = sys.argv[2] if len(sys.argv) > 2 else "docs/img/coastal_bay.png"
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


def ftle_snapshots(field, sx, sy, t_max):
    """FTLE of every seed at N_FRAMES times in (0, t_max], from one run."""
    delta = RES / 100.0
    mx, my = np.meshgrid(sx, sy)
    cx, cy = mx.ravel(), my.ravel()
    n = cx.size
    # the same stencil as lcs2d.ftle_field: up, down, left, right
    px = np.concatenate([cx, cx, cx - delta, cx + delta])
    py = np.concatenate([cy + delta, cy - delta, cy, cy])
    dt = field.cfl_dt(1.0)
    n_steps = int(np.ceil(abs(t_max) / dt))
    traj = lcs2d.advect(field, px, py, t_max, dt=abs(t_max) / n_steps,
                        record=True,
                        save_every=max(1, n_steps // N_FRAMES))
    frames = []
    for k in np.linspace(0, traj.t.size - 1, N_FRAMES + 1)[1:].round():
        k = int(k)
        ex, ey = traj.x[k], traj.y[k]
        inv2d = 1.0 / (2.0 * delta)
        j11 = (ex[3 * n:] - ex[2 * n:3 * n]) * inv2d
        j12 = (ex[:n] - ex[n:2 * n]) * inv2d
        j21 = (ey[3 * n:] - ey[2 * n:3 * n]) * inv2d
        j22 = (ey[:n] - ey[n:2 * n]) * inv2d
        a = j11 * j11 + j21 * j21
        b = j11 * j12 + j21 * j22
        d = j12 * j12 + j22 * j22
        lam = 0.5 * ((a + d) + np.sqrt((a - d) ** 2 + 4.0 * b * b))
        sigma = np.log(np.maximum(lam, 1e-300)) / (2.0 * abs(traj.t[k]))
        frames.append((abs(float(traj.t[k])), sigma.reshape(sy.size, sx.size)))
    return frames


def normalized(values, sea):
    # FTLE amplitude decays ~1/t while ridges sharpen: stretch each
    # frame between its own median and 99.5th percentile over the sea
    # (the double-gyre animation does the same with the 98th) to show
    # the structure, not the amplitude
    lo, hi = np.percentile(values[sea], [50.0, 99.5])
    return np.clip((values - lo) / (hi - lo), 0.0, 1.0) if hi > lo         else values * 0.0


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


def title(ax, text):
    ax.set_title(text, loc="left", color=INK, fontsize=10,
                 fontweight="bold")


def main() -> None:
    plt.rcParams["font.family"] = "Arial"
    field = lcs2d.coastal_bay()
    x0, x1, y0, y1 = WINDOW
    sx = np.arange(x0, x1 + RES / 2, RES)
    sy = np.arange(y0, y1 + RES / 2, RES)
    sea = np.meshgrid(sx, sy)[1] > _coastline(np.meshgrid(sx, sy)[0])

    print(f"advecting {4 * sx.size * sy.size:,} particles forward "
          f"and backward for {T_MAX:g} h...", flush=True)
    fwd = ftle_snapshots(field, sx, sy, T_MAX)
    bwd = ftle_snapshots(field, sx, sy, -T_MAX)

    extent = (sx[0] - RES / 2, sx[-1] + RES / 2,
              sy[0] - RES / 2, sy[-1] + RES / 2)
    fig, ax = plt.subplots(figsize=(8.0, 4.3))
    styled_axes(fig, ax)
    im = ax.imshow(blend(sea * 0.0, sea * 0.0), origin="lower",
                   extent=extent, interpolation="bilinear", zorder=1)
    head = ax.set_title("", loc="left", color=INK, fontsize=10,
                        fontweight="bold")
    ax.text(1.0, 1.02, "repelling", color=MAGENTA, fontsize=9,
            fontweight="bold", ha="right", va="bottom",
            transform=ax.transAxes)
    ax.text(0.86, 1.02, "attracting", color=LIGHT_BLUE, fontsize=9,
            fontweight="bold", ha="right", va="bottom",
            transform=ax.transAxes)
    fig.tight_layout()

    frames = []
    for (t, f_vals), (_, b_vals) in zip(fwd, bwd):
        im.set_data(blend(normalized(f_vals, sea) * sea,
                          normalized(b_vals, sea) * sea))
        head.set_text(f"COASTAL LCS FORMING  ·  INTEGRATION TIME {t:4.1f} h")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=100, facecolor=PAGE)
        buf.seek(0)
        frames.append(Image.open(buf).convert("RGB"))

    # quantize every frame with the final frame's palette to avoid
    # palette flicker between frames
    ref = frames[-1].quantize(colors=256, method=Image.MEDIANCUT)
    frames = [f.quantize(palette=ref) for f in frames]
    durations = [FRAME_MS] * (len(frames) - 1) + [FINAL_HOLD_MS]
    frames[0].save(OUT_GIF, save_all=True, append_images=frames[1:],
                   duration=durations, loop=0, optimize=True)
    plt.close(fig)
    print(f"written: {OUT_GIF}  ({len(frames)} frames, "
          f"{os.path.getsize(OUT_GIF) / 1e6:.1f} MB)")

    # static figure: the flow (top) and its finished LCS (bottom)
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8.0, 8.2))
    for ax in (top, bottom):
        styled_axes(fig, ax)
    cols = (field.x >= x0 - 0.5) & (field.x <= x1 + 0.5)
    rows = field.y <= y1 + 0.5
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
    t_end = fwd[-1][0]
    bottom.imshow(blend(normalized(fwd[-1][1], sea) * sea,
                        normalized(bwd[-1][1], sea) * sea),
                  origin="lower", extent=extent, interpolation="bilinear",
                  zorder=1)
    title(bottom, f"b  FORWARD AND BACKWARD FTLE  ·  t = {t_end:g} h")
    bottom.text(1.0, 1.02, "repelling", color=MAGENTA, fontsize=9,
                fontweight="bold", ha="right", va="bottom",
                transform=bottom.transAxes)
    bottom.text(0.84, 1.02, "attracting", color=LIGHT_BLUE, fontsize=9,
                fontweight="bold", ha="right", va="bottom",
                transform=bottom.transAxes)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=110, facecolor=PAGE)
    plt.close(fig)
    print(f"written: {OUT_PNG}")


if __name__ == "__main__":
    main()
