"""Render an animated GIF of the FTLE field forming on the double gyre.

Each frame is a complete FTLE computation with a longer integration time
t, so the repelling ridge sharpens out of the initially featureless
field — the modern successor of the legacy ``ftleplot.m`` video. With
the Numba backend the whole animation computes in about a minute.

Run:  python examples/make_ftle_animation.py [output.gif]
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

OUT = sys.argv[1] if len(sys.argv) > 1 else "docs/img/ftle_formation.gif"
RES = 0.005            # seed spacing (401 x 201 seeds)
T_VALUES = np.linspace(0.25, 12.0, 48)
FRAME_MS = 100         # 10 fps
FINAL_HOLD_MS = 2000   # linger on the finished field


def main() -> None:
    field = lcs2d.double_gyre(nx=401, ny=201)

    print(f"computing {T_VALUES.size} FTLE fields "
          f"({2/RES+1:.0f} x {1/RES+1:.0f} seeds each)...")
    results = []
    for k, t in enumerate(T_VALUES, 1):
        res = lcs2d.ftle_field(field, t=float(t), res=RES)
        results.append(res)
        print(f"  frame {k:2d}/{T_VALUES.size}  t={t:5.2f}  "
              f"max FTLE {res.values.max():.3f}", flush=True)

    # FTLE amplitude decays ~1/t while the ridge sharpens, so each frame
    # is normalized to its own 98th percentile (as the legacy ftleplot.m
    # did): the animation shows the *structure* emerging, not the overall
    # amplitude decay. The colorbar is dropped since units are relative.
    def normalized(values):
        scale = np.percentile(values, 98.0)
        return np.clip(values / scale, 0.0, 1.0) if scale > 0 else values

    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    im = ax.imshow(normalized(results[0].values), origin="lower",
                   extent=(0.0, 2.0, 0.0, 1.0), aspect="equal",
                   cmap="magma", vmin=0.0, vmax=1.0,
                   interpolation="bilinear")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    title = ax.set_title("")
    fig.tight_layout()

    frames = []
    for res in results:
        im.set_data(normalized(res.values))
        title.set_text("Repelling LCS forming — forward FTLE, "
                       f"integration time t = {res.params['t']:5.2f}")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=100)
        buf.seek(0)
        frames.append(Image.open(buf).convert("RGB"))
    plt.close(fig)

    # quantize every frame with the final frame's palette to avoid
    # palette flicker between frames
    ref = frames[-1].quantize(colors=256, method=Image.MEDIANCUT)
    frames = [f.quantize(palette=ref) for f in frames]
    durations = [FRAME_MS] * (len(frames) - 1) + [FINAL_HOLD_MS]
    frames[0].save(OUT, save_all=True, append_images=frames[1:],
                   duration=durations, loop=0, optimize=True)
    size_mb = os.path.getsize(OUT) / 1e6
    print(f"written: {OUT}  ({len(frames)} frames, {size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
