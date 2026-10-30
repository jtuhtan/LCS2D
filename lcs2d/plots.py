"""Matplotlib figures saved straight to PNG (no display needed).

Figures are built with :class:`matplotlib.figure.Figure` directly rather
than through ``pyplot``, so they render without a display and importing
this module never changes the caller's matplotlib backend.
"""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure

from .advect import Trajectories
from .fields import VelocityField
from .result import LCSResult

__all__ = ["plot_velocity", "plot_result", "plot_trajectories"]

_CBAR_LABEL = {"ftle": "FTLE  [1/time]", "fsle": "FSLE  [1/time]"}
_DEFAULT_CMAP = {"ftle": "magma", "fsle": "viridis"}


def _new_figure():
    fig = Figure(figsize=(9, 5))
    return fig, fig.add_subplot()


def _finish(fig, ax, path: str, dpi: int) -> None:
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)


def _color_limits(values: np.ndarray) -> tuple[float, float]:
    """Robust (2nd-98th percentile) color limits that never collapse."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0, 1.0
    vmin, vmax = (float(p) for p in np.percentile(finite, [2.0, 98.0]))
    if vmin >= vmax:
        vmin, vmax = float(finite.min()), float(finite.max())
    if vmin >= vmax:  # constant field
        pad = abs(vmin) * 0.5 or 1.0
        vmin, vmax = vmin - pad, vmax + pad
    return vmin, vmax


def plot_velocity(field: VelocityField, path: str, dpi: int = 150) -> None:
    """Speed map with a decimated quiver overlay."""
    fig, ax = _new_figure()
    pcm = ax.pcolormesh(field.x, field.y, field.speed(),
                        cmap="cividis", shading="auto")
    fig.colorbar(pcm, ax=ax, label="speed  [length/time]")
    # at most ~30 arrows along the longer axis
    step = max(1, max(field.shape) // 30)
    mx, my = np.meshgrid(field.x[::step], field.y[::step])
    ax.quiver(mx, my, field.u[::step, ::step], field.v[::step, ::step],
              color="white", width=0.002, alpha=0.8)
    ax.set_title(field.name)
    _finish(fig, ax, path, dpi)


def plot_result(result: LCSResult, path: str,
                cmap: str | None = None, dpi: int = 150) -> None:
    """Filled map of an FTLE/FSLE field with robust color limits."""
    fig, ax = _new_figure()
    vmin, vmax = _color_limits(result.values)
    pcm = ax.pcolormesh(result.x, result.y, result.values,
                        cmap=cmap or _DEFAULT_CMAP.get(result.kind, "magma"),
                        shading="auto", vmin=vmin, vmax=vmax)
    fig.colorbar(pcm, ax=ax,
                 label=_CBAR_LABEL.get(result.kind, result.kind))
    direction = result.params.get("direction", "")
    title = result.kind.upper()
    if direction:
        title += f" — {direction}"
    ax.set_title(title)
    _finish(fig, ax, path, dpi)


def plot_trajectories(traj: Trajectories, path: str,
                      field: VelocityField | None = None,
                      max_lines: int = 2000, dpi: int = 150) -> None:
    """Particle paths (thin lines) with start points marked."""
    fig, ax = _new_figure()
    if field is not None:
        ax.pcolormesh(field.x, field.y, field.speed(),
                      cmap="Greys", shading="auto", alpha=0.6)
    n = traj.x.shape[1]
    sel = np.arange(n) if n <= max_lines else \
        np.linspace(0, n - 1, max_lines).astype(int)
    ax.plot(traj.x[:, sel], traj.y[:, sel], lw=0.4, alpha=0.5)
    ax.plot(traj.x[0, sel], traj.y[0, sel], ".k", ms=2, label="start")
    ax.plot(traj.x[-1, sel], traj.y[-1, sel], ".r", ms=2, label="end")
    ax.legend(loc="upper right")
    ax.set_title(f"particle advection, t = {traj.t[-1]:.4g}")
    _finish(fig, ax, path, dpi)
