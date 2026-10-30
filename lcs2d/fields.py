"""Velocity-field container and analytic test fields.

A :class:`VelocityField` holds a steady 2D velocity field sampled on a
regular (uniformly spaced) grid.  All LCS computations in this package
sample the field with bilinear interpolation, so the grid must be
regular; scattered data is first resampled onto a grid (see
:func:`lcs2d.io.load_field` and :meth:`VelocityField.from_scattered`).

Conventions
-----------
* ``x`` has shape ``(nx,)`` and ``y`` has shape ``(ny,)``, both strictly
  increasing and uniformly spaced.
* ``u`` and ``v`` have shape ``(ny, nx)``: row index selects ``y``,
  column index selects ``x`` (the ``numpy.meshgrid`` default).
* Units are whatever your data uses (e.g. metres and metres/second);
  FTLE/FSLE then come out in 1/time of those units.
* ``NaN`` velocity samples are treated as no-flow cells (land, masked
  PIV vectors, obstacles): they are replaced by zero and counted in
  :attr:`VelocityField.n_masked`.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field

import numpy as np

__all__ = ["VelocityField", "double_gyre"]


def _check_uniform(coord: np.ndarray, name: str) -> float:
    """Validate a strictly increasing, uniformly spaced 1D coordinate array."""
    if coord.ndim != 1 or coord.size < 2:
        raise ValueError(f"{name} must be a 1D array with at least 2 entries, "
                         f"got shape {coord.shape}")
    steps = np.diff(coord)
    if np.any(steps <= 0):
        raise ValueError(f"{name} must be strictly increasing")
    mean_step = float(steps.mean())
    if not np.allclose(steps, mean_step, rtol=1e-6, atol=0.0):
        raise ValueError(
            f"{name} must be uniformly spaced (bilinear sampling assumes a "
            f"regular grid); spacing varies between {steps.min():.6g} and "
            f"{steps.max():.6g}. Resample your data first, e.g. with "
            f"VelocityField.from_scattered().")
    return mean_step


@dataclass
class VelocityField:
    """A steady 2D velocity field on a regular grid.

    Parameters
    ----------
    x, y : ndarray
        Grid coordinates, shapes ``(nx,)`` and ``(ny,)``.
    u, v : ndarray
        Velocity components, shape ``(ny, nx)``.
    name : str
        Free-form label used in printouts and plot titles.
    """

    x: np.ndarray
    y: np.ndarray
    u: np.ndarray
    v: np.ndarray
    name: str = "velocity field"
    n_masked: int = dataclass_field(default=0, init=False)

    def __post_init__(self) -> None:
        self.x = np.ascontiguousarray(self.x, dtype=np.float64)
        self.y = np.ascontiguousarray(self.y, dtype=np.float64)
        self.u = np.ascontiguousarray(self.u, dtype=np.float64)
        self.v = np.ascontiguousarray(self.v, dtype=np.float64)
        self._dx = _check_uniform(self.x, "x")
        self._dy = _check_uniform(self.y, "y")
        expected = (self.y.size, self.x.size)
        if self.u.shape != expected or self.v.shape != expected:
            raise ValueError(
                f"u and v must have shape (ny, nx) = {expected}; "
                f"got u {self.u.shape}, v {self.v.shape}")
        bad = ~(np.isfinite(self.u) & np.isfinite(self.v))
        self.n_masked = int(bad.sum())
        if self.n_masked:
            self.u[bad] = 0.0
            self.v[bad] = 0.0

    # -- basic geometry -------------------------------------------------
    @property
    def dx(self) -> float:
        return self._dx

    @property
    def dy(self) -> float:
        return self._dy

    @property
    def extent(self) -> tuple[float, float, float, float]:
        """(xmin, xmax, ymin, ymax) of the grid."""
        return float(self.x[0]), float(self.x[-1]), float(self.y[0]), float(self.y[-1])

    @property
    def shape(self) -> tuple[int, int]:
        """(ny, nx)."""
        return self.u.shape

    def speed(self) -> np.ndarray:
        """Velocity magnitude, shape ``(ny, nx)``."""
        return np.hypot(self.u, self.v)

    def max_speed(self) -> float:
        return float(self.speed().max())

    def cfl_dt(self, safety: float = 0.5) -> float:
        """Suggest an advection time step from the CFL condition.

        Returns ``safety * min(dx, dy) / max_speed`` — the largest step
        for which a particle cannot skip over a grid cell.  ``safety``
        of 0.5 is a comfortable default; use smaller values for strongly
        curved flows.
        """
        if not safety > 0.0:
            raise ValueError(f"safety must be positive, got {safety}")
        vmax = self.max_speed()
        if vmax <= 0.0:
            raise ValueError("field has zero velocity everywhere; "
                             "cannot suggest a time step")
        return safety * min(self._dx, self._dy) / vmax

    # -- constructors ---------------------------------------------------
    @classmethod
    def from_scattered(cls,
                       x: np.ndarray,
                       y: np.ndarray,
                       u: np.ndarray,
                       v: np.ndarray,
                       grid_res: float | None = None,
                       name: str = "resampled field") -> "VelocityField":
        """Build a regular-grid field from scattered ``(x, y, u, v)`` samples.

        This replaces the ``TriScatteredInterp`` workflow of the legacy
        MATLAB code: a Delaunay-based linear interpolation is done ONCE
        onto a regular grid, and the fast bilinear sampler is used from
        then on.

        Parameters
        ----------
        x, y, u, v : ndarray
            1D arrays of equal length: sample positions and velocities.
        grid_res : float, optional
            Grid spacing.  Default: an automatic estimate that gives the
            grid about as many nodes as there are scattered samples,
            ``sqrt(bounding_box_area / n_samples)``.
        """
        from scipy.interpolate import griddata

        x = np.asarray(x, dtype=np.float64).ravel()
        y = np.asarray(y, dtype=np.float64).ravel()
        u = np.asarray(u, dtype=np.float64).ravel()
        v = np.asarray(v, dtype=np.float64).ravel()
        if not (x.size == y.size == u.size == v.size):
            raise ValueError("x, y, u, v must all have the same length")
        if x.size < 4:
            raise ValueError("need at least 4 scattered samples")

        ok = np.isfinite(x) & np.isfinite(y)
        x, y, u, v = x[ok], y[ok], u[ok], v[ok]
        if x.size < 4:
            raise ValueError("need at least 4 scattered samples with "
                             "finite coordinates")

        xmin, xmax = float(x.min()), float(x.max())
        ymin, ymax = float(y.min()), float(y.max())
        if not (xmax > xmin and ymax > ymin):
            raise ValueError(
                "scattered samples span no area (all x or all y equal); "
                "a 2D field needs points spread over a plane")
        if grid_res is None:
            area = (xmax - xmin) * (ymax - ymin)
            grid_res = float(np.sqrt(area / x.size))
        if not (np.isfinite(grid_res) and grid_res > 0.0):
            raise ValueError(f"grid_res must be positive, got {grid_res}")
        nx = max(2, int(round((xmax - xmin) / grid_res)) + 1)
        ny = max(2, int(round((ymax - ymin) / grid_res)) + 1)
        gx = np.linspace(xmin, xmax, nx)
        gy = np.linspace(ymin, ymax, ny)
        mx, my = np.meshgrid(gx, gy)

        pts = np.column_stack([x, y])
        gu = griddata(pts, u, (mx, my), method="linear")
        gv = griddata(pts, v, (mx, my), method="linear")
        # outside the convex hull of the samples -> no-flow (NaN -> 0 in
        # __post_init__, counted in n_masked)
        return cls(gx, gy, gu, gv, name=name)

    # -- persistence ----------------------------------------------------
    def save(self, path: str) -> None:
        """Save to a compressed ``.npz`` with keys ``x, y, u, v, name``."""
        np.savez_compressed(path, x=self.x, y=self.y, u=self.u, v=self.v,
                            name=np.array(self.name))

    @classmethod
    def load(cls, path: str) -> "VelocityField":
        """Load a field saved with :meth:`save` (keys ``x, y, u, v``)."""
        with np.load(path, allow_pickle=False) as data:
            for key in ("x", "y", "u", "v"):
                if key not in data:
                    raise KeyError(
                        f"{path!s} is missing array {key!r}; a field .npz "
                        f"needs x (nx,), y (ny,), u (ny,nx), v (ny,nx)")
            name = str(data["name"]) if "name" in data else "velocity field"
            return cls(data["x"], data["y"], data["u"], data["v"], name=name)

    def summary(self) -> str:
        """A short human-readable description of the field."""
        xmin, xmax, ymin, ymax = self.extent
        ny, nx = self.shape
        lines = [
            f"field   : {self.name}",
            f"grid    : {nx} x {ny} nodes, spacing dx={self.dx:.6g}, dy={self.dy:.6g}",
            f"domain  : x in [{xmin:.6g}, {xmax:.6g}], y in [{ymin:.6g}, {ymax:.6g}]",
            f"speed   : max {self.max_speed():.6g}, mean {self.speed().mean():.6g}",
        ]
        if self.n_masked:
            frac = self.n_masked / self.u.size
            lines.append(f"masked  : {self.n_masked} no-flow nodes "
                         f"({frac:.1%} of grid; NaN input treated as zero velocity)")
        return "\n".join(lines)


def double_gyre(nx: int = 401,
                ny: int = 201,
                amplitude: float = 0.1) -> VelocityField:
    """The classic steady double-gyre test field on [0, 2] x [0, 1].

    Streamfunction ``psi = A sin(pi x) sin(pi y)`` gives two
    counter-rotating gyres separated by the line ``x = 1``.  The forward
    FTLE field has its ridge along that separatrix, which makes this the
    standard smoke test for LCS codes (and matches the legacy
    ``gyrefield.m``, which built the same field via a numerical gradient
    of psi; here the derivatives are analytic).

    u = -d(psi)/dy = -pi A sin(pi x) cos(pi y)
    v = +d(psi)/dx = +pi A cos(pi x) sin(pi y)
    """
    x = np.linspace(0.0, 2.0, nx)
    y = np.linspace(0.0, 1.0, ny)
    mx, my = np.meshgrid(x, y)
    u = -np.pi * amplitude * np.sin(np.pi * mx) * np.cos(np.pi * my)
    v = np.pi * amplitude * np.cos(np.pi * mx) * np.sin(np.pi * my)
    return VelocityField(x, y, u, v, name="double gyre")
