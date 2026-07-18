"""Loading velocity fields from common file formats.

Supported inputs (dispatch is by file extension):

``.npz``
    NumPy archive with arrays ``x (nx,)``, ``y (ny,)``, ``u (ny, nx)``,
    ``v (ny, nx)`` — the package's native format, written by
    :meth:`lcs2d.fields.VelocityField.save` and ``lcs2d gyre``.

``.csv`` / ``.txt`` / ``.dat``
    Plain-text table with four columns ``x  y  u  v`` (comma or
    whitespace separated; lines starting with ``#`` and a single header
    row are skipped).  If the points form a complete regular grid they
    are used directly; otherwise they are treated as scattered samples
    and linearly resampled onto a regular grid (``grid_res`` sets the
    spacing).

``.mat``
    MATLAB file containing the variables named in ``mat_vars``
    (default ``x, y, u, v``).  Grid-shaped (``meshgrid``-style 2D
    ``x``/``y``) and column-vector (scattered) layouts are both
    handled.  Old ``TriScatteredInterp`` objects cannot be read — save
    the raw coordinate and velocity arrays instead.
"""

from __future__ import annotations

import os

import numpy as np

from .fields import VelocityField

__all__ = ["load_field"]


def _table_to_field(x, y, u, v, grid_res, name) -> VelocityField:
    """Turn x,y,u,v sample columns into a field, gridding if needed."""
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    v = np.asarray(v, dtype=np.float64).ravel()

    ux = np.unique(x)
    uy = np.unique(y)
    if ux.size * uy.size == x.size:
        # complete regular grid: sort rows into (ny, nx) arrays
        order = np.lexsort((x, y))
        gu = u[order].reshape(uy.size, ux.size)
        gv = v[order].reshape(uy.size, ux.size)
        try:
            return VelocityField(ux, uy, gu, gv, name=name)
        except ValueError:
            pass  # e.g. non-uniform spacing -> fall through to resampling
    return VelocityField.from_scattered(x, y, u, v, grid_res=grid_res,
                                        name=name)


def _load_text(path: str, grid_res, name) -> VelocityField:
    try:
        table = np.genfromtxt(path, delimiter=",", comments="#",
                              skip_header=0)
        if np.isnan(table).all():
            raise ValueError
    except ValueError:
        table = None
    if table is None or table.ndim != 2 or table.shape[1] < 4:
        table = np.genfromtxt(path, comments="#")
    if table.ndim != 2 or table.shape[1] < 4:
        raise ValueError(
            f"{path}: expected at least 4 columns (x, y, u, v); "
            f"got shape {table.shape}")
    if np.isnan(table[0]).any():
        table = table[1:]  # header row failed to parse as numbers
    if np.isnan(table).any():
        n_bad = int(np.isnan(table).any(axis=1).sum())
        table = table[~np.isnan(table).any(axis=1)]
        print(f"note: dropped {n_bad} rows with non-numeric entries")
    return _table_to_field(table[:, 0], table[:, 1], table[:, 2],
                           table[:, 3], grid_res, name)


def _load_mat(path: str, mat_vars, grid_res, name) -> VelocityField:
    from scipy.io import loadmat

    data = loadmat(path)
    missing = [key for key in mat_vars if key not in data]
    if missing:
        available = [k for k in data if not k.startswith("__")]
        raise KeyError(
            f"{path}: variable(s) {missing} not found. Available: "
            f"{available}. Use --mat-vars to name the x,y,u,v variables.")
    x, y, u, v = (np.asarray(data[key], dtype=np.float64)
                  for key in mat_vars)
    if x.ndim == 2 and x.shape == u.shape and x.shape[0] > 1 and x.shape[1] > 1:
        # meshgrid-style: rows share y, columns share x
        return VelocityField(x[0, :], y[:, 0], u, v, name=name)
    return _table_to_field(x, y, u, v, grid_res, name)


def load_field(path: str,
               grid_res: float | None = None,
               mat_vars=("x", "y", "u", "v")) -> VelocityField:
    """Load a velocity field from ``.npz``, text, or ``.mat`` (see
    module docstring for the format details)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"velocity-field file not found: {path}")
    name = os.path.basename(path)
    ext = os.path.splitext(path)[1].lower()
    if ext == ".npz":
        return VelocityField.load(path)
    if ext in (".csv", ".txt", ".dat", ".xyz"):
        return _load_text(path, grid_res, name)
    if ext == ".mat":
        return _load_mat(path, mat_vars, grid_res, name)
    raise ValueError(f"unrecognized field file type {ext!r} "
                     "(expected .npz, .csv, .txt, .dat or .mat)")
