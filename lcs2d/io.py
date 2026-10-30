"""Loading velocity fields from common file formats.

Supported inputs (dispatch is by file extension):

``.npz``
    NumPy archive with arrays ``x (nx,)``, ``y (ny,)``, ``u (ny, nx)``,
    ``v (ny, nx)`` — the package's native format, written by
    :meth:`lcs2d.fields.VelocityField.save` and ``lcs2d gyre``.

``.csv`` / ``.txt`` / ``.dat`` / ``.xyz``
    Plain-text table whose first four columns are ``x  y  u  v`` (comma
    or whitespace separated; lines starting with ``#`` and a single
    header row are skipped; further columns are ignored).  If the
    points form a complete regular grid they are used directly;
    otherwise they are treated as scattered samples and linearly
    resampled onto a regular grid (``grid_res`` sets the spacing).
    ``NaN`` velocities mark no-flow cells, exactly as in a ``.npz``.

``.mat``
    MATLAB file containing the variables named in ``mat_vars``
    (default ``x, y, u, v``).  Grid-shaped (``meshgrid`` or ``ndgrid``
    style 2D ``x``/``y``), vector-axis (``x (nx)``, ``y (ny)``,
    ``u (ny, nx)``) and column-vector (scattered) layouts are all
    handled.  A name may select one column of a matrix with
    ``name:k`` (1-based, as in MATLAB), e.g.
    ``("seedXY:1", "seedXY:2", "seedUV:1", "seedUV:2")``.  Old
    ``TriScatteredInterp`` objects cannot be read — save the raw
    coordinate and velocity arrays instead.
"""

from __future__ import annotations

import os
import warnings

import numpy as np

from .fields import VelocityField

__all__ = ["load_field"]

_TEXT_EXTS = (".csv", ".txt", ".dat", ".xyz")


def _table_to_field(x, y, u, v, grid_res, name) -> VelocityField:
    """Turn x,y,u,v sample columns into a field, gridding if needed."""
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    u = np.asarray(u, dtype=np.float64).ravel()
    v = np.asarray(v, dtype=np.float64).ravel()
    if not (x.size == y.size == u.size == v.size):
        raise ValueError(f"x, y, u, v must all have the same number of "
                         f"samples; got {x.size}, {y.size}, {u.size}, "
                         f"{v.size}")

    ux = np.unique(x)
    uy = np.unique(y)
    if ux.size * uy.size == x.size:
        # complete regular grid: sort rows into (ny, nx) arrays.  The
        # count alone is not proof (a duplicated point can stand in for
        # a missing one), so check that every node appears exactly once.
        order = np.lexsort((x, y))
        gx = x[order].reshape(uy.size, ux.size)
        gy = y[order].reshape(uy.size, ux.size)
        if (gx == ux[None, :]).all() and (gy == uy[:, None]).all():
            gu = u[order].reshape(uy.size, ux.size)
            gv = v[order].reshape(uy.size, ux.size)
            try:
                return VelocityField(ux, uy, gu, gv, name=name)
            except ValueError:
                pass  # e.g. non-uniform spacing -> fall through to resampling
    return VelocityField.from_scattered(x, y, u, v, grid_res=grid_res,
                                        name=name)


def _load_text(path: str, grid_res, name) -> VelocityField:
    table = None
    try:
        table = np.genfromtxt(path, delimiter=",", comments="#", ndmin=2)
        if table.shape[1] < 4 or np.isnan(table[:, :4]).all():
            table = None  # not comma separated
    except ValueError:
        table = None
    if table is None:
        table = np.genfromtxt(path, comments="#", ndmin=2)
    if table.ndim != 2 or table.shape[1] < 4:
        raise ValueError(
            f"{path}: expected at least 4 columns (x, y, u, v); "
            f"got shape {table.shape}")
    table = table[:, :4]
    if table.shape[0] and np.isnan(table[0]).all():
        table = table[1:]  # header row failed to parse as numbers
    # rows without coordinates are unusable; NaN velocities are kept and
    # become no-flow cells, as documented
    bad = np.isnan(table[:, :2]).any(axis=1)
    if bad.any():
        warnings.warn(f"{path}: dropped {int(bad.sum())} rows with "
                      f"non-numeric x or y", stacklevel=3)
        table = table[~bad]
    return _table_to_field(table[:, 0], table[:, 1], table[:, 2],
                           table[:, 3], grid_res, name)


def _mat_var(data, path: str, spec: str) -> np.ndarray:
    """Fetch ``name`` or one column ``name:k`` (1-based) from a .mat dict."""
    key, _, col = spec.partition(":")
    arr = np.asarray(data[key], dtype=np.float64)
    if not col:
        return arr
    try:
        k = int(col)
    except ValueError:
        raise ValueError(f"{path}: bad column in {spec!r}; use name:k with "
                         f"a 1-based column number") from None
    if arr.ndim != 2 or not 1 <= k <= arr.shape[1]:
        raise ValueError(f"{path}: {spec!r} asks for column {k}, but "
                         f"{key!r} has shape {arr.shape}")
    return arr[:, k - 1]


def _load_mat(path: str, mat_vars, grid_res, name) -> VelocityField:
    from scipy.io import loadmat

    data = loadmat(path)
    missing = [spec for spec in mat_vars if spec.partition(":")[0] not in data]
    if missing:
        available = {k: np.shape(data[k]) for k in data
                     if not k.startswith("__")}
        raise KeyError(
            f"{path}: variable(s) {missing} not found. Available (with "
            f"shapes): {available}. Use --mat-vars to name the x,y,u,v "
            f"variables (name:k picks column k of a matrix).")
    x, y, u, v = (_mat_var(data, path, spec) for spec in mat_vars)

    if u.ndim == 2 and min(u.shape) > 1 and v.shape == u.shape:
        if x.shape == u.shape and y.shape == u.shape:
            if np.ptp(x[0, :]) == 0 and np.ptp(y[:, 0]) == 0:
                # ndgrid-style: x varies down rows -> transpose to meshgrid
                x, y, u, v = x.T, y.T, u.T, v.T
            # meshgrid-style: rows share y, columns share x
            return VelocityField(x[0, :], y[:, 0], u, v, name=name)
        if x.size * y.size == u.size and min(x.shape) == 1 \
                and min(y.shape) == 1:
            # separate axis vectors
            x, y = x.ravel(), y.ravel()
            if u.shape == (x.size, y.size) and x.size != y.size:
                u, v = u.T, v.T  # stored as (nx, ny)
            return VelocityField(x, y, u, v, name=name)
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
    if ext in _TEXT_EXTS:
        return _load_text(path, grid_res, name)
    if ext == ".mat":
        return _load_mat(path, mat_vars, grid_res, name)
    raise ValueError(f"unrecognized field file type {ext!r} "
                     "(expected .npz, .csv, .txt, .dat, .xyz or .mat)")
