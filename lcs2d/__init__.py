"""lcs2d — Lagrangian coherent structures from steady 2D velocity fields.

Fast, vectorized Python successor to the original LCS2D MATLAB scripts.

Typical library use::

    import lcs2d

    field = lcs2d.double_gyre()
    result = lcs2d.ftle_field(field, t=10.0)

or from the shell::

    lcs2d gyre -o gyre.npz
    lcs2d ftle gyre.npz --t 10 -o ftle_out
"""

from ._kernels import NUMBA_AVAILABLE
from .advect import Trajectories, advect
from .fields import VelocityField, coastal_bay, double_gyre
from .fsle import fsle_field
from .ftle import ftle_field, seed_grid
from .io import load_field
from .result import LCSResult

__version__ = "2.0.0"

__all__ = [
    "VelocityField", "double_gyre", "coastal_bay", "load_field",
    "advect", "Trajectories", "seed_grid",
    "ftle_field", "fsle_field", "LCSResult",
    "NUMBA_AVAILABLE", "__version__",
]
