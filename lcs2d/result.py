"""Container for computed LCS scalar fields (FTLE, FSLE)."""

from __future__ import annotations

import json
from dataclasses import dataclass, field as dataclass_field

import numpy as np

__all__ = ["LCSResult"]


@dataclass
class LCSResult:
    """A scalar LCS diagnostic on a regular seed grid.

    Attributes
    ----------
    kind : str
        ``"ftle"`` or ``"fsle"``.
    x, y : ndarray
        Seed-grid coordinates, shapes ``(nx,)`` and ``(ny,)``.
    values : ndarray, shape (ny, nx)
        The FTLE or FSLE field, units 1/time.
    params : dict
        The parameters the field was computed with (t, dt, delta, ...).
    extra : dict of ndarray
        Optional companion arrays (e.g. the FSLE separation time
        ``tau``), each shaped like ``values``.
    """

    kind: str
    x: np.ndarray
    y: np.ndarray
    values: np.ndarray
    params: dict = dataclass_field(default_factory=dict)
    extra: dict = dataclass_field(default_factory=dict)

    def save_npz(self, path: str) -> None:
        """Save to compressed ``.npz`` (arrays + JSON-encoded params)."""
        np.savez_compressed(path, kind=np.array(self.kind),
                            x=self.x, y=self.y, values=self.values,
                            params=np.array(json.dumps(self.params)),
                            **{f"extra_{k}": v for k, v in self.extra.items()})

    @classmethod
    def load_npz(cls, path: str) -> "LCSResult":
        with np.load(path, allow_pickle=False) as data:
            extra = {k[len("extra_"):]: data[k]
                     for k in data.files if k.startswith("extra_")}
            return cls(kind=str(data["kind"]),
                       x=data["x"], y=data["y"], values=data["values"],
                       params=json.loads(str(data["params"])),
                       extra=extra)

    def save_csv(self, path: str) -> None:
        """Save as plain-text CSV with columns ``x, y, value`` (one row
        per seed), easy to load into any other tool."""
        mx, my = np.meshgrid(self.x, self.y)
        table = np.column_stack([mx.ravel(), my.ravel(), self.values.ravel()])
        header = f"x,y,{self.kind}"
        np.savetxt(path, table, delimiter=",", header=header, comments="")

    def summary(self) -> str:
        finite = self.values[np.isfinite(self.values)]
        lines = [
            f"{self.kind.upper()} field: {self.values.shape[1]} x "
            f"{self.values.shape[0]} seeds",
            f"range   : min {finite.min():.6g}, max {finite.max():.6g}, "
            f"mean {finite.mean():.6g}  [1/time]",
        ]
        for key, val in self.params.items():
            lines.append(f"{key:8}: {val}")
        return "\n".join(lines)
