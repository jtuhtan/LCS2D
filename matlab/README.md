# Legacy MATLAB implementation (deprecated)

> **⚠️ Deprecated.** These scripts are the original 2013 MATLAB
> implementation of LCS2D. They are kept for reference and reproducibility,
> but they are no longer maintained and are **3–5 orders of magnitude
> slower** than the Python engine at the repository root. New work should
> use the Python package — see the [main README](../README.md) and the
> [migration guide](../docs/migration_from_matlab.md).

## Script → replacement map

| Legacy script | Purpose | Python replacement |
|---|---|---|
| `gyrefield.m`, `gyrefield_20130208.m` | double-gyre test field | `lcs2d gyre` |
| `seed.m` | seed grid + `TriScatteredInterp` interpolants | automatic: `lcs2d.load_field()` resamples scattered data; `--res` / `--margin` control seeding |
| `advect.m` | Euler advection of primary + 4 secondary seeds | `lcs2d advect` (RK4 by default; `--integrator euler` reproduces the legacy scheme) |
| `distance.m` | primary–secondary separation distances | internal to the FSLE computation |
| `ftle.m` | FTLE from the flow-map Jacobian | `lcs2d ftle` |
| `fsle.m`, `fsleMAX.m`, `fslecalc.m` | FSLE first-crossing search | `lcs2d fsle` |
| `flowmap.m`, `flowmap2d.m` | flow-map displacement plots | trajectories from `lcs2d advect` |
| `ftleplot.m`, `particleview.m`, `gertplot.m` | plots / videos | PNG output of every command (`plots.py`) |
| `groupdist.m`, `ipdm.m` | nearest-neighbor helpers for SMS export | not ported (`scipy.spatial.cKDTree` covers the same need) |

`ipdm.m` is John D'Errico's inter-point distance matrix utility from the
MATLAB File Exchange, bundled here as a dependency of the legacy scripts.

The `.mat` files are sample data used by the original scripts. Note that
`uF.mat` / `vF.mat` store `TriScatteredInterp` objects, which only MATLAB
can read; the Python engine instead reads raw `x, y, u, v` arrays from
`.mat`, `.csv`, or `.npz` files.

## Known numerical issues (fixed in the Python engine)

1. **`ftle.m` central-difference scaling.** The expression
   `(xnRt - xnLt)./2*delta` parses in MATLAB as `((xnRt - xnLt)/2) * delta`
   — it *multiplies* by the perturbation distance instead of dividing by
   `2*delta`. Ridge *positions* were unaffected (the script rescaled the
   field to [0, 1]), but absolute FTLE values were off by a factor of
   `(2*delta)^2` inside the logarithm.
2. **FTLE normalization.** `ftle.m` used `log(lambda_max)/T`; the standard
   definition (Shadden et al. 2005) is `log(lambda_max)/(2*|T|)`. The Python
   engine uses the standard form, so its values are directly comparable to
   the literature.
3. **Forward-Euler integration.** The legacy advection was first-order.
   The Python engine defaults to 4th-order Runge–Kutta.
