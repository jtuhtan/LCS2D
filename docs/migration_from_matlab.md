# Migrating from the MATLAB scripts

The Python engine implements the same method as the 2013 MATLAB code —
seed grid, four auxiliary "secondary" particles, advection, flow-map
Jacobian, FTLE/FSLE — so the migration is mostly *deleting* steps. This
page maps each legacy workflow to its replacement and lists every place
where the numbers deliberately differ.

## Workflow mapping

**Legacy FTLE pipeline** (run in order, sharing workspace variables):

```
gyrefield.m  →  seed.m  →  advect.m  →  ftle.m  →  ftleplot.m
```

**Now:**

```bash
lcs2d gyre -o gyre.npz
lcs2d ftle gyre.npz --t 10 -o ftle_out        # seeds, advects, differentiates, plots
```

**Legacy FSLE pipeline:**

```
seed.m  →  advect.m  →  distance.m  →  fsle.m / fsleMAX.m
```

**Now:**

```bash
lcs2d fsle gyre.npz --t-max 40 --alpha 2 -o fsle_out
```

## Parameter dictionary

| MATLAB variable (script) | Python equivalent | Notes |
|---|---|---|
| `resSeed` (`seed.m`) | `--res` | seed-grid spacing |
| `offset = 2*resSeed` (`seed.m`) | `--margin` | boundary inset; now optional, default 0 |
| `interpMethod` (`seed.m`) | — | scattered data is linearly resampled to a grid once; bilinear sampling thereafter |
| `TS`, `dt` (`advect.m`) | `--t`, `--dt` | specify total time directly; `dt` defaults to a CFL-safe value |
| `delta = 0.001*resSeed` (`advect.m`) | `--delta` | FTLE default is `res/100` |
| `d0`, `df = 2*d0` or `3*d0` (`fsle.m`) | `--delta`, `--alpha` | `df = alpha * delta` |
| `magSeedThresh` filtering (`seed.m`) | — | not needed; zero-velocity seeds simply produce zero FTLE/FSLE |
| CFL check (`advect.m`) | automatic | the default `dt` satisfies it by construction |

## Numerical differences (all intentional)

1. **Integrator.** Legacy: forward Euler. Now: RK4 by default.
   `--integrator euler` restores the legacy scheme if you need a
   like-for-like comparison.
2. **FTLE Jacobian scaling.** `ftle.m` computed
   `(xnRt - xnLt)./2*delta`, which MATLAB parses as
   `((xnRt-xnLt)/2)*delta` — multiplying by δ instead of dividing by 2δ.
   The Python engine divides by 2δ as the central difference requires.
   Ridge locations are unaffected; absolute values change.
3. **FTLE normalization.** `ftle.m` used `ln(λmax)/T`; the standard
   definition is `ln(λmax)/(2|T|)`. Python uses the standard form, so
   legacy values are ~2× larger (before the δ issue above).
4. **FSLE τ selection.** `fsleMAX.m` used the time of *maximum* observed
   separation; `fsle.m` used the *first threshold crossing*. The Python
   engine implements the first-crossing definition (the standard one),
   with the addition that the search stops per seed as soon as any of the
   four auxiliary particles crosses.
5. **Out-of-domain behavior.** Legacy `TriScatteredInterp` returned NaN
   outside the convex hull and NaNs were later zeroed in `ftle.m`; the
   Python sampler returns zero velocity outside the grid (particles
   freeze). Same practical effect, defined up front.

Because the legacy field was rescaled to [0, 1] before plotting, none of
these change the *structures you see* — they make the numbers match the
published definitions.

## Reading your old `.mat` data

Raw arrays (like `seedXY.mat` / `seedUV.mat`) load directly:

```bash
lcs2d info old_data.mat --mat-vars xyz,xyz,velocity,velocity   # inspect names first
```

If your variables are N×2/N×3 matrices rather than separate vectors,
export them once from MATLAB or Octave as columns:

```matlab
T = [xyz(:,1), xyz(:,2), velocity(:,1), velocity(:,2)];
writematrix(T, 'field.csv');    % or: dlmwrite('field.csv', T)
```

Saved `TriScatteredInterp` objects (`uF.mat`, `vF.mat`) cannot be read
outside MATLAB — re-export the underlying scatter points instead.

## Feature parity notes

* `flowmap.m` / `flowmap2d.m` displacement maps: compute from
  `lcs2d advect` trajectories (`hypot(x[-1]-x[0], y[-1]-y[0])`).
* `ftleplot.m` videos: save trajectory snapshots with `--save-every` and
  render with matplotlib's `FuncAnimation`, or plot per-time-step FTLE
  fields for a sequence of `--t` values.
* `groupdist.m` + `ipdm.m` (SMS export helpers): out of scope;
  `scipy.spatial.cKDTree` reproduces the nearest-neighbor mapping in a
  few lines.
