# Using your own data

The whole point of LCS2D is to run on *your* velocity fields — river
reaches, fishways, PIV planes, coastal models. This page covers the input
formats, how to choose each parameter, and what to do when the result looks
wrong.

## What you need

1. x, y coordinates of points in a plane, and
2. the velocity components u, v at those points,

for a **steady** (time-averaged or frozen) 2D flow. Units are up to you —
if positions are in metres and velocities in m/s, FTLE/FSLE come out in
1/s. (Time-resolved fields are not supported; computing LCS of an
unsteady field from one frozen snapshot is a common but real
approximation — interpret accordingly.)

## Input formats

**Text (`.csv`, `.txt`, `.dat`)** — four columns `x y u v`, comma- or
whitespace-separated; `#` comments and a single header row are fine:

```
x,y,u,v
0.00,0.00,0.012,-0.034
0.05,0.00,0.015,-0.031
...
```

If the rows form a complete regular grid it is used exactly as given.
Otherwise the points are treated as **scattered** and linearly resampled
onto a regular grid once, at spacing `--grid-res` (default: automatic,
about one grid node per sample). This replaces the old
`TriScatteredInterp` workflow — the expensive triangulation happens once,
not at every time step.

**MATLAB (`.mat`)** — provide the variable names if they differ from
`x,y,u,v`:

```bash
lcs2d ftle survey.mat --mat-vars xGrid,yGrid,uVel,vVel --t 30 -o out/ftle
```

Both `meshgrid`-style 2D arrays and column vectors (scattered) work.

**NumPy (`.npz`)** — the native format: arrays `x (nx,)`, `y (ny,)`,
`u (ny, nx)`, `v (ny, nx)`.

**Missing data:** `NaN` velocities are treated as no-flow cells (land,
masked PIV vectors, structures). Particles entering them — or leaving the
domain — simply freeze in place.

## Always start with `info`

```bash
lcs2d info myflow.csv --plot
```

This prints the grid, domain, speed range, and masked-cell count, suggests
a CFL-safe time step, and (with `--plot`) draws the field so you can spot
unit mistakes or bad vectors *before* spending compute on them.

## Choosing the parameters

| Parameter | Meaning | Guidance |
|---|---|---|
| `--t` (FTLE) | integration time | The single most important choice. Too short: no structures yet. Too long: everything saturates or leaves the domain. Start with the time it takes the mean flow to cross ~1/4 of your domain, then double/halve and compare — ridges that persist across `t` choices are the real ones. |
| `--backward` | time direction | Forward = repelling LCS (where particles will separate *from*), backward = attracting LCS (where particles collect). For transport questions you usually want both. |
| `--dt` | advection step | Default is CFL-safe (`0.5·min(dx,dy)/max_speed`) and rarely needs touching. Halve it once; if the field doesn't change, it was fine. |
| `--res` | seed spacing | Sets output resolution and cost (seeds × 5 particles). Default = velocity grid spacing. Finer `--res` sharpens ridges without changing where they are. |
| `--delta` (FTLE) | auxiliary offset | Default `res/100`. Keep it well below `res`; the result should be insensitive to it (that insensitivity is a good self-check). |
| `--delta` (FSLE) | initial separation d0 | This one is *physical*: it is the scale you probe. Default `res/2`. |
| `--alpha` (FSLE) | separation ratio | 2 is standard; larger values probe larger-scale stretching and need a longer `--t-max`. |
| `--t-max` (FSLE) | search horizon | Seeds that never separate within `t-max` get λ = 0. If the printed "never separated" count is large, increase it. |
| `--margin` | seed inset | For resampled scattered data the field near the convex-hull edge is unreliable — inset the seeds by a few grid spacings (the legacy `seed.m` used `2·res`). |

## A worked example

Say `river.csv` holds scattered ADCP-derived velocities over a 40 m reach
with speeds around 0.5 m/s:

```bash
lcs2d info river.csv --grid-res 0.25 --plot
lcs2d ftle river.csv --grid-res 0.25 --t 60  --res 0.25 --margin 1.0 -o out/river_fwd
lcs2d ftle river.csv --grid-res 0.25 --t 60  --backward --res 0.25 --margin 1.0 -o out/river_bwd
lcs2d fsle river.csv --grid-res 0.25 --t-max 240 --delta 0.25 --margin 1.0 -o out/river_fsle
```

Overlaying `river_fwd` and `river_bwd` shows the full transport skeleton:
repelling ridges upstream of obstacles, attracting ridges along
recirculation boundaries and wake shear layers.

## Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| "x must be uniformly spaced" | Your grid is stretched/irregular. Load it as scattered data (rename to `.csv`-style columns) and let `--grid-res` resample it. |
| Field is ~zero everywhere | `--t` too short for your velocities, or units mismatch (mm vs m). Check `lcs2d info`'s speed range. |
| Bright band hugging the domain edge | Particles freeze at the boundary, which mimics stretching. Use `--margin` to inset the seeds, and treat near-edge ridges with suspicion. |
| Ridges move when you change `--delta` (FTLE) | `delta` is too large relative to `--res`, or too close to machine precision. Keep `res/1000 ≤ delta ≤ res/10`. |
| FSLE full of zeros | `--t-max` too short (see the "never separated" count in the output). |
| Slow | `pip install -e .[fast]` (Numba) — typically ~50× on top of NumPy. The first Numba call compiles (a few seconds); later calls are cached. |
| Scattered load looks blocky | `--grid-res` too coarse; refine it (cost grows only at load time). |

## Getting results into other tools

* `<out>.csv` — `x,y,value` rows: Excel, R, pandas, QGIS (as delimited
  text layer), Tecplot.
* `<out>.npz` — `numpy.load()`; keys `x`, `y`, `values`, JSON `params`,
  plus `extra_tau` for FSLE.
* PNG figures are 150 dpi; regenerate at other settings from the `.npz`
  with `lcs2d.plots.plot_result` or your own matplotlib code.
