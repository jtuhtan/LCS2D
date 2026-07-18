#!/usr/bin/env bash
# Full LCS2D demo pipeline on the double-gyre test field.
# Run from anywhere:  bash examples/run_double_gyre.sh [output_dir]
set -euo pipefail

OUT="${1:-lcs2d_demo}"
mkdir -p "$OUT"
echo "== LCS2D double-gyre demo -> $OUT/"

# 1. generate the test field
lcs2d gyre -o "$OUT/gyre.npz" --plot

# 2. inspect it (prints grid, speeds, suggested dt)
lcs2d info "$OUT/gyre.npz"

# 3. forward FTLE: repelling structures (ridge on the separatrix x = 1)
lcs2d ftle "$OUT/gyre.npz" --t 10 --res 0.005 -o "$OUT/ftle_forward"

# 4. backward FTLE: attracting structures
lcs2d ftle "$OUT/gyre.npz" --t 10 --backward --res 0.005 -o "$OUT/ftle_backward"

# 5. FSLE with separation ratio 2
lcs2d fsle "$OUT/gyre.npz" --t-max 40 --alpha 2 --res 0.01 -o "$OUT/fsle"

# 6. particle trajectories, every 10th step kept
lcs2d advect "$OUT/gyre.npz" --t 10 --res 0.05 --save-every 10 -o "$OUT/trajectories"

echo
echo "== done. Figures and data are in $OUT/"
