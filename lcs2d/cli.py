"""Command-line interface: drive every LCS computation from the shell.

Run ``lcs2d --help`` (or ``python -m lcs2d --help``) for the command
list, and ``lcs2d <command> --help`` for the options of each command.

A typical first session::

    lcs2d gyre -o gyre.npz            # make the double-gyre test field
    lcs2d coast -o coast.npz          # ... or the coastal cape-and-bay field
    lcs2d ftle gyre.npz --t 10 -o out/ftle_fwd
    lcs2d fsle gyre.npz --t-max 40 -o out/fsle
    lcs2d selftest                    # verify the install end to end
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

from . import __version__
from ._kernels import NUMBA_AVAILABLE
from .advect import advect, resolve_backend
from .fields import VelocityField, coastal_bay, double_gyre
from .fsle import fsle_field
from .ftle import ftle_field
from .io import load_field
from .result import LCSResult


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def _say(args, *msg) -> None:
    if not args.quiet:
        print(*msg)


def _ensure_parent(path: str) -> None:
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)


def _load(args) -> VelocityField:
    mat_vars = tuple(s.strip() for s in args.mat_vars.split(","))
    if len(mat_vars) != 4:
        sys.exit("error: --mat-vars needs exactly 4 comma-separated names "
                 "(x,y,u,v)")
    try:
        field = load_field(args.field, grid_res=args.grid_res,
                           mat_vars=mat_vars)
    except Exception as exc:  # present file/format problems politely
        sys.exit(f"error loading {args.field}: {exc}")
    _say(args, field.summary())
    _say(args, "")
    return field


def _save_result(args, result: LCSResult) -> None:
    _ensure_parent(args.out)
    written = []
    result.save_npz(args.out + ".npz")
    written.append(args.out + ".npz")
    if not args.no_csv:
        result.save_csv(args.out + ".csv")
        written.append(args.out + ".csv")
    if not args.no_plot:
        from .plots import plot_result
        plot_result(result, args.out + ".png", cmap=args.cmap)
        written.append(args.out + ".png")
    _say(args, "")
    _say(args, result.summary())
    _say(args, "")
    _say(args, "written:")
    for path in written:
        _say(args, f"  {path}")


def _add_field_arg(p: argparse.ArgumentParser) -> None:
    p.add_argument("field",
                   help="velocity-field file (.npz, .csv/.txt/.dat with "
                        "x,y,u,v columns, or .mat)")
    p.add_argument("--grid-res", type=float, default=None,
                   help="grid spacing used when resampling scattered "
                        "input (default: automatic)")
    p.add_argument("--mat-vars", default="x,y,u,v",
                   help="for .mat input: comma-separated names of the "
                        "x,y,u,v variables (default: %(default)s)")


def _add_common_lcs_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--dt", type=float, default=None,
                   help="advection time step (default: CFL-based - one "
                        "grid cell per step for ftle, half a cell for "
                        "fsle, where dt sets the tau resolution)")
    p.add_argument("--res", type=float, default=None,
                   help="seed-grid spacing (default: the velocity grid "
                        "spacing)")
    p.add_argument("--margin", type=float, default=0.0,
                   help="inset seeds this far from the domain edge "
                        "(default: 0; useful for resampled scattered data)")
    p.add_argument("--integrator", choices=("rk4", "euler"), default="rk4",
                   help="time integrator (default: rk4; euler reproduces "
                        "the legacy MATLAB scheme)")
    p.add_argument("--backend", choices=("auto", "numpy", "numba"),
                   default="auto",
                   help="compute backend (default: auto = numba when "
                        "installed)")
    p.add_argument("-o", "--out", default="lcs_out",
                   help="output prefix; writes <out>.npz, <out>.csv and "
                        "<out>.png (default: %(default)s)")
    p.add_argument("--no-plot", action="store_true",
                   help="skip the PNG plot")
    p.add_argument("--no-csv", action="store_true",
                   help="skip the CSV table")
    p.add_argument("--cmap", default=None,
                   help="matplotlib colormap for the plot")


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------

def cmd_gyre(args) -> None:
    _write_field(args, double_gyre(nx=args.nx, ny=args.ny,
                                   amplitude=args.amplitude))


def cmd_coast(args) -> None:
    _write_field(args, coastal_bay(res=args.res, current=args.current))


def _write_field(args, field: VelocityField) -> None:
    if not args.out.lower().endswith(".npz"):
        args.out += ".npz"  # what np.savez would silently do anyway
    _ensure_parent(args.out)
    field.save(args.out)
    _say(args, field.summary())
    _say(args, "")
    _say(args, f"written:\n  {args.out}")
    if args.plot:
        from .plots import plot_velocity
        png = os.path.splitext(args.out)[0] + ".png"
        plot_velocity(field, png)
        _say(args, f"  {png}")


def cmd_info(args) -> None:
    field = _load(args)
    if field.max_speed() > 0.0:
        _say(args, f"suggested dt (CFL 0.5): {field.cfl_dt():.6g}")
    else:
        _say(args, "warning: the field has zero velocity everywhere, so "
                   "there is nothing to advect (check units and masking)")
    if args.plot:
        png = os.path.splitext(args.field)[0] + "_field.png"
        from .plots import plot_velocity
        plot_velocity(field, png)
        _say(args, f"written:\n  {png}")


def cmd_ftle(args) -> None:
    field = _load(args)
    t = -abs(args.t) if args.backward else args.t
    tic = time.perf_counter()
    result = ftle_field(field, t=t, dt=args.dt, res=args.res,
                        delta=args.delta, margin=args.margin,
                        integrator=args.integrator, backend=args.backend)
    toc = time.perf_counter() - tic
    n_seeds = result.values.size
    _say(args, f"FTLE: {n_seeds} seeds x 5 particles advected in "
               f"{toc:.2f} s ({resolve_backend(args.backend)} backend)")
    _save_result(args, result)


def cmd_fsle(args) -> None:
    field = _load(args)
    t_max = -abs(args.t_max) if args.backward else args.t_max
    tic = time.perf_counter()
    result = fsle_field(field, t_max=t_max, dt=args.dt, res=args.res,
                        delta=args.delta, alpha=args.alpha,
                        margin=args.margin, integrator=args.integrator,
                        backend=args.backend)
    toc = time.perf_counter() - tic
    _say(args, f"FSLE: {result.values.size} seeds x 5 particles in "
               f"{toc:.2f} s ({resolve_backend(args.backend)} backend)")
    never = result.params["n_never_separated"]
    if never:
        _say(args, f"note: {never} seeds never reached the separation "
                   f"threshold within t_max (lambda set to 0 there — "
                   f"increase --t-max to probe them)")
    _save_result(args, result)


def cmd_advect(args) -> None:
    field = _load(args)
    from .ftle import seed_grid
    sx, sy = seed_grid(field, res=args.res, margin=args.margin)
    mx, my = np.meshgrid(sx, sy)
    tic = time.perf_counter()
    traj = advect(field, mx.ravel(), my.ravel(), t=args.t, dt=args.dt,
                  integrator=args.integrator, backend=args.backend,
                  record=True, save_every=args.save_every)
    toc = time.perf_counter() - tic
    _say(args, f"advected {mx.size} particles for t={args.t} in {toc:.2f} s "
               f"({traj.t.size} snapshots kept)")
    _ensure_parent(args.out)
    traj.save(args.out + ".npz")
    written = [args.out + ".npz"]
    if not args.no_plot:
        from .plots import plot_trajectories
        plot_trajectories(traj, args.out + ".png", field=field)
        written.append(args.out + ".png")
    _say(args, "written:")
    for path in written:
        _say(args, f"  {path}")


def cmd_selftest(args) -> None:
    """End-to-end verification on analytic flows.  Exit code 0 = healthy."""
    failures = []

    def check(label, ok, detail=""):
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {label}" + (f"  ({detail})" if detail else ""))
        if not ok:
            failures.append(label)

    print(f"lcs2d {__version__} selftest")
    print(f"  numba available: {NUMBA_AVAILABLE}")
    print()

    # 1) saddle flow u=ax, v=-ay has uniform analytic FTLE = a
    a = 0.5
    x = np.linspace(-1.0, 1.0, 101)
    y = np.linspace(-1.0, 1.0, 101)
    mx, my = np.meshgrid(x, y)
    saddle = VelocityField(x, y, a * mx, -a * my, name="saddle")
    res = ftle_field(saddle, t=1.0, dt=0.01, res=0.1, margin=0.55,
                     backend="numpy")
    err = np.abs(res.values - a).max() / a
    check("FTLE matches analytic saddle-flow value", err < 1e-4,
          f"max relative error {err:.2e}")

    # 2) double gyre: forward-FTLE ridge sits on the separatrix x = 1
    # (the domain boundary is itself a separatrix streamline, so the
    # search is restricted to the interior)
    gyre = double_gyre(nx=201, ny=101)
    res = ftle_field(gyre, t=8.0, res=0.02, backend="numpy")
    inner = (res.x > 0.3) & (res.x < 1.7)
    ridge_x = res.x[inner][np.argmax(res.values[:, inner].mean(axis=0))]
    check("double-gyre FTLE ridge at x = 1", abs(ridge_x - 1.0) < 0.1,
          f"ridge found at x = {ridge_x:.3f}")
    check("FTLE field is finite", bool(np.isfinite(res.values).all()))

    # 3) numba and numpy agree
    if NUMBA_AVAILABLE:
        res_nb = ftle_field(gyre, t=8.0, res=0.02, backend="numba")
        diff = np.abs(res_nb.values - res.values).max()
        check("numba and numpy backends agree", diff < 1e-9,
              f"max abs difference {diff:.2e}")

    # 4) FSLE on the saddle: lambda ~ a at the probed scale
    res = fsle_field(saddle, t_max=20.0, dt=0.01, res=0.1, delta=0.01,
                     alpha=2.0, margin=0.55, backend="numpy")
    lam = res.values[np.isfinite(res.values) & (res.values > 0)]
    err = abs(np.median(lam) - a) / a
    check("FSLE matches analytic saddle-flow value", err < 0.05,
          f"median relative error {err:.2e}")

    print()
    if failures:
        sys.exit(f"selftest FAILED: {', '.join(failures)}")
    print("selftest passed - the installation is healthy.")


# ---------------------------------------------------------------------------
# parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lcs2d",
        description="Lagrangian coherent structures (FTLE / FSLE) from "
                    "steady 2D velocity fields.",
        epilog="Start with: lcs2d gyre -o gyre.npz && "
               "lcs2d ftle gyre.npz --t 10 -o ftle_out")
    parser.add_argument("--version", action="version",
                        version=f"lcs2d {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("gyre", help="generate the double-gyre test field")
    p.add_argument("--nx", type=int, default=401)
    p.add_argument("--ny", type=int, default=201)
    p.add_argument("--amplitude", type=float, default=0.1,
                   help="streamfunction amplitude A (default: %(default)s)")
    p.add_argument("-o", "--out", default="gyre.npz",
                   help="output .npz path (default: %(default)s)")
    p.add_argument("--plot", action="store_true",
                   help="also write a PNG of the field")
    p.set_defaults(func=cmd_gyre)

    p = sub.add_parser("coast",
                       help="generate the coastal cape-and-bay field "
                            "(km, km/h)")
    p.add_argument("--res", type=float, default=0.05,
                   help="grid spacing in km (default: %(default)s)")
    p.add_argument("--current", type=float, default=0.36,
                   help="alongshore current in km/h "
                        "(default: %(default)s, i.e. 0.1 m/s)")
    p.add_argument("-o", "--out", default="coast.npz",
                   help="output .npz path (default: %(default)s)")
    p.add_argument("--plot", action="store_true",
                   help="also write a PNG of the field")
    p.set_defaults(func=cmd_coast)

    p = sub.add_parser("info", help="describe a velocity-field file")
    _add_field_arg(p)
    p.add_argument("--plot", action="store_true",
                   help="write a PNG of the field next to the input file")
    p.set_defaults(func=cmd_info)

    p = sub.add_parser("ftle",
                       help="compute a finite-time Lyapunov exponent field")
    _add_field_arg(p)
    p.add_argument("--t", type=float, required=True,
                   help="integration time; use --backward (or a negative "
                        "value, e.g. --t=-10) for attracting LCS")
    p.add_argument("--backward", action="store_true",
                   help="integrate backward in time (attracting LCS)")
    p.add_argument("--delta", type=float, default=None,
                   help="auxiliary-particle offset (default: res/100)")
    _add_common_lcs_args(p)
    p.set_defaults(func=cmd_ftle)

    p = sub.add_parser("fsle",
                       help="compute a finite-size Lyapunov exponent field")
    _add_field_arg(p)
    p.add_argument("--t-max", type=float, required=True,
                   help="longest time to wait for the separation threshold")
    p.add_argument("--backward", action="store_true",
                   help="integrate backward in time")
    p.add_argument("--alpha", type=float, default=2.0,
                   help="separation ratio df/d0 (default: %(default)s)")
    p.add_argument("--delta", type=float, default=None,
                   help="initial separation d0 (default: res/2)")
    _add_common_lcs_args(p)
    p.set_defaults(func=cmd_fsle)

    p = sub.add_parser("advect",
                       help="advect a particle grid and save trajectories")
    _add_field_arg(p)
    p.add_argument("--t", type=float, required=True,
                   help="integration time (negative = backward)")
    p.add_argument("--dt", type=float, default=None,
                   help="advection time step (default: CFL-based)")
    p.add_argument("--res", type=float, default=None,
                   help="particle-grid spacing (default: velocity grid "
                        "spacing)")
    p.add_argument("--margin", type=float, default=0.0,
                   help="inset particles from the domain edge")
    p.add_argument("--integrator", choices=("rk4", "euler"), default="rk4")
    p.add_argument("--backend", choices=("auto", "numpy", "numba"),
                   default="auto")
    p.add_argument("--save-every", type=int, default=1,
                   help="record every Nth step (default: every step)")
    p.add_argument("-o", "--out", default="trajectories",
                   help="output prefix (default: %(default)s)")
    p.add_argument("--no-plot", action="store_true")
    p.set_defaults(func=cmd_advect)

    p = sub.add_parser("selftest",
                       help="verify the installation against analytic flows")
    p.set_defaults(func=cmd_selftest)

    for sp in sub.choices.values():
        sp.add_argument("--quiet", action="store_true",
                        help="suppress progress output")
    return parser


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)
    try:
        args.func(args)
    except (ValueError, RuntimeError) as exc:
        # bad parameters or data: a one-line message, not a traceback
        sys.exit(f"error: {exc}")


if __name__ == "__main__":
    main()
