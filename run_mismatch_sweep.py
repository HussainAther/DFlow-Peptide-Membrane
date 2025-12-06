#!/usr/bin/env python3
"""
run_mismatch_sweep.py

Small helper script to run the reversible DFlow simulation for a range of
mismatch–mobility couplings (MISMATCH_MOBILITY_GAMMA).

Each value of gamma gets its own output directory under:

    runs/mismatch_sweep/gamma_<value>/

This is a thin wrapper around `dflow_reversible.py` using its CLI.
"""

import argparse
import subprocess
from pathlib import Path
import shlex


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sweep MISMATCH_MOBILITY_GAMMA for DFlow reversible runs."
    )

    # values of gamma to sweep
    p.add_argument(
        "--gammas",
        type=float,
        nargs="+",
        default=[0.0, 0.1, 0.25, 0.5],
        help="List of MISMATCH_MOBILITY_GAMMA values to sweep.",
    )

    # core simulation controls (you can tweak these defaults)
    p.add_argument("--N0", type=int, default=14, help="Half-width of rhombus.")
    p.add_argument(
        "--TOTAL_EVENTS",
        type=int,
        default=2000,
        help="Total number of stochastic events per run.",
    )
    p.add_argument(
        "--SAVE_STEPS",
        type=int,
        nargs="*",
        default=None,
        help="Steps at which to save frames (default: [0, TOTAL_EVENTS/2, TOTAL_EVENTS]).",
    )

    # output root
    p.add_argument(
        "--OUT_ROOT",
        type=str,
        default="runs/mismatch_sweep",
        help="Root directory for all sweep outputs.",
    )

    # utility
    p.add_argument(
        "--dry_run",
        action="store_true",
        help="Print commands instead of running them.",
    )

    return p.parse_args()


def main() -> None:
    args = parse_args()

    out_root = Path(args.OUT_ROOT)
    out_root.mkdir(parents=True, exist_ok=True)

    # default save steps if not provided
    if args.SAVE_STEPS is None:
        mid = args.TOTAL_EVENTS // 2
        save_steps = [0, mid, args.TOTAL_EVENTS]
    else:
        save_steps = args.SAVE_STEPS

    for gamma in args.gammas:
        # make a nice directory name like gamma_0p25
        gamma_tag = f"{gamma:.3f}".replace(".", "p")
        out_dir = out_root / f"gamma_{gamma_tag}"
        out_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            "python",
            "dflow_reversible.py",
            "--N0",
            str(args.N0),
            "--TOTAL_EVENTS",
            str(args.TOTAL_EVENTS),
            "--OUT",
            str(out_dir),
            "--SAVE_STEPS",
            *[str(s) for s in save_steps],
            "--MISMATCH_MOBILITY_GAMMA",
            str(gamma),
            "--COLOR_MODE",
            "mismatch",
            "--MISMATCH_CMAP",
            "viridis",
            "--MISMATCH_NORMALIZE",
        ]

        if args.dry_run:
            print(" ".join(shlex.quote(c) for c in cmd))
        else:
            print(f"\n=== Running gamma={gamma} → {out_dir} ===")
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()

