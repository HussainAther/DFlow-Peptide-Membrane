#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Raft diffusion diagnostic: R4 drift + mismatch-coupled R2

This script runs short reversible simulations under different values of
MISMATCH_MOBILITY_GAMMA and compares the raft–centroid step distributions
(i.e. a diffusion histogram).

Goal:
  - Make sure that adding mismatch-dependent mobility to R2 (and R4)
    does not introduce any strange coupling in the raft diffusion statistics.

Usage (from repo root):

  python -m analysis.raft_diffusion_interplay \
    --out runs/diffusion_diag \
    --steps 8000 \
    --sample-every 10 \
    --gammas 0.0 0.5

This will produce:
  - diffusion_steps_gamma_*.csv
  - diffusion_hist_gamma_compare.png
"""

import os
import math
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

# Make repo root importable when running from analysis/
HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]
import sys
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Import reversible simulator pieces
from dflow_reversible import (           # type: ignore
    Config,
    Membrane,
    Diurnal,
    DiurnalEnv,
    Scheduler,
    Event,
    AX_NEI,
    rate_swap, do_swap,
    rate_thicken, do_thicken,
    rate_thin, do_thin,
    rate_polymerize, do_polymerize,
    rate_depoly, do_depoly,
    rate_insert, do_insert,
    rate_desorb, do_desorb,
    rate_pept_step, do_pept_step,
    rate_flip, do_flip,
    rate_assoc, do_assoc,
    rate_dissoc, do_dissoc,
    rate_R4_drift, do_R4_drift,
)


# ---------- helper: centroids and step-lengths ----------

def _raft_centroids_from_mem(mem):
    """
    Return {raft_id: (q_mean, r_mean)} in axial coords.
    We don't care about exact composition, just the instantaneous centroid.
    """
    buckets = defaultdict(list)
    for pid, p in mem.peptides.items():
        rid = p["raft"]
        q, r = p["cent"]
        buckets[rid].append((q, r))

    centroids = {}
    for rid, pts in buckets.items():
        qs, rs = zip(*pts)
        centroids[rid] = (np.mean(qs), np.mean(rs))
    return centroids


def _axial_to_xy(q, r):
    """Same mapping as in dflow_reversible.axial_to_xy, but inlined to avoid import cycles."""
    sqrt3 = math.sqrt(3.0)
    x = 1.0 * (3.0 / 2.0) * q
    y = 1.0 * sqrt3 * (r + q / 2.0)
    return float(x), float(y)


def run_single_diffusion_run(
    gamma: float,
    cfg_template: Config,
    n_steps: int = 8000,
    sample_every: int = 10,
    enable_R4: bool = True,
    seed_offset: int = 0,
):
    """
    Run one short reversible simulation and collect raft step lengths
    at a fixed sampling interval.

    Returns:
        np.ndarray of step lengths Δr (in xy units) for all rafts / samples.
    """
    # Clone template cfg and set mismatch gamma
    cfg = Config(**vars(cfg_template))
    cfg.MISMATCH_MOBILITY_GAMMA = gamma
    cfg.SEED = cfg.SEED + seed_offset

    # init RNG + membrane + env
    np.random.seed(cfg.SEED)
    import random as pyrand
    pyrand.seed(cfg.SEED)

    mem = Membrane(n=cfg.N0, hex_radius=cfg.HEX_RADIUS, cfg=cfg)
    env = DiurnalEnv(
        diurnal=Diurnal(cfg.DAY_STEPS, cfg.NIGHT_STEPS,
                        cfg.POLY_GAIN_DAY, cfg.POLY_GAIN_NIGHT),
        beta=cfg.BETA,
    )

    # scheduler with same events as main reversible driver
    sched = Scheduler()
    sched.register(Event("A1_swap",          rate_swap,       do_swap))
    sched.register(Event("A3_plus_thicken",  rate_thicken,    do_thicken))
    sched.register(Event("A3_minus_thin",    rate_thin,       do_thin))
    sched.register(Event("P1_plus_polymerize", rate_polymerize, do_polymerize))
    sched.register(Event("P1_minus_depoly",    rate_depoly,     do_depoly))
    sched.register(Event("P2_plus_insert",     rate_insert,     do_insert))
    sched.register(Event("P2_minus_desorb",    rate_desorb,     do_desorb))
    sched.register(Event("P3_step",          rate_pept_step,  do_pept_step))
    sched.register(Event("P4_flip",          rate_flip,       do_flip))
    sched.register(Event("R1_plus_assoc",    rate_assoc,      do_assoc))
    sched.register(Event("R1_minus_dissoc",  rate_dissoc,     do_dissoc))

    if enable_R4:
        sched.register(Event("R4_drift",     rate_R4_drift,   do_R4_drift))

    # trajectories: keep centroid at sample times and convert to step lengths
    last_centroids = None
    step_lengths = []

    for step in range(1, n_steps + 1):
        _ = sched.step(mem, env)
        env.diurnal.tick()

        if (step % sample_every) == 0:
            centroids = _raft_centroids_from_mem(mem)

            if last_centroids is not None:
                # match rafts that exist at both times
                common_ids = set(centroids.keys()) & set(last_centroids.keys())
                for rid in common_ids:
                    q0, r0 = last_centroids[rid]
                    q1, r1 = centroids[rid]
                    x0, y0 = _axial_to_xy(q0, r0)
                    x1, y1 = _axial_to_xy(q1, r1)
                    dr = math.hypot(x1 - x0, y1 - y0)
                    step_lengths.append(dr)

            last_centroids = centroids

    return np.array(step_lengths, dtype=float)


# ---------- CLI + plotting ----------

def main():
    ap = argparse.ArgumentParser(
        description="Compare raft diffusion histograms for R4 drift + mismatch-coupled rates."
    )
    ap.add_argument(
        "--out",
        type=str,
        default="runs/diffusion_diag",
        help="Output directory for CSV + PNG.",
    )
    ap.add_argument(
        "--steps",
        type=int,
        default=8000,
        help="Number of KMC steps per run.",
    )
    ap.add_argument(
        "--sample-every",
        type=int,
        default=10,
        help="Sampling stride for centroid logging (in KMC steps).",
    )
    ap.add_argument(
        "--gammas",
        type=float,
        nargs="+",
        default=[0.0, 0.5],
        help="List of MISMATCH_MOBILITY_GAMMA values to compare.",
    )
    ap.add_argument(
        "--replicates",
        type=int,
        default=3,
        help="Number of independent runs per gamma.",
    )
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # base config template (shares all other settings with main driver)
    cfg_template = Config()

    all_hist_data = {}

    for gamma in args.gammas:
        all_steps = []
        for rep in range(args.replicates):
            sl = run_single_diffusion_run(
                gamma=gamma,
                cfg_template=cfg_template,
                n_steps=args.steps,
                sample_every=args.sample_every,
                enable_R4=True,
                seed_offset=rep * 1000,
            )
            all_steps.append(sl)

        steps_gamma = np.concatenate(all_steps) if all_steps else np.zeros(0)
        all_hist_data[gamma] = steps_gamma

        # save raw step lengths for this gamma
        np.savetxt(
            outdir / f"diffusion_steps_gamma_{gamma:.3f}.csv",
            steps_gamma,
            delimiter=",",
            header="dr_xy",
            comments="",
        )

    # --------- plot overlayed histograms ---------
    plt.figure(figsize=(6, 4))
    bins = 40

    for gamma, steps in all_hist_data.items():
        if steps.size == 0:
            continue
        plt.hist(
            steps,
            bins=bins,
            density=True,
            alpha=0.4,
            label=f"γ = {gamma}",
            histtype="stepfilled",
        )

    plt.xlabel("Raft centroid step length Δr (lattice units)")
    plt.ylabel("Probability density")
    plt.title("R4 drift + mismatch coupling: raft diffusion histogram")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(outdir / "diffusion_hist_gamma_compare.png", dpi=300)
    plt.close()

    print(f"[diffusion_interplay] Saved results to: {outdir}")


if __name__ == "__main__":
    main()
