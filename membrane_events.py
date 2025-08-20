#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Membrane events on a hexagonal lattice have:
- Amphiphile swaps (liquid behavior)
- Peptide-triad insertion (3 touching hexes)
- Membrane growth: when displaced amphiphiles >= 4n+2 -> grow to (n+2)x(n+2)
- Save figures only (no plt.show)

Coordinate system: axial (q, r). No float drift.

Author: you
"""

from __future__ import annotations
import math
import os
import random
import argparse
from dataclasses import dataclass, field
from typing import Dict, Set, Tuple, List, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")  # ensure we only save figures
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

# --------------------------- Hex math (axial) -----------------------------

AXIAL_DIRS = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]  # cube-q,r neighbors

def add(a: Tuple[int, int], b: Tuple[int, int]) -> Tuple[int, int]:
    return (a[0] + b[0], a[1] + b[1])

def neighbors(ax: Tuple[int, int]) -> List[Tuple[int, int]]:
    q, r = ax
    return [(q + dq, r + dr) for (dq, dr) in AXIAL_DIRS]

def triangle_cells(center: Tuple[int, int], dir_index: int) -> List[Tuple[int, int]]:
    """
    A compact 3-hex triangle: center + neighbor dir d + neighbor dir (d+1).
    All three touch pairwise.
    """
    d1 = AXIAL_DIRS[dir_index % 6]
    d2 = AXIAL_DIRS[(dir_index + 1) % 6]
    a = center
    b = add(center, d1)
    c = add(center, d2)
    return [a, b, c]

# --------------------------- Data structures ------------------------------

@dataclass
class Peptide:
    pid: int
    cells: List[Tuple[int, int]]           # 3 axial coords
    anchor_len: int                        # A#
    rod_len: int                           # T#
    orientation: str                       # "inside" | "outside"
    color_index: int                       # color bucket for plotting

@dataclass
class MembraneState:
    n: int                                          # logical half-width (rhombus n x n)
    amph: Dict[Tuple[int, int], int] = field(default_factory=dict)   # axial -> carbon count
    peptides: Dict[int, Peptide] = field(default_factory=dict)
    occupied: Set[Tuple[int, int]] = field(default_factory=set)      # cells used by peptides
    displaced_count: int = 0                               # amphiphiles displaced by peptides
    next_pid: int = 0

# --------------------------- Grid builders --------------------------------

def rhombus_coords(n: int) -> List[Tuple[int, int]]:
    """
    Return axial coords for an n x n rhombus (q in [0..n-1], r in [0..n-1]) shifted so it looks centered.
    We'll offset to keep old content when growing.
    """
    return [(q, r) for q in range(n) for r in range(n)]

def init_membrane(n: int, rng: random.Random, carbon_min: int = 8, carbon_max: int = 20) -> MembraneState:
    st = MembraneState(n=n)
    for ax in rhombus_coords(n):
        st.amph[ax] = rng.randint(carbon_min, carbon_max)
    return st

# --------------------------- Events ---------------------------------------

def swap_amphiphiles(st: MembraneState, rng: random.Random, attempts: int = 500) -> None:
    """
    Randomly swaps amphiphile carbon counts between neighboring membrane cells
    to simulate liquid mixing. Skips cells occupied by peptides.
    """
    amph = st.amph
    occ = st.occupied
    keys = list(amph.keys())
    if not keys:
        return
    for _ in range(attempts):
        a = rng.choice(keys)
        if a in occ:  # peptide occupies (amph removed)
            continue
        neighs = [nb for nb in neighbors(a) if nb in amph and nb not in occ]
        if not neighs:
            continue
        b = rng.choice(neighs)
        amph[a], amph[b] = amph[b], amph[a]

def find_triangular_slot(st: MembraneState, rng: random.Random, max_tries: int = 2000) -> Optional[List[Tuple[int, int]]]:
    """
    Find three touching cells (triangle) that are all present in amph and not occupied.
    """
    amph = st.amph
    occ = st.occupied
    candidates = [ax for ax in amph.keys() if ax not in occ]
    if not candidates:
        return None
    for _ in range(max_tries):
        center = rng.choice(candidates)
        d = rng.randrange(6)
        tri = triangle_cells(center, d)
        if all((c in amph) and (c not in occ) for c in tri):
            return tri
    return None

def insert_peptide(st: MembraneState, rng: random.Random,
                   anchor_range=(8, 20), rod_range=(1, 20),
                   inside_prob: float = 0.5) -> Optional[Peptide]:
    """
    Remove 3 amphiphiles in a triangular slot and place a peptide triad there.
    """
    tri = find_triangular_slot(st, rng)
    if tri is None:
        return None

    # Displace amphiphiles
    displaced_here = 0
    for c in tri:
        if c in st.amph:
            st.amph.pop(c, None)
            displaced_here += 1
        st.occupied.add(c)

    st.displaced_count += displaced_here

    pid = st.next_pid
    st.next_pid += 1

    anchor_len = rng.randint(*anchor_range)
    rod_len = rng.randint(*rod_range)
    orientation = "inside" if rng.random() < inside_prob else "outside"
    pep = Peptide(
        pid=pid,
        cells=tri,
        anchor_len=anchor_len,
        rod_len=rod_len,
        orientation=orientation,
        color_index=pid % 20
    )
    st.peptides[pid] = pep
    return pep

def grow_membrane_if_needed(st: MembraneState, rng: random.Random,
                            carbon_min: int = 8, carbon_max: int = 20) -> bool:
    """
    If displaced amphiphiles >= 4n+2, expand rhombus from n x n to (n+2) x (n+2).
    We keep existing coords where they are and 'add a ring' around the rhombus.
    Returns True if growth occurred.
    """
    threshold = 4 * st.n + 2
    if st.displaced_count < threshold:
        return False

    old_n = st.n
    new_n = st.n + 2
    st.n = new_n

    # Build new coordinate set
    new_coords = set(rhombus_coords(new_n))
    old_coords = set(rhombus_coords(old_n))

    # Everything already in amph stays as-is; fill only the newly-added ring cells
    added = new_coords - old_coords

    for ax in added:
        if ax in st.occupied:  # safety (unlikely)
            continue
        st.amph[ax] = rng.randint(carbon_min, carbon_max)

    # Reduce displaced count by one threshold “unit”
    st.displaced_count -= threshold
    return True

# --------------------------- Plotting --------------------------------------

def axial_to_xy(ax: Tuple[int, int], radius: float) -> Tuple[float, float]:
    """
    Pointy-top axial to pixel. orientation for RegularPolygon will be pi/6 (30°) so edges are horizontal-ish.
    """
    q, r = ax
    x = radius * math.sqrt(3) * (q + r / 2.0)
    y = radius * 1.5 * r
    return x, y

def plot_membrane(st: MembraneState,
                  outpath: str,
                  figsize=(7.0, 11.0),
                  hex_radius: float = 9.0,
                  font_scale: float = 0.8) -> None:
    """
    Save a figure of the current membrane + peptides.
    """
    plt.close("all")
    fig = plt.figure(figsize=figsize, dpi=150)
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.set_axis_off()

    # Color map for peptides
    cmap = plt.colormaps["tab20"]

    # Plot amphiphiles: light gray, label carbon #
    for axc, ccount in st.amph.items():
        x, y = axial_to_xy(axc, hex_radius)
        patch = RegularPolygon((x, y), numVertices=6, radius=hex_radius,
                               orientation=math.pi/6, facecolor="#f0f0f0",
                               edgecolor="black", linewidth=0.8)
        ax.add_patch(patch)
        ax.text(x, y, f"{ccount}", ha="center", va="center",
                fontsize=8 * font_scale, color="#333333")

    # Plot peptides
    for pep in st.peptides.values():
        color = cmap(pep.color_index / 20.0)
        for i, cell in enumerate(pep.cells):
            x, y = axial_to_xy(cell, hex_radius)
            patch = RegularPolygon((x, y), numVertices=6, radius=hex_radius,
                                   orientation=math.pi/6, facecolor=color,
                                   edgecolor="black", linewidth=1.2)
            ax.add_patch(patch)
            # Label once (on the "center" hex = first in list)
            if i == 0:
                ax.text(x, y, f"P{pep.pid}\n{pep.anchor_len}/{pep.rod_len}\n{pep.orientation}",
                        ha="center", va="center", fontsize=7 * font_scale, color="white")

    # Tidy bounds
    if st.amph:
        xs, ys = zip(*(axial_to_xy(axc, hex_radius) for axc in st.amph.keys()))
        pad = 4 * hex_radius
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)

    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)

# --------------------------- Public API (for Mathematica) ------------------

def api_swap(state: MembraneState, attempts: int, seed: Optional[int] = None) -> MembraneState:
    rng = random.Random(seed)
    swap_amphiphiles(state, rng, attempts=attempts)
    return state

def api_insert(state: MembraneState,
               anchor_min=8, anchor_max=20,
               rod_min=1, rod_max=20,
               inside_prob=0.5,
               seed: Optional[int] = None) -> Optional[int]:
    rng = random.Random(seed)
    pep = insert_peptide(state, rng, anchor_range=(anchor_min, anchor_max),
                         rod_range=(rod_min, rod_max),
                         inside_prob=inside_prob)
    return None if pep is None else pep.pid

def api_grow(state: MembraneState, seed: Optional[int] = None) -> bool:
    rng = random.Random(seed)
    return grow_membrane_if_needed(state, rng)

# --------------------------- CLI runner ------------------------------------

def run_sim(n: int,
            steps: int,
            swaps_per_step: int,
            insert_every: int,
            outdir: str,
            seed: Optional[int]) -> None:
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)

    state = init_membrane(n, rng)
    # Optionally plot initial frame
    plot_membrane(state, os.path.join(outdir, f"frame_000.png"))

    for t in range(1, steps + 1):
        # liquid mixing
        swap_amphiphiles(state, rng, attempts=swaps_per_step)

        # peptide insertion cadence
        if insert_every > 0 and (t % insert_every == 0):
            insert_peptide(state, rng)

        # growth check
        grew = grow_membrane_if_needed(state, rng)

        # save frame
        out = os.path.join(outdir, f"frame_{t:03d}.png")
        plot_membrane(state, out)

    # Also save a snapshot of scalar data for inspection (simple text)
    with open(os.path.join(outdir, "state_summary.txt"), "w") as f:
        f.write(f"n={state.n}\n")
        f.write(f"peptides={len(state.peptides)}\n")
        f.write(f"displaced_count={state.displaced_count}\n")
        f.write(f"amph_count={len(state.amph)}\n")

def parse_args():
    p = argparse.ArgumentParser(description="Hex membrane events simulator")
    p.add_argument("--n", type=int, default=20, help="initial rhombus size n (n x n)")
    p.add_argument("--steps", type=int, default=50, help="number of time steps")
    p.add_argument("--swaps", type=int, default=800, help="swap attempts per step")
    p.add_argument("--insert-every", type=int, default=10, help="insert peptide every k steps (0 = never)")
    p.add_argument("--outdir", type=str, default="frames", help="output directory for PNGs")
    p.add_argument("--seed", type=int, default=None, help="random seed")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    run_sim(n=args.n,
            steps=args.steps,
            swaps_per_step=args.swaps,
            insert_every=args.insert_every,
            outdir=args.outdir,
            seed=args.seed)

