# membrane_events_case_loop_fixed.py
# -*- coding: utf-8 -*-
import os
import json
import math
import random
import numpy as np
import matplotlib as mpl

mpl.use("Agg")  # save-only
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

# -----------------------------
# Config
# -----------------------------
SEED = 42
random.seed(SEED)
np.random.seed(SEED)

OUT_DIR = "frames"
os.makedirs(OUT_DIR, exist_ok=True)

N0 = 20  # initial rhombus width/height (axial coords range)
HEX_RADIUS = 1.0  # RegularPolygon radius
BORDER = 1.0  # figure padding
LABEL_CARBONS = False  # set True to print carbon count on each amphiphile hex
SAVE_EVERY = 100  # save every k events
TOTAL_EVENTS = 5000

EVENT_WEIGHTS = {
    "amph_swap": 1.0,           # swap neighboring amphiphiles (liquid membrane)
    "pept_insert": 0.15,        # insert a 3-hex peptide (triad)
    "pept_diffuse": 0.8,        # try diffusive step for a random peptide
    "raft_adhere_merge": 0.6,   # if two rafts touch, merge ids
    "pept_flip": 0.2,           # toggle inside/outside flag
}

# -----------------------------
# Hex helpers (axial coords, flat-top)
# axial (q,r) -> pixel (x,y)
# -----------------------------
SQRT3 = math.sqrt(3.0)


def axial_to_xy(q, r, s=HEX_RADIUS):
    x = s * (3 / 2) * q
    y = s * SQRT3 * (r + q / 2)
    return x, y


# 6 axial neighbors (point-to-point touch)
AX_NEI = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]

# For a triad (three hexes meeting at a vertex), we use a canonical chevron:
# (0,0), (1,0), (0,1) rotated by 6 possible axial rotations via symmetry.
TRIAD_SHAPES = [
    [(0, 0), (1, 0), (0, 1)],
    [(0, 0), (0, 1), (-1, 1)],
    [(0, 0), (-1, 1), (-1, 0)],
    [(0, 0), (-1, 0), (0, -1)],
    [(0, 0), (0, -1), (1, -1)],
    [(0, 0), (1, -1), (1, 0)],
]


def rotate_triad(anchor, orient):
    aq, ar = anchor
    offs = TRIAD_SHAPES[orient % 6]
    return [(aq + dq, ar + dr) for (dq, dr) in offs]


# -----------------------------
# Colormap helper (robust for old/new Matplotlib)
# -----------------------------
def get_palette(n, name="tab20"):
    try:
        cmap = mpl.colormaps[name]  # Matplotlib >= 3.7
    except AttributeError:
        cmap = plt.get_cmap(name)  # fallback older

    # listed colormap? (tab10/tab20/etc)
    if hasattr(cmap, "colors"):
        base = list(cmap.colors)
        if n <= len(base):
            return base[:n]
        reps = (n + len(base) - 1) // len(base)
        return (base * reps)[:n]

    # continuous
    if n == 1:
        return [cmap(0.0)]
    return [cmap(i / (n - 1)) for i in range(n)]


# -----------------------------
# Membrane state
# -----------------------------
class Membrane:
    def __init__(self, n=N0):
        # Rhombus region in axial coords: -n..+n with |q|, |r|, |q+r| <= n (hex range)
        self.n = n
        self.amph = {}  # (q,r) -> dict(carbon:int, pep:bool)
        self.peptides = {}  # pid -> dict(cent:(q,r), orient:int, raft:int, inside:bool)
        self.pid_counter = 0
        self.raft_counter = 0
        self.displaced_amph = 0  # global counter to trigger growth
        self.event_log = []
        self._init_membrane()

    def _axial_in_rhombus(self, q, r):
        # Use hex range: |q| <= n, |r| <= n, |q+r| <= n
        n = self.n
        return (abs(q) <= n) and (abs(r) <= n) and (abs(q + r) <= n)

    def _init_membrane(self):
        # Fill with amphiphiles (random carbon length 10..20 for variety)
        for q in range(-self.n, self.n + 1):
            for r in range(-self.n, self.n + 1):
                if self._axial_in_rhombus(q, r):
                    self.amph[(q, r)] = {"carbon": random.randint(10, 20), "pep": False}

    # ----- Queries -----
    def triad_cells(self, center, orient):
        return rotate_triad(center, orient)

    def cells_free_for_triad(self, cells):
        # All within region & not occupied by peptide (any amph is removable)
        for c in cells:
            if c not in self.amph:  # outside current region
                return False
            if self.amph[c]["pep"]:
                return False
        return True

    def random_free_anchor(self, max_tries=1000):
        keys = list(self.amph.keys())
        for _ in range(max_tries):
            anchor = random.choice(keys)
            orient = random.randrange(6)
            cells = self.triad_cells(anchor, orient)
            if self.cells_free_for_triad(cells):
                return anchor, orient, cells
        return None, None, None

    # ----- Events -----
    def event_amph_swap(self):
        """Pick random amph cell, swap carbon count with a random neighbor amph (liquid-like)."""
        cells = list(self.amph.keys())
        if not cells:
            return
        a = random.choice(cells)
        # pick neighbor that exists in amph
        q, r = a
        random.shuffle(AX_NEI)
        for dq, dr in AX_NEI:
            b = (q + dq, r + dr)
            if b in self.amph:
                # swap regardless of peptide overlay; amph pool is liquid
                ca, cb = self.amph[a]["carbon"], self.amph[b]["carbon"]
                self.amph[a]["carbon"], self.amph[b]["carbon"] = cb, ca
                self.event_log.append({"evt": "amph_swap", "a": a, "b": b})
                return

    def event_peptide_insert(self):
        """Insert a triad, displacing amphiphiles (count toward growth)."""
        if len(self.amph) == 0:
            return
        anchor, orient, cells = self.random_free_anchor()
        if anchor is None:
            return
        # mark amph displaced
        displaced = 0
        for c in cells:
            if c in self.amph and not self.amph[c]["pep"]:
                displaced += 1
        self.displaced_amph += displaced

        # occupy cells by peptide
        pid = self.pid_counter
        self.pid_counter += 1
        rid = self.raft_counter
        self.raft_counter += 1
        self.peptides[pid] = {
            "cent": anchor,
            "orient": orient,
            "raft": rid,
            "inside": bool(random.getrandbits(1)),
        }
        for c in cells:
            # leave amph entry (carbon remains) but set pep True
            self.amph[c]["pep"] = True

        self.event_log.append(
            {"evt": "pept_insert", "pid": pid, "cells": cells, "displaced_inc": displaced}
        )

        # growth trigger: if displaced >= 4n+2 -> expand to (n+2)
        needed = 4 * self.n + 2
        if self.displaced_amph >= needed:
            self._grow_membrane()
            self.displaced_amph = 0
            self.event_log.append({"evt": "grow", "new_n": self.n})

    def _grow_membrane(self):
        """Expand region by +2 (n -> n+2), adding an extra ring of amphiphiles."""
        new_n = self.n + 2
        for q in range(-new_n, new_n + 1):
            for r in range(-new_n, new_n + 1):
                if self._axial_in_rhombus(q, r) and (q, r) not in self.amph:
                    self.amph[(q, r)] = {"carbon": random.randint(10, 20), "pep": False}
        self.n = new_n

    def event_peptide_diffuse(self):
        """Attempt to move one random peptide by ± one axial step if triad cells are free."""
        if not self.peptides:
            return
        pid = random.choice(list(self.peptides.keys()))
        p = self.peptides[pid]
        q, r = p["cent"]
        dq, dr = random.choice(AX_NEI)
        new_cent = (q + dq, r + dr)
        new_cells = self.triad_cells(new_cent, p["orient"])
        if not self.cells_free_for_triad(new_cells):
            return
        # free old cells
        old_cells = self.triad_cells(p["cent"], p["orient"])
        for c in old_cells:
            if c in self.amph:
                self.amph[c]["pep"] = False
        # occupy new
        for c in new_cells:
            self.amph[c]["pep"] = True
        p["cent"] = new_cent
        self.event_log.append({"evt": "pept_diffuse", "pid": pid, "from": (q, r), "to": new_cent})

    def event_raft_adhere_merge(self):
        """If two rafts touch (share edge among any triad hexes), merge smaller into larger raft_id."""
        if len(self.peptides) < 2:
            return
        # Build map: cell -> pid
        cell_to_pid = {}
        for pid, p in self.peptides.items():
            for c in self.triad_cells(p["cent"], p["orient"]):
                cell_to_pid[c] = pid

        # Find touching raft pairs via neighbors
        pairs = set()
        for c, pid in cell_to_pid.items():
            q, r = c
            for dq, dr in AX_NEI:
                nb = (q + dq, r + dr)
                if nb in cell_to_pid:
                    pid2 = cell_to_pid[nb]
                    if pid2 != pid:
                        r1 = self.peptides[pid]["raft"]
                        r2 = self.peptides[pid2]["raft"]
                        if r1 != r2:
                            pairs.add(tuple(sorted((r1, r2))))
        if not pairs:
            return
        # merge one random pair
        rA, rB = random.choice(list(pairs))
        # choose target raft (smaller id wins for determinism)
        target, source = (rA, rB) if rA < rB else (rB, rA)
        for pid, p in self.peptides.items():
            if p["raft"] == source:
                p["raft"] = target
        self.event_log.append({"evt": "raft_merge", "from": source, "into": target})

    def event_pept_flip(self):
        """Toggle inside/outside property (no geometry change, metadata only)."""
        if not self.peptides:
            return
        pid = random.choice(list(self.peptides.keys()))
        p = self.peptides[pid]
        p["inside"] = not p["inside"]
        self.event_log.append({"evt": "pept_flip", "pid": pid, "inside": p["inside"]})

    # -----------------------------
    # Drawing / Saving
    # -----------------------------
    def _raft_colors(self):
        raft_ids = sorted({p["raft"] for p in self.peptides.values()}) if self.peptides else []
        palette = get_palette(max(1, len(raft_ids)), name="tab20")
        mapping = {rid: palette[i] for i, rid in enumerate(raft_ids)}
        # default color if no peptides
        return mapping

    def save_frame(self, idx):
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_aspect("equal")
        ax.axis("off")

        # Amphiphile background (light gray tint by carbon length)
        minC, maxC = 8.0, 22.0
        for (q, r), info in self.amph.items():
            x, y = axial_to_xy(q, r, HEX_RADIUS)
            # lightness by carbon count (longer = darker)
            t = (info["carbon"] - minC) / (maxC - minC)
            t = min(max(t, 0.0), 1.0)
            base = 0.92 - 0.25 * t  # light gray → medium
            color = (base, base, base)
            patch = RegularPolygon(
                (x, y),
                numVertices=6,
                radius=HEX_RADIUS,
                orientation=math.radians(30),
                facecolor=color,
                edgecolor="k",
                linewidth=0.5,
            )
            ax.add_patch(patch)
            if LABEL_CARBONS:
                ax.text(x, y, str(info["carbon"]), ha="center", va="center", fontsize=6)

        # Peptides (triads) overlaid, colored by raft id; inside/outside shown by edge width
        raft_color = self._raft_colors()
        for pid, p in self.peptides.items():
            cells = self.triad_cells(p["cent"], p["orient"])
            fc = raft_color.get(p["raft"], (0.8, 0.2, 0.2))
            lw = 1.5 if p["inside"] else 0.6
            for (q, r) in cells:
                x, y = axial_to_xy(q, r, HEX_RADIUS)
                patch = RegularPolygon(
                    (x, y),
                    numVertices=6,
                    radius=HEX_RADIUS,
                    orientation=math.radians(30),
                    facecolor=fc,
                    edgecolor="black",
                    linewidth=lw,
                )
                ax.add_patch(patch)

        # Axis limits
        # q,r in [-n..n], use a bit of margin
        qmin, qmax = -self.n, self.n
        rmin, rmax = -self.n, self.n
        # extremes in pixel space
        corners = [
            axial_to_xy(qmin, rmin),
            axial_to_xy(qmin, rmax),
            axial_to_xy(qmax, rmin),
            axial_to_xy(qmax, rmax),
            axial_to_xy(0, -self.n),
            axial_to_xy(0, self.n),
            axial_to_xy(-self.n, 0),
            axial_to_xy(self.n, 0),
        ]
        xs = [c[0] for c in corners]
        ys = [c[1] for c in corners]
        xmin, xmax = min(xs) - BORDER, max(xs) + BORDER
        ymin, ymax = min(ys) - BORDER, max(ys) + BORDER
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        fig.savefig(os.path.join(OUT_DIR, f"frame_{idx:04d}.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)


# -----------------------------
# Case loop
# -----------------------------
def pick_event():
    # weighted random event name
    items = list(EVENT_WEIGHTS.items())
    names = [k for k, _ in items]
    w = np.array([v for _, v in items], dtype=float)
    w = w / w.sum()
    return np.random.choice(names, p=w)


def run_sim():
    mem = Membrane(n=N0)
    mem.save_frame(0)

    for step in range(1, TOTAL_EVENTS + 1):
        evt = pick_event()
        if evt == "amph_swap":
            mem.event_amph_swap()
        elif evt == "pept_insert":
            mem.event_peptide_insert()
        elif evt == "pept_diffuse":
            mem.event_peptide_diffuse()
        elif evt == "raft_adhere_merge":
            mem.event_raft_adhere_merge()
        elif evt == "pept_flip":
            mem.event_pept_flip()

        if step % SAVE_EVERY == 0:
            mem.save_frame(step)

    # final save + event log
    mem.save_frame(TOTAL_EVENTS)
    with open(os.path.join(OUT_DIR, "event_log.json"), "w") as f:
        json.dump(mem.event_log, f, indent=2)


if __name__ == "__main__":
    run_sim()
