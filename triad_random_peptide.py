"""
triad_random_motion.py
----------------------
Minimal demo:  random Brownian-like motion of peptide-triad
clusters on a flat-top hexagonal grid.

Updated: Saves all PNG frames into a subfolder ./adh_sel/
"""

import os
import random
import math
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

# ─────────────────────────────── SETTINGS ─────────────────────────────────────
HEX_SIZE = 1.0
N_FRAMES = 10
OUTDIR = "adh_sel"
os.makedirs(OUTDIR, exist_ok=True)

AXIAL_NEIGHBOURS = [(+1, 0), (+1, -1), (0, -1),
                    (-1, 0), (-1, +1), (0, +1)]

def axial_to_cart(q, r, s=HEX_SIZE):
    x = 1.5 * q * s
    y = math.sqrt(3) * (r + q/2) * s
    return x, y


# ────────────────────────────  TRIAD DATA STRUCTURE  ───────────────────────────
class Triad:
    RELATIVE_HEXES = [(0,0), (1,0), (0,1)]

    def __init__(self, q, r, colour):
        self.q, self.r = q, r
        self.colour = colour

    def hexes(self):
        return [(self.q+dq, self.r+dr)
                for dq, dr in self.RELATIVE_HEXES]

    def try_random_move(self, occupied):
        random_dirs = random.sample(AXIAL_NEIGHBOURS, len(AXIAL_NEIGHBOURS))
        for dq, dr in random_dirs:
            new_q = self.q + dq
            new_r = self.r + dr
            new_cells = [(new_q+ddq, new_r+ddr)
                         for ddq, ddr in self.RELATIVE_HEXES]
            if all(cell not in occupied for cell in new_cells):
                self.q, self.r = new_q, new_r
                return


# ────────────────────────────  INITIAL SET-UP  ────────────────────────────────
COLOURS = ['tab:blue','tab:orange','tab:green','tab:red',
           'tab:purple','tab:brown','tab:pink','tab:gray',
           'tab:olive','tab:cyan','gold','lime','violet','coral',
           'navy','firebrick','darkgreen','slategray','peru','teal']

random.seed(0)

triads = []
for i in range(25):
    q = random.randint(-7, 7)
    r = random.randint(-7, 7)
    triads.append(Triad(q, r, COLOURS[i % len(COLOURS)]))


# ───────────────────────────────  MAIN LOOP  ───────────────────────────────────
for frame in range(N_FRAMES):
    occupied = set()
    for t in triads:
        occupied.update(t.hexes())

    for t in triads:
        for cell in t.hexes():
            occupied.discard(cell)
        t.try_random_move(occupied)
        occupied.update(t.hexes())

    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_aspect('equal')
    ax.axis('off')

    for t in triads:
        for (q,r) in t.hexes():
            cx, cy = axial_to_cart(q, r)
            hex_patch = RegularPolygon(
                (cx, cy), numVertices=6, radius=HEX_SIZE*0.95,
                orientation=0, edgecolor='k', facecolor=t.colour, lw=0.8)
            ax.add_patch(hex_patch)

    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    plt.tight_layout()
    outpath = os.path.join(OUTDIR, f'adh_sel_{frame:02d}.png')
    plt.savefig(outpath, dpi=150)
    plt.close()

print(f"{N_FRAMES} frames saved to folder: {OUTDIR}/")

