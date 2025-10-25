import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend

import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import os
import random
from tqdm import tqdm

# ── Configuration ────────────────────────────────────────
GRID_W, GRID_H = 30, 30
N_TRIADS = 50
N_STEPS = 50
HEX_SIZE = 1.0
TRIAD_OFFSETS = [(0, 0), (1, 0), (0, 1)]  # Chevron
COLORS = ["#%06x" % random.randint(0, 0xFFFFFF) for _ in range(N_TRIADS)]

SAVE_DIR = "frames_random"
os.makedirs(SAVE_DIR, exist_ok=True)

# ── Utility functions ─────────────────────────────────────
def axial_to_cart(q, r, size=HEX_SIZE):
    x = size * 3/2 * q
    y = size * np.sqrt(3) * (r + q / 2)
    return x, y

def wrap_coord(q, r):
    return q % GRID_W, r % GRID_H

# ── Triad class ───────────────────────────────────────────
class Triad:
    def __init__(self, q, r, color):
        self.anchor = (q, r)
        self.color = color
        self.hexes = [wrap_coord(q + dq, r + dr) for dq, dr in TRIAD_OFFSETS]

# ── Grid + simulation state ───────────────────────────────
class HexGrid:
    def __init__(self, w, h):
        self.w = w
        self.h = h
        self.occupied = set()
        self.triads = []

    def add_random_triads(self, n):
        attempts = 0
        while len(self.triads) < n and attempts < 1000:
            q, r = random.randint(0, self.w - 1), random.randint(0, self.h - 1)
            triad = Triad(q, r, COLORS[len(self.triads)])
            if all((hq, hr) not in self.occupied for (hq, hr) in triad.hexes):
                self.triads.append(triad)
                for pos in triad.hexes:
                    self.occupied.add(pos)
            attempts += 1

    def save_frame(self, step):
        fig, ax = plt.subplots(figsize=(8, 8))
        for triad in self.triads:
            for q, r in triad.hexes:
                x, y = axial_to_cart(q, r)
                hex_patch = RegularPolygon(
                    (x, y), numVertices=6, radius=HEX_SIZE,
                    orientation=np.radians(30),
                    facecolor=triad.color, edgecolor='k'
                )
                ax.add_patch(hex_patch)

        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_xlim(-HEX_SIZE, self.w * HEX_SIZE * 2)
        ax.set_ylim(-HEX_SIZE, self.h * HEX_SIZE * 2)
        plt.tight_layout()
        fname = os.path.join(SAVE_DIR, f"frame_{step:03}.png")
        plt.savefig(fname, dpi=100)
        plt.close()

    def step(self):
        # No replication in this logic (non-replicative peptides)
        pass

# ── Main simulation ───────────────────────────────────────
def main():
    grid = HexGrid(GRID_W, GRID_H)
    grid.add_random_triads(N_TRIADS)

    for t in tqdm(range(N_STEPS), desc="Simulating"):
        grid.save_frame(t)
        grid.step()

if __name__ == "__main__":
    main()

