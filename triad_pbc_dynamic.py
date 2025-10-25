import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import random
import os

# ── Config ─────────────────────────────────────────────────────────────
HEX_SIZE = 1.0
GRID_PADDING = 1
INITIAL_SIZE = 10
MAX_STEPS = 10
TRIAD_OFFSETS = [(0, 0), (1, 0), (0, 1)]
DIRECTIONS = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]

COLOR_MAP = {
    "triad_inside": "#ff9999",
    "triad_outside": "#99ccff",
    "amphiphile": "#cccccc"
}

SAVE_DIR = "frames_pbc_dynamic"
os.makedirs(SAVE_DIR, exist_ok=True)

# ── Hex Utilities ───────────────────────────────────────────────────────
def axial_to_cart(q, r, size=HEX_SIZE):
    x = size * 3/2 * q
    y = size * np.sqrt(3) * (r + q/2)
    return x, y

def wrap(q, r, w, h):
    return q % w, r % h

# ── Grid Class ──────────────────────────────────────────────────────────
class HexGrid:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.triads = {}         # (q, r) -> {'type': str, 'label': '+' or '-'}
        self.amphiphiles = set() # set of (q, r)

    def is_occupied(self, q, r):
        return (q, r) in self.triads or (q, r) in self.amphiphiles

    def grow_if_needed(self, q, r):
        expanded = False
        while not (0 <= q < self.width and 0 <= r < self.height):
            self.width += 2
            self.height += 2
            # Shift all positions to center new grid
            self.triads = {(qq+1, rr+1): v for (qq, rr), v in self.triads.items()}
            self.amphiphiles = {(qq+1, rr+1) for (qq, rr) in self.amphiphiles}
            expanded = True
        return expanded

    def place_triad(self, q, r, label="+"):
        coords = [(q+dx, r+dy) for dx, dy in TRIAD_OFFSETS]
        for cq, cr in coords:
            self.grow_if_needed(cq, cr)
            if self.is_occupied(cq, cr):
                return False
        for cq, cr in coords:
            self.triads[(cq, cr)] = {"type": "triad", "label": label}
        return True

    def place_random_nearby(self):
        triad_positions = [k for k in self.triads]
        random.shuffle(triad_positions)
        for q, r in triad_positions:
            dq, dr = random.choice(DIRECTIONS)
            q2, r2 = q + dq, r + dr
            label = self.triads[(q, r)]["label"]
            if self.place_triad(q2, r2, label=label):
                self.place_amphiphiles_around(q2, r2)
                break

    def place_amphiphiles_around(self, q, r):
        triad_coords = [(q+dx, r+dy) for dx, dy in TRIAD_OFFSETS]
        for tq, tr in triad_coords:
            for dq, dr in DIRECTIONS:
                nq, nr = tq + dq, tr + dr
                if not self.is_occupied(nq, nr):
                    self.amphiphiles.add((nq, nr))

    def save_frame(self, step):
        fig, ax = plt.subplots(figsize=(8, 8))
        for (q, r), info in self.triads.items():
            x, y = axial_to_cart(q, r)
            color = COLOR_MAP["triad_inside"] if info["label"] == "-" else COLOR_MAP["triad_outside"]
            edge = 'k'
            lw = 2 if info["label"] == "-" else 1
            hex_patch = RegularPolygon(
                xy=(x, y), numVertices=6, radius=HEX_SIZE,
                orientation=np.radians(30),
                facecolor=color, edgecolor=edge, linewidth=lw
            )
            ax.add_patch(hex_patch)

        for q, r in self.amphiphiles:
            x, y = axial_to_cart(q, r)
            hex_patch = RegularPolygon(
                xy=(x, y), numVertices=6, radius=HEX_SIZE,
                orientation=np.radians(30),
                facecolor=COLOR_MAP["amphiphile"], edgecolor='gray', linewidth=0.5
            )
            ax.add_patch(hex_patch)

        ax.set_xlim(-HEX_SIZE, self.width * HEX_SIZE * 1.6)
        ax.set_ylim(-HEX_SIZE, self.height * HEX_SIZE * 1.6)
        ax.set_aspect('equal')
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(f"{SAVE_DIR}/frame_{step:02}.png")
        plt.close()

# ── Simulation Runner ───────────────────────────────────────────────────
def main():
    grid = HexGrid(INITIAL_SIZE, INITIAL_SIZE)
    grid.place_triad(5, 5, label='-')  # Start with inside peptide

    for t in range(MAX_STEPS):
        grid.save_frame(t)
        grid.place_random_nearby()

if __name__ == "__main__":
    main()

