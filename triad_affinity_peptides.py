import matplotlib
matplotlib.use("Agg")  # Headless backend
print("Using matplotlib backend:", matplotlib.get_backend())

import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import os
import random

# Grid and simulation parameters
GRID_WIDTH = 30
GRID_HEIGHT = 30
HEX_SIZE = 1
N_TRIADS = 20
N_STEPS = 10
OUTPUT_DIR = "adh_affinity"

# Create output folder
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────
# Peptide metadata
class Peptide:
    def __init__(self, peptide_id, color):
        self.peptide_id = peptide_id
        self.color = color
        self.chirality = random.choice([-1, 1])
        self.hydrophobicity = round(random.uniform(0, 1), 2)
        self.affinities = {}  # Filled later

# Triad data
class Triad:
    def __init__(self, center, orientation, peptide):
        self.center = center  # (i, j)
        self.orientation = orientation  # 0–5 hex directions
        self.peptide = peptide

    def get_hexes(self):
        i, j = self.center
        offsets = [
            (0, 0),
            (1, 0),
            (0, 1),
            (-1, 1),
            (-1, 0),
            (0, -1),
            (1, -1)
        ]
        orientation_offsets = [offsets[o % 6] for o in range(self.orientation, self.orientation + 3)]
        return [(i + dx, j + dy) for dx, dy in orientation_offsets]

    def move_random(self):
        dx, dy = random.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
        self.center = ((self.center[0] + dx) % GRID_WIDTH,
                       (self.center[1] + dy) % GRID_HEIGHT)
        self.orientation = (self.orientation + random.choice([-1, 0, 1])) % 6

# ─────────────────────────────────────────────────────────────
# Initialization
def init_triads():
    triads = []
    occupied = set()
    colors = plt.get_cmap('tab20').colors
    peptides = []

    for idx in range(N_TRIADS):
        while True:
            i, j = random.randint(0, GRID_WIDTH - 1), random.randint(0, GRID_HEIGHT - 1)
            temp_triad = Triad((i, j), random.randint(0, 5), None)
            positions = temp_triad.get_hexes()
            if all(p not in occupied for p in positions):
                for p in positions:
                    occupied.add(p)
                peptide = Peptide(f"P{idx:02}", colors[idx % len(colors)])
                triad = Triad((i, j), temp_triad.orientation, peptide)
                triads.append(triad)
                peptides.append(peptide)
                break

    # Build symmetric affinity matrix
    for p1 in peptides:
        for p2 in peptides:
            if p2.peptide_id not in p1.affinities:
                score = random.uniform(0, 1)
                p1.affinities[p2.peptide_id] = score
                p2.affinities[p1.peptide_id] = score

    return triads

# ─────────────────────────────────────────────────────────────
# Save frame to PNG
def save_frame(triads, step):
    print(f"Saving frame {step}...")
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    ax.axis('off')

    for triad in triads:
        for (i, j) in triad.get_hexes():
            x = i + 0.5 * (j % 2)
            y = j * np.sqrt(3) / 2
            patch = RegularPolygon(
                (x, y), numVertices=6, radius=HEX_SIZE / np.sqrt(3),
                orientation=np.radians(30),
                facecolor=triad.peptide.color, edgecolor='k', linewidth=1
            )
            ax.add_patch(patch)

        # Label peptide ID
        i, j = triad.center
        x = i + 0.5 * (j % 2)
        y = j * np.sqrt(3) / 2
        ax.text(x, y, triad.peptide.peptide_id, fontsize=6, ha='center', va='center', color='white')

    # Fix plot bounds
    ax.set_xlim(-1, GRID_WIDTH + 1)
    ax.set_ylim(-1, GRID_HEIGHT * np.sqrt(3)/2 + 1)

    fig.savefig(os.path.join(OUTPUT_DIR, f"adh_affinity_{step:02}.png"),
                bbox_inches='tight', pad_inches=0.1)
    plt.close()

# ─────────────────────────────────────────────────────────────
# Apply affinity rules
def apply_affinity_interactions(triads):
    for triad in triads:
        neighbors = [t for t in triads if t != triad]
        for other in neighbors:
            dist = np.linalg.norm(np.array(triad.center) - np.array(other.center))
            if dist < 2:
                affinity = triad.peptide.affinities[other.peptide.peptide_id]
                if random.random() < affinity:
                    dx = np.sign(other.center[0] - triad.center[0])
                    dy = np.sign(other.center[1] - triad.center[1])
                    triad.center = ((triad.center[0] + dx) % GRID_WIDTH,
                                    (triad.center[1] + dy) % GRID_HEIGHT)

# ─────────────────────────────────────────────────────────────
# Main
def main():
    triads = init_triads()
    for step in range(N_STEPS):
        for triad in triads:
            triad.move_random()
        apply_affinity_interactions(triads)
        save_frame(triads, step)

if __name__ == "__main__":
    main()

