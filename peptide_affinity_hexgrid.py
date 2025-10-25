import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import math
import os
import random

# --- CONFIGURATION ---
NUM_TRIADS = 10
NUM_STEPS = 10
SAVE_DIR = "affinity_frames"

USE_REPLICATION_BIAS = True
USE_ADHESION_BIAS = True

np.random.seed(0)
random.seed(0)

# --- DATA STRUCTURES ---
class Peptide:
    def __init__(self, id, chirality):
        self.id = id
        self.chirality = chirality
        self.affinity = set()

class Triad:
    def __init__(self, pid, center):
        self.pid = pid
        self.center = center
        self.angle = np.random.uniform(0, 2 * np.pi)
        self.positions = self.compute_positions()
        self.color = plt.cm.tab20(pid % 20)
        self.chirality = random.choice(['L', 'D'])

    def compute_positions(self):
        positions = []
        for i in range(3):
            angle = self.angle + i * 2 * np.pi / 3
            x = self.center[0] + math.cos(angle)
            y = self.center[1] + math.sin(angle)
            positions.append((x, y))
        return positions

# --- FUNCTIONS ---
def adhesion_affinity(p1, p2):
    """Return 1 if peptides can bind based on chirality rules."""
    if USE_ADHESION_BIAS and p1.chirality != p2.chirality:
        return 0
    return 1

def replicate_peptide(parent):
    """Create a new peptide from parent, inheriting chirality if bias is on."""
    chirality = parent.chirality if USE_REPLICATION_BIAS else random.choice(['L', 'D'])
    return Peptide(parent.id, chirality)

def random_position():
    return (np.random.uniform(0, 10), np.random.uniform(0, 10))

def init_world():
    triads = []
    peptides = []
    for pid in range(NUM_TRIADS):
        triad = Triad(pid, random_position())
        triads.append(triad)
        pep = Peptide(pid, triad.chirality)
        peptides.append(pep)

    # Set affinity randomly, but respect chirality constraints
    for a in peptides:
        for b in peptides:
            if a.id != b.id and adhesion_affinity(a, b):
                if random.random() < 0.5:
                    a.affinity.add(b.id)

    return triads, peptides

def save_frame(triads, step):
    fig, ax = plt.subplots(figsize=(6, 6))
    for triad in triads:
        for x, y in triad.positions:
            hex_patch = RegularPolygon((x, y), numVertices=6, radius=0.5, orientation=np.radians(30),
                                       facecolor=triad.color, edgecolor='k', lw=1)
            ax.add_patch(hex_patch)

        # Add label in the center
        ax.text(triad.center[0], triad.center[1], str(triad.pid), ha='center', va='center', fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=0.5))

    ax.set_xlim(0, 12)
    ax.set_ylim(0, 12)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout()

    os.makedirs(SAVE_DIR, exist_ok=True)
    plt.savefig(f"{SAVE_DIR}/affinity_{step:02d}.png", dpi=150)
    plt.close()

def main():
    triads, peptides = init_world()

    for step in range(NUM_STEPS):
        save_frame(triads, step)

        # Update world: triads move slightly & try to adhere to others
        for triad in triads:
            dx, dy = np.random.uniform(-0.3, 0.3), np.random.uniform(-0.3, 0.3)
            triad.center = (triad.center[0] + dx, triad.center[1] + dy)
            triad.positions = triad.compute_positions()

if __name__ == "__main__":
    main()

