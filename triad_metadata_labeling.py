import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import random
import os

# ─── Config ──────────────────────────────────────────────────────────────
HEX_SIZE = 1.0
PADDING = 0.6
N_TRIADS = 9
N_AMPHIPHILES = 12
STEPS = 10
OUTDIR = "frames_affinity"
os.makedirs(OUTDIR, exist_ok=True)

TRIAD_OFFSETS = [(0, 0), (1, 0), (0, 1)]  # Chevron shape
NEIGHBOR_DIRS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]  # Hex dirs

# ─── Axial to Cartesian ──────────────────────────────────────────────────
def axial_to_cart(q, r, size=HEX_SIZE):
    x = size * 3/2 * q
    y = size * np.sqrt(3) * (r + q / 2)
    return x, y

# ─── Movement step ───────────────────────────────────────────────────────
def random_hex_step():
    return random.choice(NEIGHBOR_DIRS)

# ─── Triad + Amphiphile Init ─────────────────────────────────────────────
def initialize_triads():
    triads = []
    used = set()
    for i in range(N_TRIADS):
        base_q = (i % 3) * 4
        base_r = (i // 3) * 4
        peptide_id = f"P{i+1:02d}"
        chirality = random.choice([-1, 1])
        hexes = []
        for dq, dr in TRIAD_OFFSETS:
            q, r = base_q + dq, base_r + dr
            used.add((q, r))
            hexes.append((q, r))
        triads.append({"id": peptide_id, "chirality": chirality, "hexes": hexes})
    return triads

def initialize_amphiphiles(triads):
    aphs = set()
    while len(aphs) < N_AMPHIPHILES:
        triad = random.choice(triads)
        q, r = random.choice(triad["hexes"])
        dq, dr = random_hex_step()
        coord = (q + dq, r + dr)
        aphs.add(coord)
    return list(aphs)

# ─── Triad Movement Toward Similar Chirality ─────────────────────────────
def biased_step(triad, neighbors):
    same_chirals = [t for t in neighbors if t["chirality"] == triad["chirality"]]
    if not same_chirals:
        return random_hex_step()
    target = random.choice(same_chirals)
    tq, tr = target["hexes"][0]
    cq, cr = triad["hexes"][0]
    dq = np.sign(tq - cq)
    dr = np.sign(tr - cr)
    return (dq, dr)

# ─── Frame Saving ────────────────────────────────────────────────────────
def save_frame(triads, amphiphiles, step):
    fig, ax = plt.subplots(figsize=(8, 8))
    coords = []

    for triad in triads:
        color = "#ccffff" if triad["chirality"] == 1 else "#ffcccc"
        for q, r in triad["hexes"]:
            x, y = axial_to_cart(q, r)
            coords.append((x, y))
            patch = RegularPolygon((x, y), numVertices=6, radius=HEX_SIZE,
                                   orientation=np.radians(30), facecolor=color, edgecolor='k')
            ax.add_patch(patch)

        # Label triad
        q_vals = [q for q, r in triad["hexes"]]
        r_vals = [r for q, r in triad["hexes"]]
        cq = sum(q_vals) / len(q_vals)
        cr = sum(r_vals) / len(r_vals)
        cx, cy = axial_to_cart(cq, cr)
        ax.text(cx, cy, triad["id"], ha='center', va='center', fontsize=7, weight='bold')

    # Amphiphiles (single-hex, dark blue)
    for q, r in amphiphiles:
        x, y = axial_to_cart(q, r)
        patch = RegularPolygon((x, y), numVertices=6, radius=HEX_SIZE,
                               orientation=np.radians(30), facecolor="#336699", edgecolor='k')
        ax.add_patch(patch)
        coords.append((x, y))

    coords = np.array(coords)
    xmin, ymin = coords.min(axis=0) - PADDING
    xmax, ymax = coords.max(axis=0) + PADDING
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, f"step_{step:02d}.png"), dpi=300)
    plt.close()

# ─── Main Loop ───────────────────────────────────────────────────────────
triads = initialize_triads()
amphiphiles = initialize_amphiphiles(triads)

for step in range(STEPS):
    save_frame(triads, amphiphiles, step)

    # Update triads
    new_triads = []
    for triad in triads:
        all_others = [t for t in triads if t != triad]
        dq, dr = biased_step(triad, all_others)
        new_hexes = [(q + dq, r + dr) for q, r in triad["hexes"]]
        triad["hexes"] = new_hexes
        new_triads.append(triad)
    triads = new_triads

    # Move amphiphiles
    amphiphiles = [(q + dq, r + dr) if random.random() < 0.5 else (q, r)
                   for q, r in amphiphiles
                   for dq, dr in [random_hex_step()]]

