"""
Hex-grid peptide-triad simulation
– Adds minimal “adhesion” logic.
– Saves each frame as PNG (no GUI).
"""

import matplotlib
matplotlib.use("Agg")          # headless PNG output
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
from pathlib import Path

# ── CONFIG ──────────────────────────────────────────────────────────
SIZE        = 1.0                # hex radius
N_STEPS     = 10                 # frames to save
OUTDIR      = Path("adhesion_steps")
OUTDIR.mkdir(exist_ok=True)
COLORS      = {
    "yellow": "#ffffb2", "orange": "#ff8800", "green": "#b2ffb2",
    "pink":   "#ffb2b2", "red":    "#ff0000", "blue":  "#9999ff",
    "gray":   "#cfc7b7"
}
ANCHORS = {                      # seed triads (axial q,r)
    "yellow": (0, 0), "orange": (3, 0), "green": (6, 0),
    "pink":   (1, 3), "red":    (4, 3), "blue":  (0, 6),
    "gray":   (3, 6),
}
OFFSETS = [(0, 0), (1, 0), (0, 1)]   # chevron triad layout

# ── HELPERS ─────────────────────────────────────────────────────────
def axial_to_cart(q, r, s=SIZE):
    return (s * 1.5 * q,
            s * np.sqrt(3) * (r + q/2))

def axial_neighbors(q, r):
    """Six axial neighbors for flat-top hexes."""
    for dq, dr in [(1,0), (-1,0), (0,1), (0,-1), (1,-1), (-1,1)]:
        yield (q + dq, r + dr)

# ── DATA STRUCTURES ─────────────────────────────────────────────────
hex_taken = set()      # every occupied hex coordinate
triads = []            # list of dicts: {coords, color, adhered}

def add_triad(color, anchor_qr):
    """Attempt to place a triad; return True if successful."""
    coords = [(anchor_qr[0]+dq, anchor_qr[1]+dr) for dq,dr in OFFSETS]
    if any(c in hex_taken for c in coords):
        return False                     # spot blocked
    hex_taken.update(coords)
    triads.append({"coords": coords, "color": color, "adhered": False})
    return True

# ── INITIAL SEEDING ─────────────────────────────────────────────────
for col, qr in ANCHORS.items():
    add_triad(col, qr)

# ── MAIN LOOP ───────────────────────────────────────────────────────
for step in range(N_STEPS):
    # 1) Attempt replication for each non-adhered triad
    new_triads = []
    for T in triads:
        if T["adhered"]:
            continue
        # pick a random offset direction to replicate
        dq, dr = np.random.choice([1,-1]), np.random.choice([0,1,-1])
        aq, ar = T["coords"][0]          # use triad's anchor as base
        target_anchor = (aq+dq*2, ar+dr*2)  # skip two to avoid overlap
        if add_triad(T["color"], target_anchor):
            new_triads.append(triads[-1])    # last added

    # 2) Adhesion check: if any hex of new triad touches existing triad
    for T in new_triads:
        for q,r in T["coords"]:
            if any(nbr in hex_taken and      # occupied neighbor
                   all(nbr not in T2["coords"] for T2 in [T])  # not itself
                   for nbr in axial_neighbors(q,r)):
                T["adhered"] = True
                break    # one touch is enough

    # 3) Render frame
    fig, ax = plt.subplots(figsize=(6,6))
    for T in triads:
        for q,r in T["coords"]:
            x,y = axial_to_cart(q,r)
            ax.add_patch(RegularPolygon(
                xy=(x, y),
                numVertices=6,
                radius=SIZE,
                orientation=np.radians(30),
                facecolor=COLORS[T["color"]],
                edgecolor='k'
            ))
    ax.set_aspect('equal'); ax.axis('off')
    xmin,ymin = np.min([axial_to_cart(*c) for c in hex_taken], axis=0)-SIZE*1.2
    xmax,ymax = np.max([axial_to_cart(*c) for c in hex_taken], axis=0)+SIZE*1.2
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymin,ymax)
    fig.tight_layout()
    fig.savefig(OUTDIR / f"adh_{step:02d}.png", dpi=150)
    plt.close(fig)

