# ── triads_wall.py ─────────────────────────────────────────────────────────────
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
from pathlib import Path

# ── 1. Parameters you may tweak ───────────────────────────────────────────────
HEX_RADIUS   = 1.0          # size of each hex
MARGIN       = 0.6          # white space around the drawing
ROWS, COLS   = 4, 5         # logical rows × triads per row
OUTPUT_FILE  = Path("triads_wall.png")

# nice, high-contrast colours
COLOURS = [
    "#cfc7b7", "#ffb2b2", "#ff0000", "#b2ffb2", "#ff8800",
    "#ffffb2", "#9999ff"
]

# ── 2. Geometry helpers ───────────────────────────────────────────────────────
OFFSETS = [(0, 0), (1, 0), (0, 1)]            # chevron / “triad” footprint

def axial_to_cart(q, r, s=HEX_RADIUS):
    """Flat-top hex axial → Cartesian"""
    x = s * 1.5 * q
    y = s * np.sqrt(3) * (r + q/2)
    return x, y

# ── 3. Build anchor coordinates on a hole-free rhombus lattice ────────────────
anchors = {}                                   # colour → (q,r)
row_q, row_r = 0, 0                           # first anchor of first row
colour_idx   = 0

for _row in range(ROWS):
    q, r = row_q, row_r
    for _col in range(COLS):
        colour = COLOURS[colour_idx % len(COLOURS)]
        anchors[(colour, colour_idx)] = (q, r)
        q += 2; r -= 1                        # stride (2,-1) seals all gaps
        colour_idx += 1
    row_q += 1; row_r += 2                    # next anchor *row*

# ── 4. Expand anchors → individual hexes, collect for plotting ────────────────
hexes = []
for (colour, _cid), (aq, ar) in anchors.items():
    for dq, dr in OFFSETS:
        q, r = aq + dq, ar + dr
        hexes.append(((q, r), colour))

# ── 5. Plot & save ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 8))
cart_coords = []

for (q, r), col in hexes:
    x, y = axial_to_cart(q, r)
    cart_coords.append((x, y))
    ax.add_patch(
        RegularPolygon(
            (x, y), numVertices=6, radius=HEX_RADIUS,
            orientation=np.radians(30),           # flat-top
            facecolor=col, edgecolor="k", lw=1
        )
    )

cart_coords = np.array(cart_coords)
xmin, ymin = cart_coords.min(axis=0) - MARGIN
xmax, ymax = cart_coords.max(axis=0) + MARGIN
ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
ax.set_aspect("equal");  ax.axis("off")
plt.tight_layout()
fig.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight")
print(f"✔ Saved to {OUTPUT_FILE.resolve()}")
# plt.show()    # uncomment if you also want a pop-up window

