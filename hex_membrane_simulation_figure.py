import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np

# Hex size and padding
size = 1.0
pad = 0.5

# Color palette for each triad
colors = [
    "#ffffb2", "#ff8800", "#ffb2b2", "#ff0000",
    "#b2ffb2", "#9999ff", "#cfc7b7", "#b2d1ff",
    "#ffd1dc", "#aaffaa", "#ffcc99", "#e6e6fa"
]

# Triangular triads (each 3 hexes share a common vertex)
# Triads defined using axial coordinates (q, r)
def make_triangle_triad(q, r, direction):
    # Triangle options in axial: Upward or downward triangle formation
    if direction == 'up':
        return [(q, r), (q+1, r), (q, r+1)]
    elif direction == 'down':
        return [(q, r), (q-1, r), (q, r-1)]

# Build multiple adjacent triads in triangle formation
hexes = []
directions = ['up', 'down']
triad_origins = [
    (0, 0), (2, 0), (4, 0),
    (1, 2), (3, 2), (5, 2),
    (0, 4), (2, 4), (4, 4),
    (1, 6), (3, 6), (5, 6),
]

for i, (q, r) in enumerate(triad_origins):
    color = colors[i % len(colors)]
    direction = directions[i % 2]
    triad = make_triangle_triad(q, r, direction)
    for hex_coord in triad:
        hexes.append((hex_coord, color))

# Convert axial to x, y for flat-topped hexes
def axial_to_cart(q, r, size):
    x = size * 3/2 * q
    y = size * np.sqrt(3) * (r + q / 2)
    return x, y

# Plotting
fig, ax = plt.subplots(figsize=(8, 8))
coords = []

for (q, r), col in hexes:
    x, y = axial_to_cart(q, r, size)
    coords.append((x, y))
    hex_patch = RegularPolygon(
        (x, y), numVertices=6, radius=size,
        orientation=np.radians(30),
        facecolor=col, edgecolor='k'
    )
    ax.add_patch(hex_patch)

coords = np.array(coords)
xmin, ymin = coords.min(axis=0) - pad
xmax, ymax = coords.max(axis=0) + pad

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_aspect('equal')
ax.axis('off')
plt.tight_layout()
plt.show()

