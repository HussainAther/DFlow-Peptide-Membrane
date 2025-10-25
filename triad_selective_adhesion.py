import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import os
import random

# Parameters
GRID_RADIUS = 15
HEX_SIZE = 1
N_TRIADS = 30
N_STEPS = 10
OUTPUT_FOLDER = "adh_sel"

# Canonical triad offset pattern (chevron)
TRIAD_OFFSETS = [(0, 0), (1, 0), (0, 1)]

# Utility: axial to pixel
def axial_to_pixel(q, r):
    x = HEX_SIZE * 3/2 * q
    y = HEX_SIZE * np.sqrt(3) * (r + q/2)
    return x, y

# Generate distinct colors
def get_colors(n):
    cmap = plt.colormaps['tab20']
    return [cmap(i % 20) for i in range(n)]

# Initialize triads randomly without overlap
def init_triads():
    triads = []
    occupied = set()
    colors = get_colors(N_TRIADS)

    attempts = 0
    while len(triads) < N_TRIADS and attempts < 500:
        base_q = random.randint(-GRID_RADIUS, GRID_RADIUS)
        base_r = random.randint(-GRID_RADIUS, GRID_RADIUS)

        cells = [(base_q + dq, base_r + dr) for dq, dr in TRIAD_OFFSETS]

        if any((q, r) in occupied for q, r in cells):
            attempts += 1
            continue

        for q, r in cells:
            occupied.add((q, r))

        triads.append({
            'id': f"P{len(triads)}",
            'cells': cells,
            'color': colors[len(triads)],
            'stuck': False
        })

    return triads, occupied

# Try to move free triads randomly
def step_triads(triads, occupied):
    directions = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]

    for triad in triads:
        if triad['stuck']:
            continue

        dq, dr = random.choice(directions)
        new_cells = [(q + dq, r + dr) for q, r in triad['cells']]

        # Check for collision
        if any((q, r) in occupied for q, r in new_cells):
            triad['stuck'] = True
            continue

        # Update occupied set
        for q, r in triad['cells']:
            occupied.discard((q, r))
        for q, r in new_cells:
            occupied.add((q, r))

        triad['cells'] = new_cells

# Save frame
def save_frame(triads, t):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    ax.axis('off')

    for triad in triads:
        for q, r in triad['cells']:
            x, y = axial_to_pixel(q, r)
            hex_patch = RegularPolygon(
                (x, y), numVertices=6, radius=HEX_SIZE,
                orientation=np.radians(30), facecolor=triad['color'],
                edgecolor='black'
            )
            ax.add_patch(hex_patch)

        # Add peptide label at centroid
        cx = np.mean([axial_to_pixel(q, r)[0] for q, r in triad['cells']])
        cy = np.mean([axial_to_pixel(q, r)[1] for q, r in triad['cells']])
        ax.text(cx, cy, triad['id'], ha='center', va='center', fontsize=8)

    ax.relim()
    ax.autoscale_view()
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
    fig.savefig(f"{OUTPUT_FOLDER}/adh_sel_{t:02d}.png")
    plt.close(fig)

# Main loop
def main():
    triads, occupied = init_triads()
    for t in range(N_STEPS):
        step_triads(triads, occupied)
        save_frame(triads, t)

if __name__ == "__main__":
    main()

