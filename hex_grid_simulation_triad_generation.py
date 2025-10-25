import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import random
import os
import json

# Hex grid configuration
HEX_SIZE = 1.0
GRID_WIDTH = 20
GRID_HEIGHT = 20
N_STEPS = 10
SAVE_DIR = "hex_steps"
os.makedirs(SAVE_DIR, exist_ok=True)

# Colors
colors = {
    "yellow": "#ffffb2",
    "orange": "#ff8800",
    "green":  "#b2ffb2",
    "pink":   "#ffb2b2",
    "red":    "#ff0000",
    "blue":   "#9999ff",
    "gray":   "#cfc7b7"
}

# Triad definitions (3 adjacent hexes forming a triangle)
triad_offsets = [(0, 0), (1, 0), (0, 1)]  # Q, R axial

# Convert axial coordinates to cartesian for flat-top hexes
def axial_to_cart(q, r):
    x = HEX_SIZE * 3/2 * q
    y = HEX_SIZE * np.sqrt(3) * (r + q / 2)
    return x, y

# Generate initial triads
def generate_initial_triads():
    triads = []
    used_positions = set()
    for color in list(colors.keys()):
        placed = False
        while not placed:
            base_q = random.randint(0, GRID_WIDTH - 2)
            base_r = random.randint(0, GRID_HEIGHT - 2)
            triad = [(base_q + dq, base_r + dr) for dq, dr in triad_offsets]
            if all((q, r) not in used_positions for q, r in triad):
                triads.append({"hexes": triad, "color": color})
                used_positions.update(triad)
                placed = True
    return triads

# Move triads randomly

def move_triads(triads):
    directions = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    new_triads = []
    for triad in triads:
        dq, dr = random.choice(directions)
        new_hexes = [(q + dq, r + dr) for q, r in triad["hexes"]]
        new_triads.append({"hexes": new_hexes, "color": triad["color"]})
    return new_triads

# Divide a random triad

def divide_triads(triads):
    if len(triads) < 15 and random.random() < 0.3:
        parent = random.choice(triads)
        offset = random.choice([(2, 0), (0, 2), (-2, 0), (0, -2)])
        new_hexes = [(q + offset[0], r + offset[1]) for q, r in parent["hexes"]]
        triads.append({"hexes": new_hexes, "color": parent["color"]})
    return triads

# Plot and save

def plot_triads(triads, step):
    fig, ax = plt.subplots(figsize=(6, 6))
    for triad in triads:
        for q, r in triad["hexes"]:
            x, y = axial_to_cart(q, r)
            hex_patch = RegularPolygon(
                (x, y), numVertices=6, radius=HEX_SIZE,
                orientation=np.radians(30),
                facecolor=colors[triad["color"]], edgecolor='k'
            )
            ax.add_patch(hex_patch)
    ax.set_xlim(-5, 25)
    ax.set_ylim(-5, 25)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.savefig(f"{SAVE_DIR}/step_{step}.png")
    plt.close(fig)

# Save triad metadata per step
def save_metadata(triads, step):
    with open(f"{SAVE_DIR}/step_{step}.json", "w") as f:
        json.dump(triads, f)

# Simulation loop
def simulate():
    triads = generate_initial_triads()
    for step in range(N_STEPS):
        plot_triads(triads, step)
        save_metadata(triads, step)
        triads = move_triads(triads)
        triads = divide_triads(triads)

simulate()

