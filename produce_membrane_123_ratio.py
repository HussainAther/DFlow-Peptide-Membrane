import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

# Hex grid dimensions
GRID_WIDTH = 12
GRID_HEIGHT = 8

# Define types with size ratios and colors
MOLECULE_TYPES = {
    "bacteria": {"size": 1, "color": "green"},
    "archaea": {"size": 2, "color": "orange"},
    "helix": {"size": 3, "color": "blue"}
}

def assign_molecules(grid_height, grid_width):
    grid = np.full((grid_height, grid_width), None)
    choices = list(MOLECULE_TYPES.keys())
    weights = [1, 1, 1]  # Equal weights initially; tweak for realism

    for r in range(grid_height):
        for c in range(grid_width):
            molecule = random.choices(choices, weights)[0]
            grid[r, c] = molecule

    return grid

def draw_hex_grid(grid):
    dx = np.sqrt(3)
    fig, ax = plt.subplots(figsize=(10, 6))

    for r in range(grid.shape[0]):
        for c in range(grid.shape[1]):
            molecule = grid[r, c]
            props = MOLECULE_TYPES[molecule]
            radius = 0.4 + 0.15 * (props["size"] - 1)  # size scaling
            x = dx * c + (dx / 2 if r % 2 else 0)
            y = 1.5 * r
            hexagon = patches.RegularPolygon(
                (x, y),
                numVertices=6,
                radius=radius,
                orientation=np.radians(30),
                edgecolor='k',
                facecolor=props["color"]
            )
            ax.add_patch(hexagon)

    ax.set_xlim(-1, dx * (grid.shape[1] + 1))
    ax.set_ylim(-1, 1.5 * (grid.shape[0] + 1))
    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig("hex_membrane_123_ratio.png", dpi=200)
    plt.close()
    print("Saved: hex_membrane_123_ratio.png")

# Run simulation
grid = assign_molecules(GRID_HEIGHT, GRID_WIDTH)
draw_hex_grid(grid)

