import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

GRID_WIDTH = 12
GRID_HEIGHT = 10
NUM_LIPIDS = 50
NUM_RAFTS = 10

# Triangular raft shapes (three connected hexes)
TRIANGLE_PATTERNS = [
    [(0, 0), (0, 1), (1, 0)],
    [(0, 0), (0, -1), (1, 0)],
    [(0, 0), (1, 0), (1, 1)],
    [(0, 0), (-1, 0), (-1, 1)],
    [(0, 0), (0, 1), (-1, 0)],
    [(0, 0), (1, 0), (0, -1)],
]

def is_valid(triangle, grid):
    return all(
        0 <= r < GRID_HEIGHT and 0 <= c < GRID_WIDTH and grid[r, c] == 0
        for r, c in triangle
    )

def place_lipids(grid, count):
    placed = 0
    while placed < count:
        r = random.randint(0, GRID_HEIGHT - 1)
        c = random.randint(0, GRID_WIDTH - 1)
        if grid[r, c] == 0:
            grid[r, c] = 1  # lipid marker
            placed += 1

def place_rafts(grid, count, raft_list):
    placed = 0
    while placed < count:
        base_r = random.randint(1, GRID_HEIGHT - 2)
        base_c = random.randint(1, GRID_WIDTH - 2)
        pattern = random.choice(TRIANGLE_PATTERNS)
        triangle = [(base_r + dr, base_c + dc) for dr, dc in pattern]
        if is_valid(triangle, grid):
            for r, c in triangle:
                grid[r, c] = 3  # alpha helix marker
            raft_list.append(triangle)
            placed += 1

def draw_hex_grid(grid, raft_list, filename="hex_membrane_13_ratio.png"):
    fig, ax = plt.subplots(figsize=(10, 8))
    dx = np.sqrt(3)
    color_map = {0: "white", 1: "green", 3: "blue"}

    for r in range(GRID_HEIGHT):
        for c in range(GRID_WIDTH):
            val = grid[r, c]
            x = dx * c + (dx / 2 if r % 2 else 0)
            y = 1.5 * r
            hexagon = patches.RegularPolygon(
                (x, y), numVertices=6, radius=1,
                orientation=np.radians(30),
                edgecolor='k', facecolor=color_map.get(val, "grey")
            )
            ax.add_patch(hexagon)

    ax.set_xlim(-1, dx * (GRID_WIDTH + 1))
    ax.set_ylim(-1, 1.5 * (GRID_HEIGHT + 1))
    ax.set_aspect('equal')
    ax.axis('off')
    plt.title("Membrane with 1-unit lipids and 3-unit alpha-helix rafts")
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.close()
    print(f"Saved: {filename}")

if __name__ == "__main__":
    grid = np.zeros((GRID_HEIGHT, GRID_WIDTH), dtype=int)
    rafts = []

    place_lipids(grid, NUM_LIPIDS)
    place_rafts(grid, NUM_RAFTS, rafts)
    draw_hex_grid(grid, rafts)

