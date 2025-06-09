import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random

GRID_WIDTH = 12
GRID_HEIGHT = 8
MAX_RAFTS = 20

TRIANGLE_PATTERNS = [
    [(0, 0), (0, 1), (1, 0)],
    [(0, 0), (0, -1), (1, 0)],
    [(0, 0), (1, 0), (1, 1)],
    [(0, 0), (-1, 0), (-1, 1)],
    [(0, 0), (0, 1), (-1, 0)],
    [(0, 0), (1, 0), (0, -1)],
]

grid = np.zeros((GRID_HEIGHT, GRID_WIDTH), dtype=int)

rafts = []  # Each raft is a dict: {"id": int, "tiles": [(r, c)], "chirality": "L"/"D"}

def is_valid(triangle):
    return all(
        0 <= r < GRID_HEIGHT and 0 <= c < GRID_WIDTH and grid[r, c] == 0
        for r, c in triangle
    )

def get_neighbors(cell):
    r, c = cell
    offsets = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0)]
    neighbors = [(r + dr, c + dc) for dr, dc in offsets]
    return [(nr, nc) for nr, nc in neighbors if 0 <= nr < GRID_HEIGHT and 0 <= nc < GRID_WIDTH]

def place_random_raft(raft_id):
    for _ in range(100):
        base_r = random.randint(1, GRID_HEIGHT - 2)
        base_c = random.randint(1, GRID_WIDTH - 2)
        pattern = random.choice(TRIANGLE_PATTERNS)
        triangle = [(base_r + dr, base_c + dc) for dr, dc in pattern]
        if is_valid(triangle):
            for r, c in triangle:
                grid[r, c] = raft_id
            chirality = random.choice(["L", "D"])
            rafts.append({"id": raft_id, "tiles": triangle, "chirality": chirality})
            return True
    return False

def should_desorb(raft):
    neighbor_count = 0
    for r, c in raft["tiles"]:
        for nr, nc in get_neighbors((r, c)):
            if grid[nr, nc] != 0 and (nr, nc) not in raft["tiles"]:
                neighbor_count += 1
    return neighbor_count == 0

def move_raft(raft):
    raft_id = raft["id"]
    dr, dc = random.choice([(-1,0), (1,0), (0,-1), (0,1)])
    new_positions = [(r+dr, c+dc) for r, c in raft["tiles"]]
    if all(
        0 <= r < GRID_HEIGHT and 0 <= c < GRID_WIDTH and (grid[r, c] == 0 or (r, c) in raft["tiles"])
        for r, c in new_positions
    ):
        for r, c in raft["tiles"]:
            grid[r, c] = 0
        for r, c in new_positions:
            grid[r, c] = raft_id
        raft["tiles"] = new_positions
    return raft

def simulate_step():
    global rafts
    updated_rafts = []
    for raft in rafts:
        if should_desorb(raft):
            for r, c in raft["tiles"]:
                grid[r, c] = 0
        else:
            updated_raft = move_raft(raft)
            updated_rafts.append(updated_raft)
    rafts = updated_rafts
    if len(rafts) < MAX_RAFTS:
        place_random_raft(len(rafts) + 1)

def draw_grid():
    fig, ax = plt.subplots(figsize=(10, 6))
    dx = 3 ** 0.5
    color_map = {"L": "#1f77b4", "D": "#d62728"}  # blue vs red

    for r in range(GRID_HEIGHT):
        for c in range(GRID_WIDTH):
            val = grid[r, c]
            x = dx * c + (dx / 2 if r % 2 else 0)
            y = 1.5 * r
            face_color = "lightgrey"
            if val != 0:
                for raft in rafts:
                    if raft["id"] == val:
                        face_color = color_map[raft["chirality"]]
            hexagon = patches.RegularPolygon((x, y), numVertices=6, radius=1,
                                             orientation=np.radians(30),
                                             edgecolor='k', facecolor=face_color)
            ax.add_patch(hexagon)

    ax.set_xlim(-1, dx * (GRID_WIDTH + 1))
    ax.set_ylim(-1, 1.5 * (GRID_HEIGHT + 1))
    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig("raft_hex_chiral.png", dpi=200)
    plt.close()

# --------------------- Run Simulation -----------------------
if __name__ == "__main__":
    for _ in range(30):
        simulate_step()
    draw_grid()

