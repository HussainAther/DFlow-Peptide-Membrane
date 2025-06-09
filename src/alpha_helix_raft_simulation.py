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
class RaftSimulation:
    def __init__(self, grid_height=8, grid_width=12, max_rafts=20):
        self.GRID_HEIGHT = grid_height
        self.GRID_WIDTH = grid_width
        self.MAX_RAFTS = max_rafts
        self.grid = np.zeros((self.GRID_HEIGHT, self.GRID_WIDTH), dtype=int)
        self.rafts = []
        self.patterns = TRIANGLE_PATTERNS

    def is_valid(self, triangle):
        return all(
            0 <= r < self.GRID_HEIGHT and 0 <= c < self.GRID_WIDTH and self.grid[r, c] == 0
            for r, c in triangle
        )

    def get_neighbors(self, cell):
        r, c = cell
        offsets = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0)]
        neighbors = [(r + dr, c + dc) for dr, dc in offsets]
        return [(nr, nc) for nr, nc in neighbors if 0 <= nr < self.GRID_HEIGHT and 0 <= nc < self.GRID_WIDTH]
 
    def place_peptides(self):
        """
        Attempt to place up to MAX_RAFTS peptide rafts on the hex grid,
       with each raft composed of 3 connected hexagons in one of the allowed triangular configurations.
        """
        raft_id = 1
        attempts = 0
        max_attempts = 500  # Prevent infinite loops

        while raft_id <= MAX_RAFTS and attempts < max_attempts:
            base_r = random.randint(1, GRID_HEIGHT - 2)
            base_c = random.randint(1, GRID_WIDTH - 2)
            pattern = random.choice(TRIANGLE_PATTERNS)
            triangle = [(base_r + dr, base_c + dc) for dr, dc in pattern]
 
            if self.is_valid(triangle):
                chirality = random.choice(["L", "D"])
                for r, c in triangle:
                    self.grid[r, c] = raft_id
                self.rafts.append({
                    "id": raft_id,
                    "tiles": triangle,
                    "chirality": chirality
                 })
                raft_id += 1
            attempts += 1

    def place_random_raft(self, raft_id=None):
        for _ in range(100):
            base_r = random.randint(1, self.GRID_HEIGHT - 2)
            base_c = random.randint(1, self.GRID_WIDTH - 2)
            pattern = random.choice(self.patterns)
            triangle = [(base_r + dr, base_c + dc) for dr, dc in pattern]
            if self.is_valid(triangle):
                for r, c in triangle:
                    self.grid[r, c] = raft_id or len(self.rafts) + 1
                chirality = random.choice(["L", "D"])
                self.rafts.append({
                    "id": raft_id or len(self.rafts) + 1,
                    "tiles": triangle,
                    "chirality": chirality
                })
                return True
        return False

    def should_desorb(self, raft):
        neighbor_count = 0
        for r, c in raft["tiles"]:
            for nr, nc in self.get_neighbors((r, c)):
                if self.grid[nr, nc] != 0 and (nr, nc) not in raft["tiles"]:
                    neighbor_count += 1
        return neighbor_count == 0

    def move_raft(self, raft):
        raft_id = raft["id"]
        dr, dc = random.choice([(-1, 0), (1, 0), (0, -1), (0, 1)])
        new_positions = [(r + dr, c + dc) for r, c in raft["tiles"]]
        if all(
            0 <= r < self.GRID_HEIGHT and 0 <= c < self.GRID_WIDTH and
            (self.grid[r, c] == 0 or (r, c) in raft["tiles"])
            for r, c in new_positions
        ):
            for r, c in raft["tiles"]:
                self.grid[r, c] = 0
            for r, c in new_positions:
                self.grid[r, c] = raft_id
            raft["tiles"] = new_positions
        return raft

    def simulate_step(self):
        updated_rafts = []
        for raft in self.rafts:
            if self.should_desorb(raft):
                for r, c in raft["tiles"]:
                    self.grid[r, c] = 0
            else:
                updated_raft = self.move_raft(raft)
                updated_rafts.append(updated_raft)
        self.rafts = updated_rafts
        while len(self.rafts) < self.MAX_RAFTS:
            self.place_random_raft()

    def count_rafts(self):
        return len(self.rafts)

    def chirality_fraction(self, chirality_type):
        if not self.rafts:
            return 0.0
        return sum(1 for raft in self.rafts if raft["chirality"] == chirality_type) / len(self.rafts)

    def save_raft_distribution(self, filename="raft_hex_chiral.png"):
        fig, ax = plt.subplots(figsize=(10, 6))
        dx = 3 ** 0.5
        color_map = {"L": "#1f77b4", "D": "#d62728"}

        for r in range(self.GRID_HEIGHT):
            for c in range(self.GRID_WIDTH):
                val = self.grid[r, c]
                x = dx * c + (dx / 2 if r % 2 else 0)
                y = 1.5 * r
                face_color = "lightgrey"
                if val != 0:
                    for raft in self.rafts:
                        if raft["id"] == val:
                            face_color = color_map[raft["chirality"]]
                hexagon = patches.RegularPolygon((x, y), numVertices=6, radius=1,
                                                 orientation=np.radians(30),
                                                 edgecolor='k', facecolor=face_color)
                ax.add_patch(hexagon)

        ax.set_xlim(-1, dx * (self.GRID_WIDTH + 1))
        ax.set_ylim(-1, 1.5 * (self.GRID_HEIGHT + 1))
        ax.set_aspect('equal')
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(filename, dpi=200)
        plt.close()
