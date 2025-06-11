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
    def __init__(self, grid_height=8, grid_width=12, max_rafts=20, enable_catalysis=True):
        self.GRID_HEIGHT = grid_height
        self.GRID_WIDTH = grid_width
        self.MAX_RAFTS = max_rafts
        self.grid = np.zeros((self.GRID_HEIGHT, self.GRID_WIDTH), dtype=int)
        self.rafts = []
        self.patterns = TRIANGLE_PATTERNS
        self.enable_catalysis = enable_catalysis

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

    def catalytic_bias(self, r, c):
        neighbors = self.get_neighbors((r, c))
        l_count = 0
        d_count = 0
        for nr, nc in neighbors:
            val = self.grid[nr, nc]
            if val != 0:
                for raft in self.rafts:
                    if raft["id"] == val:
                        if raft["chirality"] == "L":
                            l_count += 1
                        elif raft["chirality"] == "D":
                            d_count += 1
        total = l_count + d_count + 1e-5
        l_prob = l_count / total
        d_prob = d_count / total
        return l_prob, d_prob

    def place_random_raft(self, raft_id=None):
        for _ in range(100):
            base_r = random.randint(1, self.GRID_HEIGHT - 2)
            base_c = random.randint(1, self.GRID_WIDTH - 2)
            pattern = random.choice(self.patterns)
            triangle = [(base_r + dr, base_c + dc) for dr, dc in pattern]
            if self.is_valid(triangle):
                if self.enable_catalysis:
                    l_count = sum(1 for raft in self.rafts if raft["chirality"] == "L")
                    d_count = sum(1 for raft in self.rafts if raft["chirality"] == "D")
                    total = l_count + d_count
                    if total > 0:
                        l_prob = l_count / total
                        d_prob = d_count / total
                    else:
                        l_prob = d_prob = 0.5  # fallback when no rafts yet
                    chirality = random.choices(["L", "D"], weights=[l_prob, d_prob])[0]
                else:
                    chirality = random.choice(["L", "D"])
                new_id = raft_id or len(self.rafts) + 1
                for r, c in triangle:
                    self.grid[r, c] = new_id
                self.rafts.append({
                    "id": new_id,
                    "tiles": triangle,
                    "chirality": chirality
                })
                return True
        return False

    def place_peptides(self):
        raft_id = 1
        attempts = 0
        max_attempts = 500
        while raft_id <= self.MAX_RAFTS and attempts < max_attempts:
            if self.place_random_raft(raft_id):
                raft_id += 1
            attempts += 1

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

if __name__ == "__main__":
    sim = RaftSimulation(enable_catalysis=True)  # toggle True/False
    sim.place_peptides()
    for _ in range(30):
        sim.simulate_step()
    sim.save_raft_distribution("raft_hex_chiral_final.png")
    print(f"Rafts placed: {sim.count_rafts()}")
    print(f"L-fraction: {sim.chirality_fraction('L'):.2f}")
    print(f"D-fraction: {sim.chirality_fraction('D'):.2f}")

