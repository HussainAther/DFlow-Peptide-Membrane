import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import random
import os

# Ensure output directory
os.makedirs("division_steps", exist_ok=True)

# Hex geometry
HEX_RADIUS = 1
HEX_HEIGHT = np.sqrt(3) * HEX_RADIUS
HEX_WIDTH = 2 * HEX_RADIUS

# Triad structure offsets (relative positions for triangle of 3 hexes)
TRIAD_OFFSETS = [(0, 0), (1, 0), (0.5, np.sqrt(3)/2)]

# Grid occupied positions
occupied = set()

# Define Triad class
class Triad:
    def __init__(self, center, color):
        self.center = center  # center (q, r)
        self.color = color
        self.hexes = self.calculate_hexes(center)
    
    def calculate_hexes(self, center):
        q, r = center
        return [(q + dq, r + dr) for dq, dr in TRIAD_OFFSETS]
    
    def can_divide(self):
        # Try to place a new triad in adjacent directions
        directions = [(1.5, 0), (-1.5, 0), (0.75, 1.3), (-0.75, 1.3), (0.75, -1.3), (-0.75, -1.3)]
        random.shuffle(directions)
        for dx, dy in directions:
            new_center = (self.center[0] + dx, self.center[1] + dy)
            new_hexes = [(new_center[0] + dq, new_center[1] + dr) for dq, dr in TRIAD_OFFSETS]
            if all(h not in occupied for h in new_hexes):
                return Triad(new_center, self.color)
        return None

    def add_to_occupied(self):
        for h in self.hexes:
            occupied.add(h)

# Create initial triads
colors = ['red', 'orange', 'yellow', 'lightgreen', 'lightblue', 'violet', 'pink', 'tan']
triads = []

# Start with 5 random triads
for i in range(5):
    center = (random.uniform(0, 10), random.uniform(0, 10))
    triad = Triad(center, colors[i % len(colors)])
    if all(h not in occupied for h in triad.hexes):
        triads.append(triad)
        triad.add_to_occupied()

# Simulation loop
for step in range(10):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')
    ax.axis('off')

    # Draw each hexagon
    for triad in triads:
        for x, y in triad.hexes:
            hx = x * HEX_RADIUS * 1.5
            hy = y * HEX_HEIGHT * 0.5
            hex_patch = RegularPolygon((hx, hy), numVertices=6, radius=HEX_RADIUS * 0.95,
                                       orientation=np.radians(30), facecolor=triad.color,
                                       edgecolor='k')
            ax.add_patch(hex_patch)
    
    plt.xlim(-2, 20)
    plt.ylim(-2, 20)
    plt.tight_layout()
    plt.savefig(f"division_steps/step_{step}.png", dpi=100)
    plt.close()

    # Division: some triads try to divide
    new_triads = []
    for triad in triads:
        if random.random() < 0.3:  # division probability
            new_triad = triad.can_divide()
            if new_triad:
                new_triad.add_to_occupied()
                new_triads.append(new_triad)
    triads.extend(new_triads)

