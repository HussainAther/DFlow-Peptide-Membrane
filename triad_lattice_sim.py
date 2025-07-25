import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import matplotlib.colors as mcolors
import numpy as np
import os
import random

# Constants
TRIAD_COUNT = 12
HEX_SIZE = 1.0
SAVE_DIR = "frames"

# Grid spacing for hex layout (pointy-top)
def hex_to_xy(q, r):
    """Converts axial (q, r) coordinates to x, y in a pointy-top hex grid."""
    x = HEX_SIZE * 3/2 * q
    y = HEX_SIZE * np.sqrt(3) * (r + q/2)
    return x, y

# Triad structure: base + two neighbors forming a triangle
# Using pointy-top offsets: (0,0) base, (1,0) right neighbor, (1,-1) top-right neighbor
# This forms a triangle pointing upwards.
TRIAD_OFFSETS = [(0, 0), (1, 0), (1, -1)] # Corrected to a valid connected triad

# Neighbors for adjacency check (all 6 neighbors in axial coordinates for pointy-top)
HEX_NEIGHBORS = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]

# Colors
colors = list(mcolors.TABLEAU_COLORS.values())

class PeptideTriad:
    def __init__(self, q, r, color):
        self.q = q
        self.r = r
        self.color = color
        self.anchor = random.randint(8, 14)  # anchor_length: hydrophobic portion
        self.tail = random.randint(2, 8)     # tail_length: remainder of peptide
        self.chirality = random.choice(['L', 'D'])  # Chirality: L or D
        self.charge = random.choice([-1, 0, 1])    # Charge bias: -1, 0, +1

    def hex_cells(self):
        return [(self.q + dq, self.r + dr) for dq, dr in TRIAD_OFFSETS]

    def labels(self):
        return [f"T{self.anchor}", f"A{self.anchor}", f"C{self.tail}"]

def get_all_neighbors(q, r):
    """Get coordinates of all 6 neighbors of a hex."""
    return [(q + dq, r + dr) for dq, dr in HEX_NEIGHBORS]

def create_triads(count):
    """
    Creates triads ensuring each new triad touches at least one hex
    from a previously placed triad.
    """
    if count <= 0:
        return []

    triads = []
    placed_hex_coords = set()  # Set of all individual hex coordinates placed so far
    available_neighbor_hexes = set()  # Set of coordinates adjacent to placed hexes where a new triad base can go

    # Place the first triad at origin
    color = colors[len(triads) % len(colors)]
    first_triad = PeptideTriad(0, 0, color)
    triads.append(first_triad)
    placed_hex_coords.update(first_triad.hex_cells())
    
    # Add neighbors of the first triad's hexes to the available set
    for hq, hr in first_triad.hex_cells():
        for nq, nr in get_all_neighbors(hq, hr):
            if (nq, nr) not in placed_hex_coords:
                 available_neighbor_hexes.add((nq, nr))

    # Place subsequent triads
    for _ in range(1, count):
        max_attempts = 1000 # Prevent infinite loop
        placed = False
        attempts = 0
        while not placed and attempts < max_attempts:
            attempts += 1
            if not available_neighbor_hexes:
                print(f"Warning: Could not place triad after {attempts-1} attempts. Stopping.")
                return triads # Return what we have so far

            # Pick a random available neighbor hex to try as the base of the new triad
            base_q, base_r = random.choice(list(available_neighbor_hexes))
            
            # Create a trial triad at this base position
            trial_triad = PeptideTriad(base_q, base_r, colors[len(triads) % len(colors)])
            trial_cells = trial_triad.hex_cells()

            # Check if the entire trial triad fits without overlapping existing hexes
            if all(cell not in placed_hex_coords for cell in trial_cells):
                # Success! Place the triad
                triads.append(trial_triad)
                placed_hex_coords.update(trial_cells)
                placed = True
                
                # Remove the base cell from available neighbors (it's now occupied)
                available_neighbor_hexes.discard((base_q, base_r))
                
                # Add new neighbors of the newly placed hexes to available set
                for hq, hr in trial_cells:
                    for nq, nr in get_all_neighbors(hq, hr):
                        if (nq, nr) not in placed_hex_coords:
                            available_neighbor_hexes.add((nq, nr))
                
                # Optimization: Remove any newly occupied cells from available set
                # (In case a neighbor of the new triad was also a neighbor of an old one)
                available_neighbor_hexes.difference_update(trial_cells)

        if not placed:
             print(f"Warning: Could not place triad after {max_attempts} attempts. Stopping.")
             return triads # Return what we have so far
            
    return triads

def draw_triads(triads, filename="frame_00.png"):
    fig, ax = plt.subplots(figsize=(10, 10)) # Slightly larger figure

    for triad in triads:
        for (q, r), label in zip(triad.hex_cells(), triad.labels()):
            x, y = hex_to_xy(q, r)
            hex_patch = RegularPolygon((x, y), numVertices=6, radius=HEX_SIZE / np.sqrt(3),
                                       orientation=np.radians(30), # Pointy top
                                       facecolor=triad.color, edgecolor='black', linewidth=0.5)
            ax.add_patch(hex_patch)
            ax.text(x, y, label, ha='center', va='center', fontsize=6, color='white', weight='bold')

    ax.set_aspect('equal')
    ax.axis('off')
    # Dynamically adjust limits based on placed triads if needed, or keep fixed for consistency
    ax.set_xlim(-3, 20)
    ax.set_ylim(-10, 5)

    os.makedirs(SAVE_DIR, exist_ok=True)
    plt.savefig(os.path.join(SAVE_DIR, filename), bbox_inches='tight', dpi=150) # Higher DPI
    plt.close()

def main():
    triads = create_triads(TRIAD_COUNT)
    print(f"Placed {len(triads)} triads.")
    draw_triads(triads)

if __name__ == "__main__":
    main()
