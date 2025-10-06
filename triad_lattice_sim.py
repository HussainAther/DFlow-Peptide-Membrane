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
# Using pointy-top offsets: (0,0) base, (1,0) right neighbor, (0,1) bottom-left neighbor
# This forms a triangle pointing upwards, suitable for tight packing.
TRIAD_OFFSETS = [(0, 0), (1, 0), (0, 1)]

# Neighbors for adjacency check (all 6 neighbors in axial coordinates for pointy-top)
HEX_NEIGHBORS = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]

# Define how a triad can be placed relative to an existing hex
# For each neighbor direction of an existing hex, define the base (q,r) offset for a new triad
# such that the new triad's base occupies that neighbor position.
# The key is to map a neighbor direction to the (dq, dr) for the new triad's base.
TRIAD_ATTACHMENT_RULES = {
    (1, 0): (1, 0),   # Attach to right neighbor
    (1, -1): (1, -1), # Attach to top-right neighbor
    (0, -1): (0, -1), # Attach to top-left neighbor
    (-1, 0): (-1, 0), # Attach to left neighbor
    (-1, 1): (-1, 1), # Attach to bottom-left neighbor
    (0, 1): (0, 1),   # Attach to bottom-right neighbor
}

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
    Creates triads ensuring they form a single connected cluster with no gaps.
    """
    if count <= 0:
        return []

    triads = []
    placed_hex_coords = set()  # Set of all individual hex coordinates placed so far
    perimeter_edges = {}  # Map of (hex_coord, neighbor_direction) -> list of triad indices that could attach there

    # Place the first triad at origin
    color = colors[len(triads) % len(colors)]
    first_triad = PeptideTriad(0, 0, color)
    triads.append(first_triad)
    first_hexes = first_triad.hex_cells()
    placed_hex_coords.update(first_hexes)
    
    # Initialize perimeter: for each hex in the first triad, add its neighbors to potential attachment points
    for hq, hr in first_hexes:
        for n_dir in HEX_NEIGHBORS:
            neighbor_coord = (hq + n_dir[0], hr + n_dir[1])
            # Only track neighbors that are not part of the initial triad
            if neighbor_coord not in placed_hex_coords:
                key = (neighbor_coord, (-n_dir[0], -n_dir[1])) # Store the reverse direction
                if key not in perimeter_edges:
                    perimeter_edges[key] = []
                perimeter_edges[key].append(0) # Index of the first triad


    # Place subsequent triads
    for i in range(1, count):
        if not perimeter_edges:
            print(f"Warning: No perimeter edges available to place triad {i+1}. Stopping.")
            break

        # Get a list of all valid attachment points
        valid_attachment_points = list(perimeter_edges.keys())
        
        if not valid_attachment_points:
            print(f"Warning: No valid attachment points found for triad {i+1}. Stopping.")
            break
            
        # Choose a random valid attachment point
        attach_neighbor_coord, attach_direction_from_neighbor = random.choice(valid_attachment_points)
        
        # Calculate the base (q,r) for the new triad based on the attachment rule
        # The base of the new triad should be at attach_neighbor_coord
        base_q, base_r = attach_neighbor_coord[0], attach_neighbor_coord[1]
        
        # Ensure the rule exists for this direction
        if attach_direction_from_neighbor not in TRIAD_ATTACHMENT_RULES:
            print(f"Warning: No attachment rule for direction {attach_direction_from_neighbor}. Skipping.")
            # Remove this edge from perimeter as it's problematic
            del perimeter_edges[(attach_neighbor_coord, attach_direction_from_neighbor)]
            continue
            
        # The TRIAD_ATTACHMENT_RULES maps the direction from the neighbor to where the new triad's base should be.
        # Since we already have the target base (attach_neighbor_coord), we don't need to adjust base_q, base_r further.
        # But we need to verify that the full triad will fit.
        
        # Create a trial triad at this base position
        trial_triad = PeptideTriad(base_q, base_r, colors[i % len(colors)])
        trial_cells = trial_triad.hex_cells()

        # Check if the entire trial triad fits without overlapping existing hexes
        if all(cell not in placed_hex_coords for cell in trial_cells):
            # Success! Place the triad
            triads.append(trial_triad)
            placed_hex_coords.update(trial_cells)
            
            # Update perimeter_edges:
            # 1. Remove the edge where we just attached
            del perimeter_edges[(attach_neighbor_coord, attach_direction_from_neighbor)]
            
            # 2. Add new perimeter edges created by this triad
            for hq, hr in trial_cells:
                for n_dir in HEX_NEIGHBORS:
                    neighbor_coord = (hq + n_dir[0], hr + n_dir[1])
                    if neighbor_coord not in placed_hex_coords:
                        key = (neighbor_coord, (-n_dir[0], -n_dir[1]))
                        if key not in perimeter_edges:
                            perimeter_edges[key] = []
                        perimeter_edges[key].append(i) # Index of the new triad
                        
            # 3. Remove any edges from perimeter that are now occupied by the new triad
            for cell in trial_cells:
                for n_dir in HEX_NEIGHBORS:
                    # Check if this cell+direction was a perimeter edge pointing to another triad
                    reverse_key = (cell, n_dir)
                    if reverse_key in perimeter_edges:
                         del perimeter_edges[reverse_key]
                         
        else:
            # This specific attachment point was invalid, remove it from options
            print(f"Warning: Trial placement failed for triad {i+1}, trying another point.")
            del perimeter_edges[(attach_neighbor_coord, attach_direction_from_neighbor)]
            
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
