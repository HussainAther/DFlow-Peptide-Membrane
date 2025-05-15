import numpy as np
import matplotlib.pyplot as plt

def generate_initial_grid(anchor_lengths, grid_size):
    """Place anchor lengths randomly in a 2D grid."""
    padded = anchor_lengths + [0] * (grid_size**2 - len(anchor_lengths))
    np.random.shuffle(padded)
    return np.array(padded).reshape((grid_size, grid_size))

def calculate_energy(grid):
    """Calculate energy as sum of differences between neighboring peptides."""
    energy = 0
    rows, cols = grid.shape
    for i in range(rows):
        for j in range(cols):
            center = grid[i, j]
            for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:  # 4-neighbor connectivity
                ni, nj = i + dx, j + dy
                if 0 <= ni < rows and 0 <= nj < cols:
                    energy += abs(center - grid[ni, nj])
    return energy / 2  # avoid double counting

def swap_and_calc_energy(grid, i1, j1, i2, j2):
    """Try swapping two grid values and calculate new energy."""
    new_grid = grid.copy()
    new_grid[i1, j1], new_grid[i2, j2] = new_grid[i2, j2], new_grid[i1, j1]
    return new_grid, calculate_energy(new_grid)

def optimize_layout(grid, iterations=1000):
    """Greedy optimizer that swaps elements to minimize total energy."""
    current_grid = grid.copy()
    current_energy = calculate_energy(current_grid)
    rows, cols = current_grid.shape

    for _ in range(iterations):
        i1, j1 = np.random.randint(rows), np.random.randint(cols)
        i2, j2 = np.random.randint(rows), np.random.randint(cols)
        new_grid, new_energy = swap_and_calc_energy(current_grid, i1, j1, i2, j2)
        if new_energy < current_energy:
            current_grid, current_energy = new_grid, new_energy

    return current_grid, current_energy

def visualize_grid(grid, title="Optimized Peptide Raft Layout"):
    """Plot grid with anchor length labels."""
    fig, ax = plt.subplots(figsize=(6, 6))
    cax = ax.matshow(grid, cmap='viridis')
    for (i, j), val in np.ndenumerate(grid):
        ax.text(j, i, str(int(val)), ha='center', va='center',
                color='white' if val > 5 else 'black')
    plt.title(title)
    plt.colorbar(cax, label="Anchor Length")
    plt.axis('off')
    plt.tight_layout()
    plt.show()

# ---- Example Usage ----
if __name__ == "__main__":
    anchor_lengths = [7, 6, 6, 5, 5, 4, 4, 3, 3]
    grid_size = 4
    init_grid = generate_initial_grid(anchor_lengths, grid_size)
    opt_grid, opt_energy = optimize_layout(init_grid)
    visualize_grid(opt_grid, f"Min-Energy Raft (Energy={opt_energy:.1f})")

