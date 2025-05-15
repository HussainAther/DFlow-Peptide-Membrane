import numpy as np
import matplotlib.pyplot as plt
from analysis.layout_min_energy import generate_initial_grid, optimize_layout, calculate_energy

def spiral_layout(anchor_lengths, grid_size):
    """Arrange peptides in spiral from longest (center) to shortest (edges)."""
    anchor_lengths = sorted(anchor_lengths, reverse=True)
    grid = np.zeros((grid_size, grid_size), dtype=int)
    cx, cy = grid_size // 2, grid_size // 2  # start in center

    directions = [(0,1), (1,0), (0,-1), (-1,0)]  # right, down, left, up
    steps = 1
    x, y = cx, cy
    idx = 0

    while idx < len(anchor_lengths):
        for d in directions:
            for _ in range(steps):
                if 0 <= x < grid_size and 0 <= y < grid_size and idx < len(anchor_lengths):
                    grid[x, y] = anchor_lengths[idx]
                    idx += 1
                x += d[0]
                y += d[1]
            if d in [(0,-1), (0,1)]:
                steps += 1  # increase steps every other turn
    return grid

def plot_comparison(grid1, grid2, energy1, energy2):
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    
    for ax, grid, title, energy in zip(
        axs,
        [grid1, grid2],
        ["Spiral Layout", "Min Energy Layout"],
        [energy1, energy2]
    ):
        cax = ax.matshow(grid, cmap='viridis')
        for (i, j), val in np.ndenumerate(grid):
            ax.text(j, i, str(val), ha='center', va='center',
                    color='white' if val > 5 else 'black')
        ax.set_title(f"{title}\n(Energy = {energy:.1f})")
        ax.axis('off')
    fig.colorbar(cax, ax=axs, orientation='vertical')
    plt.tight_layout()
    plt.show()

# ---- Example Usage ----
if __name__ == "__main__":
    anchor_lengths = [7, 6, 6, 5, 5, 4, 4, 3, 3]
    grid_size = 4

    spiral = spiral_layout(anchor_lengths, grid_size)
    spiral_energy = calculate_energy(spiral)

    initial = generate_initial_grid(anchor_lengths, grid_size)
    min_energy, min_energy_val = optimize_layout(initial)

    plot_comparison(spiral, min_energy, spiral_energy, min_energy_val)

