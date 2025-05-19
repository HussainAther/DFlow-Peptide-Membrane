import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path

def spiral_layout(anchor_lengths, grid_size):
    grid = np.zeros((grid_size, grid_size), dtype=int)
    center = grid_size // 2
    x, y = center, center
    dx, dy = 0, -1
    i = 0

    for step in range(grid_size ** 2):
        if 0 <= x < grid_size and 0 <= y < grid_size:
            if i < len(anchor_lengths):
                grid[y][x] = anchor_lengths[i]
                i += 1
        if (x == y) or (x < center and x == -y) or (x > center and x == 1 - y):
            dx, dy = -dy, dx
        x += dx
        y += dy
    return grid

def generate_initial_grid(anchor_lengths, grid_size):
    shuffled = np.random.permutation(anchor_lengths)
    grid = np.zeros((grid_size, grid_size), dtype=int)
    indices = np.argwhere(grid == 0)
    for i, (y, x) in enumerate(indices):
        if i < len(shuffled):
            grid[y, x] = shuffled[i]
    return grid

def calculate_energy(grid):
    energy = 0
    rows, cols = grid.shape
    for i in range(rows):
        for j in range(cols):
            current = grid[i, j]
            if current == 0:
                continue
            for di, dj in [(-1,0),(1,0),(0,-1),(0,1)]:
                ni, nj = i + di, j + dj
                if 0 <= ni < rows and 0 <= nj < cols:
                    neighbor = grid[ni, nj]
                    if neighbor != 0:
                        energy += abs(current - neighbor)
    return energy / 2  # Each pair counted twice

def optimize_layout(grid, iterations=1000):
    best_grid = grid.copy()
    best_energy = calculate_energy(best_grid)

    for _ in range(iterations):
        a, b = np.random.randint(0, grid.shape[0], 2), np.random.randint(0, grid.shape[1], 2)
        c, d = np.random.randint(0, grid.shape[0], 2), np.random.randint(0, grid.shape[1], 2)
        grid[a[0], a[1]], grid[c[0], c[1]] = grid[c[0], c[1]], grid[a[0], a[1]]
        energy = calculate_energy(grid)
        if energy < best_energy:
            best_energy = energy
            best_grid = grid.copy()
        else:
            grid[a[0], a[1]], grid[c[0], c[1]] = grid[c[0], c[1]], grid[a[0], a[1]]
    return best_grid, best_energy

def find_perimeter_coords(grid):
    rows, cols = grid.shape
    perimeter = set()
    for i in range(rows):
        for j in range(cols):
            if grid[i, j] != 0:
                for di, dj in [(-1,0),(1,0),(0,-1),(0,1)]:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < rows and 0 <= nj < cols and grid[ni, nj] == 0:
                        perimeter.add((i, j))
    return list(perimeter)

def calculate_mismatch(grid, membrane_thickness):
    perimeter = find_perimeter_coords(grid)
    total_mismatch = sum(abs(grid[i, j] - membrane_thickness) for i, j in perimeter)
    return total_mismatch, perimeter

def plot_comparison(grid1, grid2, energy1, energy2, perimeter1, perimeter2,
                    save_path="plots/layout_comparison_mismatch.png",
                    membrane_thickness=5):
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    grids = [grid1, grid2]
    titles = [f"Spiral Layout\nEnergy={energy1:.1f}", f"Min-Energy Layout\nEnergy={energy2:.1f}"]
    perims = [perimeter1, perimeter2]

    for ax, grid, title, perim in zip(axs, grids, titles, perims):
        cax = ax.matshow(grid, cmap='viridis')
        for (i, j), val in np.ndenumerate(grid):
            ax.text(j, i, str(val), ha='center', va='center', color='white' if val > 5 else 'black')
        for (i, j) in perim:
            ax.add_patch(Rectangle((j-0.5, i-0.5), 1, 1, fill=False, edgecolor='red', linewidth=2))
        ax.set_title(title)
        ax.axis('off')

    fig.colorbar(cax, ax=axs.ravel().tolist(), shrink=0.8)
    plt.tight_layout()

    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path, dpi=300)
    plt.close()

    print(f"✅ Plot saved: {save_path.resolve()}")

# Run the comparison
if __name__ == "__main__":
    anchor_lengths = [7, 6, 6, 5, 5, 4, 4, 3, 3]
    grid_size = 4
    membrane_thickness = 5

    spiral = spiral_layout(anchor_lengths, grid_size)
    spiral_energy = calculate_energy(spiral)
    spiral_mismatch, perimeter1 = calculate_mismatch(spiral, membrane_thickness)

    initial = generate_initial_grid(anchor_lengths, grid_size)
    min_energy, min_energy_val = optimize_layout(initial)
    min_mismatch, perimeter2 = calculate_mismatch(min_energy, membrane_thickness)

    print(f"Spiral mismatch: {spiral_mismatch}")
    print(f"Min-Energy mismatch: {min_mismatch}")

    plot_comparison(spiral, min_energy, spiral_energy, min_energy_val, perimeter1, perimeter2)

