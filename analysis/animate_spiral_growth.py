import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from pathlib import Path

def spiral_order_coords(n):
    """Yield grid coordinates in spiral order for an n x n grid."""
    x, y = n // 2, n // 2
    dx, dy = 0, -1
    for _ in range(n ** 2):
        if 0 <= x < n and 0 <= y < n:
            yield y, x
        if (x == y) or (x < n//2 and x == -y) or (x > n//2 and x == 1 - y):
            dx, dy = -dy, dx
        x += dx
        y += dy

def generate_spiral_frames(anchor_lengths, grid_size=None, save_frames=False, save_dir="spiral_frames"):
    import math

    num_peptides = len(anchor_lengths)

    if grid_size is None:
        grid_size = math.ceil(math.sqrt(num_peptides))

    spiral_coords = list(spiral_order_coords(grid_size))
    if len(spiral_coords) < num_peptides:
        raise ValueError(f"Grid size {grid_size} too small for {num_peptides} peptides.")

    if save_frames:
        os.makedirs(save_dir, exist_ok=True)

    frames = []
    grid = np.zeros((grid_size, grid_size), dtype=int)

    for i in range(num_peptides):
        temp_grid = grid.copy()
        for j in range(i + 1):
            y, x = spiral_coords[j]
            temp_grid[y, x] = anchor_lengths[j]
        frames.append(temp_grid)

        if save_frames:
            fig, ax = plt.subplots(figsize=(4, 4))
            im = ax.imshow(temp_grid, cmap='viridis', vmin=0, vmax=max(anchor_lengths))
            for y in range(grid_size):
                for x in range(grid_size):
                    val = temp_grid[y, x]
                    if val > 0:
                        ax.text(x, y, f'{val}', ha='center', va='center', color='white')
            ax.axis('off')
            plt.tight_layout()
            plt.savefig(f"{save_dir}/frame_{i:03d}.png")
            plt.close()

    return frames

def animate_growth(frames, save_path="plots/spiral_growth.gif"):
    fig, ax = plt.subplots()
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=0, vmax=max(max(f.flatten()) for f in frames))

    def update(frame):
        ax.clear()
        im = ax.matshow(frame, cmap=cmap, norm=norm)
        for (i, j), val in np.ndenumerate(frame):
            if val > 0:
                ax.text(j, i, str(val), ha='center', va='center', color='white' if val > 5 else 'black')
        ax.set_title(f"Step {np.count_nonzero(frame)} / {len(frames)}")
        ax.axis('off')
        return [im]

    ani = animation.FuncAnimation(fig, update, frames=frames, interval=600)

    # Save
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    ani.save(save_path, dpi=150, writer='pillow')
    print(f"✅ Animation saved: {save_path.resolve()}")

if __name__ == "__main__":
    anchor_lengths = [7, 6, 6, 5, 5, 4, 4, 3, 3]
    grid_size = 5

    frames = generate_spiral_frames(anchor_lengths, grid_size)
    animate_growth(frames, save_path="plots/spiral_growth.gif")

