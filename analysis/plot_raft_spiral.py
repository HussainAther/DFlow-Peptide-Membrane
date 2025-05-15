import matplotlib.pyplot as plt
import numpy as np


def spiral_coordinates(n):
    """Generate (x, y) coordinates in a spiral layout for n points."""
    coords = [(0, 0)]
    angle = 0
    radius = 0.5
    step = 0.5

    for i in range(1, n):
        angle += np.pi / 3  # 60° per step gives a hexagonal feel
        radius += step / (2 * np.pi)
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        coords.append((x, y))

    return coords


def plot_peptide_raft(peptides, anchor_lengths, chirality):
    """
    Plot a 2D peptide raft using a spiral layout.
    
    Parameters:
    - peptides: list of peptide IDs
    - anchor_lengths: list of anchor lengths (int)
    - chirality: list of 'L' or 'D'
    """
    # Sort by anchor length descending
    sorted_data = sorted(zip(peptides, anchor_lengths, chirality), key=lambda x: -x[1])
    coords = spiral_coordinates(len(sorted_data))
    xs, ys = zip(*coords)

    fig, ax = plt.subplots(figsize=(6, 6))
    for (pid, anchor, ch), x, y in zip(sorted_data, xs, ys):
        color = "blue" if ch == 'L' else "red"
        ax.scatter(x, y, s=100, color=color, alpha=0.7, edgecolor='black')
        ax.text(x, y + 0.1, f"{anchor}", ha='center', fontsize=8)
        ax.text(x, y - 0.15, f"{pid}", ha='center', fontsize=7, color='gray')

    ax.set_aspect('equal')
    ax.set_title("Peptide Raft Layout (Anchor Length Spiral)")
    ax.axis('off')
    plt.tight_layout()
    plt.show()


# Example usage
if __name__ == "__main__":
    peptides = [f"P{i}" for i in range(10)]
    anchor_lengths = [7, 6, 6, 5, 5, 4, 4, 3, 3, 2]
    chirality = ['L', 'L', 'D', 'L', 'D', 'D', 'L', 'D', 'L', 'D']

    plot_peptide_raft(peptides, anchor_lengths, chirality)

