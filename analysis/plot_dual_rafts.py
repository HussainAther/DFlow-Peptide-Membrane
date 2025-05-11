import matplotlib.pyplot as plt
from src.hex_raft import Peptide, HexRaft
import random
import numpy as np
from pathlib import Path

def axial_to_cartesian(q, r):
    x = 3 / 2 * q
    y = np.sqrt(3) * (r + q / 2)
    return x, y

def plot_raft_subplot(ax, raft: HexRaft, title=None, color="blue"):
    layout = raft.sort_and_place()
    avg = raft.get_perimeter_avg()

    for (q, r), peptide in layout.items():
        x, y = axial_to_cartesian(q, r)
        ax.add_patch(plt.Circle((x, y), 0.45, color=color, alpha=0.25))
        ax.text(x, y, str(peptide.anchor_len), ha="center", va="center", fontsize=8)

    ax.set_title(f"{title} (Perimeter Avg: {avg:.2f})", fontsize=10)
    ax.set_aspect("equal")
    ax.axis("off")

def generate_random_peptides(chirality, count=40):
    raft = HexRaft(chirality)
    for _ in range(count):
        anchor = random.randint(5, 20)
        raft.add_peptide(Peptide(chirality, anchor))
    return raft

if __name__ == "__main__":
    l_raft = generate_random_peptides("L", 40)
    d_raft = generate_random_peptides("D", 40)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    plot_raft_subplot(ax1, l_raft, title="L-Raft", color="blue")
    plot_raft_subplot(ax2, d_raft, title="D-Raft", color="red")

    plt.suptitle("Peptide Spiral Layouts in L and D Rafts", fontsize=14)
    plt.tight_layout()
    Path("plots").mkdir(exist_ok=True)
    plt.savefig("plots/dual_rafts_spiral.png", dpi=300)
    plt.show()

