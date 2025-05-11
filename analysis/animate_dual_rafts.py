import matplotlib.pyplot as plt
import matplotlib.animation as animation
from src.hex_raft import Peptide, HexRaft
import numpy as np
import random
from pathlib import Path

def axial_to_cartesian(q, r):
    x = 3 / 2 * q
    y = np.sqrt(3) * (r + q / 2)
    return x, y

def generate_peptides(chirality, total=40, min_anchor=5, max_anchor=20):
    return [Peptide(chirality, random.randint(min_anchor, max_anchor)) for _ in range(total)]

def build_spiral_frames(peptides):
    raft = HexRaft(peptides[0].chirality)
    frames = []

    for i in range(1, len(peptides) + 1):
        raft.peptides = peptides[:i]
        layout = raft.sort_and_place()
        frames.append(dict(layout=layout.copy(), perimeter_avg=raft.get_perimeter_avg(), count=i))

    return frames

def animate_dual_rafts(l_frames, d_frames, save_path="plots/dual_raft_growth.gif"):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    def update(i):
        ax1.clear()
        ax2.clear()

        def plot_raft(ax, frame, chirality, color):
            layout = frame["layout"]
            for (q, r), peptide in layout.items():
                x, y = axial_to_cartesian(q, r)
                ax.add_patch(plt.Circle((x, y), 0.45, color=color, alpha=0.3))
                ax.text(x, y, str(peptide.anchor_len), ha="center", va="center", fontsize=7)

            ax.set_title(f"{chirality}-Raft\nPeptides: {frame['count']} | Avg: {frame['perimeter_avg']:.2f}")
            ax.axis("off")
            ax.set_aspect("equal")

        idx = min(i, len(l_frames) - 1, len(d_frames) - 1)
        plot_raft(ax1, l_frames[idx], "L", "blue")
        plot_raft(ax2, d_frames[idx], "D", "red")

    ani = animation.FuncAnimation(fig, update, frames=len(l_frames), interval=400, repeat=False)
    Path("plots").mkdir(exist_ok=True)
    ani.save(save_path, writer='pillow', fps=2)
    print(f"Saved animation to {save_path}")

if __name__ == "__main__":
    l_peptides = generate_peptides("L", total=40)
    d_peptides = generate_peptides("D", total=40)

    l_frames = build_spiral_frames(l_peptides)
    d_frames = build_spiral_frames(d_peptides)

    animate_dual_rafts(l_frames, d_frames, save_path="plots/dual_raft_growth.gif")

