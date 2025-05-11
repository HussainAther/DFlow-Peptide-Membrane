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

    for i in range(1, len(peptides)+1):
        raft.peptides = peptides[:i]
        layout = raft.sort_and_place()
        frames.append(dict(layout=layout.copy(), perimeter_avg=raft.get_perimeter_avg()))

    return frames

def animate_raft(frames, chirality, save_path="plots/raft_growth.gif"):
    fig, ax = plt.subplots(figsize=(6, 6))

    def update(frame_data):
        ax.clear()
        layout = frame_data["layout"]
        avg = frame_data["perimeter_avg"]

        for (q, r), peptide in layout.items():
            x, y = axial_to_cartesian(q, r)
            ax.add_patch(plt.Circle((x, y), 0.45, color="blue" if chirality == "L" else "red", alpha=0.3))
            ax.text(x, y, str(peptide.anchor_len), ha="center", va="center", fontsize=7)

        ax.set_title(f"{chirality}-Raft Growth\nPeptides: {len(layout)} | Perimeter Avg: {avg:.2f}")
        ax.axis("off")
        ax.set_aspect("equal")

    ani = animation.FuncAnimation(fig, update, frames=frames, interval=400, repeat=False)
    Path("plots").mkdir(exist_ok=True)
    ani.save(save_path, writer='pillow', fps=2)
    print(f"Animation saved to {save_path}")

if __name__ == "__main__":
    peptides = generate_peptides("L", total=40)
    frames = build_spiral_frames(peptides)
    animate_raft(frames, chirality="L", save_path="plots/l_raft_growth.gif")

