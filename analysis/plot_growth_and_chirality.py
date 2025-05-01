import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

# -------------- ARGPARSE ----------------
parser = argparse.ArgumentParser(description="Plot membrane growth + chirality drift")
parser.add_argument('logfile', type=str, help="Path to membrane_log.csv")
parser.add_argument('--outdir', type=str, default="plots", help="Directory to save the figure")
args = parser.parse_args()

logfile = Path(args.logfile)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# -------------- LOAD DATA ----------------
df = pd.read_csv(logfile)

# -------------- SETUP PLOT ----------------
fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# Subplot 1: Membrane thickness growth
axs[0].plot(df['Cycle'], df['MembraneThickness'], color='teal', linewidth=2)
axs[0].set_ylabel("Membrane Thickness (a.a.)")
axs[0].set_title("Membrane Growth Over Time")
axs[0].grid(True)

# Subplot 2: Chirality Drift
axs[1].plot(df['Cycle'], df['L_to_D_Ratio'], color='purple', linewidth=2)
axs[1].axhline(1, color='gray', linestyle='--', linewidth=0.7)
axs[1].set_ylabel("L/D Raft Ratio")
axs[1].set_xlabel("Simulation Cycle")
axs[1].set_title("Chirality Drift Over Time")
axs[1].grid(True)

plt.tight_layout()

# -------------- SAVE ----------------
basename = logfile.parent.name
fig_path = outdir / f"growth_chirality_{basename}.png"
plt.savefig(fig_path)
print(f"✅ Saved multi-panel plot: {fig_path}")

