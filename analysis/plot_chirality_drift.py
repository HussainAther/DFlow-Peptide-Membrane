import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

# -------------- ARGPARSE ----------------
parser = argparse.ArgumentParser(description="Plot Chirality Drift Over Time")
parser.add_argument('logfile', type=str, help="Path to membrane_log.csv")
parser.add_argument('--outdir', type=str, default="plots", help="Where to save the plot")
args = parser.parse_args()

logfile = Path(args.logfile)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# -------------- LOAD LOG ----------------
df = pd.read_csv(logfile)

# -------------- PLOT ----------------
plt.figure(figsize=(10, 6))
plt.plot(df['Cycle'], df['L_to_D_Ratio'], label='L/D Ratio', color='purple')
plt.axhline(y=1, color='gray', linestyle='--', linewidth=0.7)
plt.xlabel("Cycle")
plt.ylabel("L-to-D Insertion Ratio")
plt.title("🧬 Chirality Drift Over Time")
plt.legend()
plt.tight_layout()

# -------------- SAVE ----------------
plot_path = outdir / f"chirality_drift_{logfile.parent.name}.png"
plt.savefig(plot_path)
print(f"✅ Saved: {plot_path}")

