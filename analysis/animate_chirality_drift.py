import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from pathlib import Path
import argparse

# -------------- ARGPARSE ----------------
parser = argparse.ArgumentParser(description="Animate chirality drift over time")
parser.add_argument('logfile', type=str, help="Path to membrane_log.csv")
parser.add_argument('--outdir', type=str, default="plots", help="Where to save the animation")
parser.add_argument('--format', type=str, default="gif", choices=["gif", "mp4"], help="Output format")
args = parser.parse_args()

logfile = Path(args.logfile)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(logfile)

# -------------- SETUP FIGURE ----------------
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(0, df['Cycle'].max())
ax.set_ylim(0, max(df['L_to_D_Ratio'].max(), 2))
ax.set_xlabel("Cycle")
ax.set_ylabel("L:D Raft Ratio")
ax.set_title("Chirality Drift Over Time")

line, = ax.plot([], [], color='purple', linewidth=2)
text = ax.text(0.95, 0.95, '', ha='right', va='top', transform=ax.transAxes)

# -------------- ANIMATION FUNCTION ----------------
def update(frame):
    x = df['Cycle'][:frame]
    y = df['L_to_D_Ratio'][:frame]
    line.set_data(x, y)
    text.set_text(f"Cycle: {x.iloc[-1] if not x.empty else 0}")
    return line, text

ani = FuncAnimation(fig, update, frames=len(df), interval=20, blit=True)

# -------------- SAVE ----------------
basename = logfile.parent.name
out_path = outdir / f"chirality_drift_{basename}.{args.format}"

if args.format == "gif":
    ani.save(out_path, writer=PillowWriter(fps=30))
else:
    ani.save(out_path, writer=FFMpegWriter(fps=30))

print(f"✅ Saved animation: {out_path}")

