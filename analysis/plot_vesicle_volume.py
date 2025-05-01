import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

# ---------- ARGUMENT PARSING ----------
parser = argparse.ArgumentParser(description="Plot vesicle volume vs cycles")
parser.add_argument("logfile", type=str, help="Path to membrane_log.csv")
parser.add_argument("--outdir", type=str, default="plots", help="Output directory")
args = parser.parse_args()

logfile = Path(args.logfile)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(logfile)

# ---------- PLOTTING ----------
plt.figure(figsize=(10, 6))
plt.plot(df["Cycle"], df["VesicleVolume_nm3"], color="darkgreen", linewidth=2)
plt.xlabel("Cycle")
plt.ylabel("Vesicle Volume (nm³)")
plt.title("📈 Vesicle Volume Over Simulation")
plt.grid(True)

filename = f"vesicle_volume_{logfile.parent.name}.png"
plt.tight_layout()
plt.savefig(outdir / filename)
print(f"✅ Saved: {filename}")

