import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

# --------- ARGUMENTS ---------
parser = argparse.ArgumentParser(description="Plot peptide crowding density over time")
parser.add_argument("logfile", type=str, help="Path to membrane_log.csv")
parser.add_argument("--outdir", type=str, default="plots", help="Where to save plot")
args = parser.parse_args()

logfile = Path(args.logfile)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# --------- LOAD LOG ---------
df = pd.read_csv(logfile)

# --------- PLOT ---------
plt.figure(figsize=(10, 6))
plt.plot(df["Cycle"], df["PeptideDensity_per_nm2"], color="darkorange", linewidth=2)
plt.xlabel("Cycle")
plt.ylabel("Peptide Density (residues / nm²)")
plt.title("📊 Peptide Crowding Density on Vesicle Membrane")
plt.grid(True)

# --------- SAVE ---------
filename = f"density_drift_{logfile.parent.name}.png"
plt.tight_layout()
plt.savefig(outdir / filename)
print(f"✅ Saved peptide density plot: {outdir / filename}")

