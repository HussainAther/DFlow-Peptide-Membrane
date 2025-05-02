import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

# --------- ARGS ---------
parser = argparse.ArgumentParser(description="Compare flattened vs spherical vesicle runs")
parser.add_argument("spherical_log", type=str, help="Path to spherical membrane_log.csv")
parser.add_argument("flattened_log", type=str, help="Path to flattened membrane_log.csv")
parser.add_argument("--outdir", type=str, default="plots", help="Where to save output")
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# --------- LOAD ---------
df_sph = pd.read_csv(args.spherical_log)
df_flat = pd.read_csv(args.flattened_log)

# --------- FIGURE ---------
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle("🫧 Vesicle Shape Comparison: Spherical vs Flattened", fontsize=16)

# 1. Volume
axs[0, 0].plot(df_sph["Cycle"], df_sph["VesicleVolume_nm3"], label="Spherical", color="blue")
axs[0, 0].plot(df_flat["Cycle"], df_flat["VesicleVolume_nm3"], label="Flattened", color="orange")
axs[0, 0].set_title("Vesicle Volume")
axs[0, 0].set_ylabel("Volume (nm³)")
axs[0, 0].legend()
axs[0, 0].grid(True)

# 2. Peptide Density
axs[0, 1].plot(df_sph["Cycle"], df_sph["PeptideDensity_per_nm2"], label="Spherical", color="blue")
axs[0, 1].plot(df_flat["Cycle"], df_flat["PeptideDensity_per_nm2"], label="Flattened", color="orange")
axs[0, 1].set_title("Peptide Density")
axs[0, 1].set_ylabel("Density (residues / nm²)")
axs[0, 1].legend()
axs[0, 1].grid(True)

# 3. L:D Raft Ratio
axs[1, 0].plot(df_sph["Cycle"], df_sph["L_to_D_Ratio"], label="Spherical", color="blue")
axs[1, 0].plot(df_flat["Cycle"], df_flat["L_to_D_Ratio"], label="Flattened", color="orange")
axs[1, 0].set_title("L:D Raft Ratio")
axs[1, 0].set_ylabel("L:D")
axs[1, 0].legend()
axs[1, 0].grid(True)

# 4. Membrane Thickness
axs[1, 1].plot(df_sph["Cycle"], df_sph["MembraneThickness"], label="Spherical", color="blue")
axs[1, 1].plot(df_flat["Cycle"], df_flat["MembraneThickness"], label="Flattened", color="orange")
axs[1, 1].set_title("Membrane Thickness")
axs[1, 1].set_ylabel("Amino Acid Units")
axs[1, 1].legend()
axs[1, 1].grid(True)

for ax in axs.flat:
    ax.set_xlabel("Cycle")

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
fname = outdir / "vesicle_shape_comparison.png"
plt.savefig(fname)
print(f"✅ Comparison saved to: {fname}")

