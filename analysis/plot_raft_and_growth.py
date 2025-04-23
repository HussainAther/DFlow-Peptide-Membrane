import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# === Config ===
log_dir = Path("experiments/exp_latest/logs")
peptide_path = log_dir / "peptide_log.csv"
vesicle_path = log_dir / "vesicle_state_log.csv"
out_dir = Path("analysis/plots")
out_dir.mkdir(parents=True, exist_ok=True)

# === Load logs ===
peptides = pd.read_csv(peptide_path)
vesicle = pd.read_csv(vesicle_path)

# === Infer Raft Counts (L vs D insertions) ===
def dominant_chirality(seq):
    l = seq.count("L")
    d = seq.count("D")
    return "L" if l > d else "D" if d > l else "0"

inserted = peptides[peptides["inserted"] == True].copy()
inserted["raft"] = inserted["sequence"].apply(dominant_chirality)

# Cumulative L/D raft count
inserted["L_cumsum"] = (inserted["raft"] == "L").cumsum()
inserted["D_cumsum"] = (inserted["raft"] == "D").cumsum()

# === Plot 1: Raft Asymmetry ===
plt.figure(figsize=(8, 5))
plt.plot(inserted["id"], inserted["L_cumsum"], label="L-Rafts", color="blue")
plt.plot(inserted["id"], inserted["D_cumsum"], label="D-Rafts", color="red")
plt.title("🏗️ Raft Asymmetry Over Time")
plt.xlabel("Peptide Insertions")
plt.ylabel("Cumulative Count")
plt.legend()
plt.tight_layout()
plt.savefig(out_dir / "raft_asymmetry.png", dpi=300)
plt.close()

# === Plot 2: Membrane Growth Phases ===
plt.figure(figsize=(8, 5))
plt.plot(vesicle["step"], vesicle["membrane_thickness"], label="Thickness", color="green")
plt.title("🧱 Membrane Growth Over Time")
plt.xlabel("Simulation Step")
plt.ylabel("Membrane Thickness (nm)")
plt.tight_layout()
plt.savefig(out_dir / "membrane_growth.png", dpi=300)
plt.close()

print("[✔] Plots saved to:", out_dir)

