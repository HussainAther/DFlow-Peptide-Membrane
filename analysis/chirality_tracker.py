import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Load the peptide log
log_path = Path("experiments/exp_latest/logs/peptide_log.csv")  # Adjust if needed
df = pd.read_csv(log_path)

# Count residues per chirality
def count_chirality(seq):
    return {
        "L": seq.count("L"),
        "D": seq.count("D"),
        "0": seq.count("0")
    }

chirality_data = df["sequence"].apply(count_chirality)
chirality_df = pd.DataFrame(list(chirality_data))
df = pd.concat([df, chirality_df], axis=1)

# Normalize per peptide
df["total"] = df[["L", "D", "0"]].sum(axis=1)
df["L_frac"] = df["L"] / df["total"]
df["D_frac"] = df["D"] / df["total"]
df["0_frac"] = df["0"] / df["total"]

# Add cumulative step index
df["step"] = df.index

# Moving average to smooth trends
window = 100
df["L_avg"] = df["L_frac"].rolling(window=window).mean()
df["D_avg"] = df["D_frac"].rolling(window=window).mean()
df["0_avg"] = df["0_frac"].rolling(window=window).mean()

# Plot
plt.figure(figsize=(10, 6))
plt.plot(df["step"], df["L_avg"], label="L-frac", color="blue")
plt.plot(df["step"], df["D_avg"], label="D-frac", color="red")
plt.plot(df["step"], df["0_avg"], label="0-frac", color="gray")
plt.xlabel("Peptide Index (Step)")
plt.ylabel(f"Fraction per Peptide (rolling avg {window})")
plt.title("Chirality Drift Over Time")
plt.legend()
plt.tight_layout()

output_dir = Path("analysis/plots")
output_dir.mkdir(parents=True, exist_ok=True)
plt.savefig(output_dir / "chirality_drift.png", dpi=300)
plt.close()

print("[✔] Chirality tracker plot saved to analysis/plots/chirality_drift.png")

