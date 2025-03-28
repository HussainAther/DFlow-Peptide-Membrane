import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Setup
log_dir = Path("logs_v2")
output_dir = log_dir / "plots"
output_dir.mkdir(parents=True, exist_ok=True)

# Load logs
growth_df = pd.read_csv(log_dir / "membrane_growth_log.csv")
peptides_df = pd.read_csv(log_dir / "peptide_log.csv")
vesicle_df = pd.read_csv(log_dir / "vesicle_state_log.csv")

# Style
sns.set(style="whitegrid")

# ---------- Plot 1: Membrane Thickness Over Time ----------
plt.figure(figsize=(10, 4))
plt.plot(growth_df["step"], growth_df["thickness_nm"], color="black")
plt.title("Membrane Thickness Over Time")
plt.xlabel("Peptide Index")
plt.ylabel("Thickness (nm)")
plt.tight_layout()
plt.savefig(output_dir / "membrane_thickness.png")
plt.close()

# ---------- Plot 2: Hydrophobicity Histogram by Fate ----------
plt.figure(figsize=(8, 5))
sns.histplot(data=peptides_df, x="hydrophobicity", hue="location", kde=True, bins=30)
plt.title("Peptide Hydrophobicity Distribution by Fate")
plt.xlabel("Avg Hydrophobicity (Kyte-Doolittle)")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig(output_dir / "hydrophobicity_by_fate.png")
plt.close()

# ---------- Plot 3: Insertion Type Breakdown ----------
counts = peptides_df["location"].value_counts()
plt.figure(figsize=(6, 4))
sns.barplot(x=counts.index, y=counts.values, palette="pastel")
plt.title("Peptide Fate Breakdown")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig(output_dir / "fate_breakdown.png")
plt.close()

# ---------- Plot 4: Chirality Breakdown ----------
chirality_counts = peptides_df["dominant_chirality"].value_counts(dropna=True)
plt.figure(figsize=(4, 4))
plt.pie(chirality_counts, labels=chirality_counts.index, autopct="%1.1f%%", startangle=150)
plt.title("Dominant Chirality of Inserted Peptides")
plt.savefig(output_dir / "chirality_pie.png")
plt.close()

# ---------- Plot 5: Raft L vs D Over Time ----------
plt.figure(figsize=(10, 4))
plt.plot(vesicle_df["step"], vesicle_df["raft_L"], label="Raft L", color="blue")
plt.plot(vesicle_df["step"], vesicle_df["raft_D"], label="Raft D", color="red")
plt.title("Growth of L vs D Peptide Rafts")
plt.xlabel("Step")
plt.ylabel("Raft Peptides Count")
plt.legend()
plt.tight_layout()
plt.savefig(output_dir / "raft_growth.png")
plt.close()

print(f"[✔] All v2 plots saved to: {output_dir}")

