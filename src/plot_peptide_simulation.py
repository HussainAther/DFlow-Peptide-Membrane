import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Load logs
log_dir = Path("logs")
growth_df = pd.read_csv(log_dir / "membrane_growth_log.csv")
peptides_df = pd.read_csv(log_dir / "peptide_log.csv")

# Style
sns.set(style="whitegrid", palette="muted")

# ---------- Plot 1: Membrane Thickness Over Time ----------
plt.figure(figsize=(10, 4))
plt.plot(growth_df["step"], growth_df["thickness_nm"], label="Membrane Thickness")
plt.xlabel("Peptide Index")
plt.ylabel("Thickness (nm)")
plt.title("Membrane Thickness Over Time")
plt.tight_layout()
plt.savefig(log_dir / "membrane_thickness_plot.png")
plt.close()

# ---------- Plot 2: Hydrophobicity Histogram ----------
plt.figure(figsize=(6, 4))
sns.histplot(data=peptides_df, x="hydrophobicity", hue="location", kde=True, bins=30)
plt.title("Hydrophobicity Distribution by Peptide Fate")
plt.xlabel("Hydrophobicity Score (Kyte-Doolittle)")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig(log_dir / "hydrophobicity_histogram.png")
plt.close()

# ---------- Plot 3: Chirality Composition ----------
from collections import Counter

chirals = ''.join(''.join(seq.split('-')) for seq in peptides_df["sequence"])
counts = Counter(chirals)
total = sum(counts.values())

labels = ['L', 'D', '0']
sizes = [counts.get(k, 0)/total for k in labels]

plt.figure(figsize=(4, 4))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
plt.title("Chirality Composition (All Peptides)")
plt.savefig(log_dir / "chirality_pie.png")
plt.close()

# ---------- Plot 4: Hydrophobicity vs Insertion ----------
plt.figure(figsize=(6, 4))
sns.boxplot(x="location", y="hydrophobicity", data=peptides_df)
plt.title("Hydrophobicity vs Peptide Fate")
plt.ylabel("Avg Hydrophobicity")
plt.xlabel("Peptide Fate")
plt.tight_layout()
plt.savefig(log_dir / "hydrophobicity_vs_fate.png")
plt.close()

print("[✔] Plots generated in /logs/")

