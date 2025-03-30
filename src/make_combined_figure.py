import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

log_dir = Path("logs_v2")
plot_dir = log_dir / "plots"
plot_dir.mkdir(exist_ok=True)

growth_df = pd.read_csv(log_dir / "membrane_growth_log.csv")
peptides_df = pd.read_csv(log_dir / "peptide_log.csv")
vesicle_df = pd.read_csv(log_dir / "vesicle_state_log.csv")

sns.set(style="whitegrid")
fig, axs = plt.subplots(2, 2, figsize=(12, 8))

# Plot 1: Membrane Thickness
axs[0, 0].plot(growth_df["step"], growth_df["thickness_nm"], color="black")
axs[0, 0].set_title("Membrane Thickness Over Time")
axs[0, 0].set_xlabel("Step")
axs[0, 0].set_ylabel("Thickness (nm)")

# Plot 2: Hydrophobicity Histogram
sns.histplot(data=peptides_df, x="hydrophobicity", hue="location",
             ax=axs[0, 1], bins=30, kde=True)
axs[0, 1].set_title("Hydrophobicity by Peptide Fate")

# Plot 3: Raft Growth
axs[1, 0].plot(vesicle_df["step"], vesicle_df["raft_L"], label="Raft L", color="blue")
axs[1, 0].plot(vesicle_df["step"], vesicle_df["raft_D"], label="Raft D", color="red")
axs[1, 0].set_title("Raft L vs D Population")
axs[1, 0].set_xlabel("Step")
axs[1, 0].set_ylabel("Peptides")
axs[1, 0].legend()

# Plot 4: Fate Breakdown
fate_counts = peptides_df["location"].value_counts()
sns.barplot(x=fate_counts.index, y=fate_counts.values, ax=axs[1, 1], palette="pastel")
axs[1, 1].set_title("Peptide Fate Distribution")
axs[1, 1].set_ylabel("Count")

plt.tight_layout()
fig.savefig(plot_dir / "combined_figure.png", dpi=300)
fig.savefig(plot_dir / "combined_figure.svg")  # For LaTeX
plt.close()
print("[✔] Combined figure created.")

