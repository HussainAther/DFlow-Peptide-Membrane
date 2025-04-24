import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

batch_dir = Path("experiments/batch_runs")
run_folders = [f for f in batch_dir.iterdir() if f.is_dir()]

for run in run_folders:
    log_file = run / "membrane_log.csv"
    if not log_file.exists():
        continue

    df = pd.read_csv(log_file)

    fig, ax1 = plt.subplots(figsize=(10,6))

    ax1.plot(df['Cycle'], df['MembraneThickness'], label='Membrane Thickness', color='black')
    ax1.set_xlabel('Cycle')
    ax1.set_ylabel('Membrane Thickness', color='black')
    ax1.tick_params(axis='y', labelcolor='black')

    ax2 = ax1.twinx()
    ax2.plot(df['Cycle'], df['L_Raft_Size'], label='L-Raft', color='blue')
    ax2.plot(df['Cycle'], df['D_Raft_Size'], label='D-Raft', color='red')
    ax2.set_ylabel('Raft Sizes', color='gray')
    ax2.tick_params(axis='y', labelcolor='gray')

    plt.title(f"Simulation: {run.name}")
    fig.legend(loc="upper left", bbox_to_anchor=(0.1,0.9))
    plt.tight_layout()

    plt.savefig(run / "growth_raft_plot.png")
    plt.close()

print("✅ Plots generated for all batch runs!")

