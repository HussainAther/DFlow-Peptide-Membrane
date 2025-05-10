import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def plot_drift(results):
    max_len = max(len(run["L_history"]) for run in results)
    drift_matrix = []

    for run in results:
        l_hist = run["L_history"]
        d_hist = run["D_history"]
        drift = [l - d for l, d in zip(l_hist, d_hist)]
        padded = drift + [0] * (max_len - len(drift))
        drift_matrix.append(padded)

    drift_matrix = np.array(drift_matrix)

    plt.imshow(drift_matrix, aspect='auto', cmap='coolwarm', interpolation='nearest')
    plt.colorbar(label="L - D Thickness (nm)")
    plt.xlabel("Cycle")
    plt.ylabel("Run #")
    plt.title("Heatmap of Raft Thickness Drift (L - D)")
    plt.tight_layout()
    plt.savefig("plots/drift_heatmap.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    result_path = sorted(Path("experiments/lbias_60").glob("dual_raft_results_*.json"))[-1]

    with open(result_path) as f:
        results = json.load(f)

    plot_drift(results)

