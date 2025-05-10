import json
import matplotlib.pyplot as plt
from pathlib import Path

def plot_crowded(results):
    L = [r["L_crowded"] for r in results]
    D = [r["D_crowded"] for r in results]

    plt.hist(L, bins=20, alpha=0.6, label="L-Crowded", color="blue")
    plt.hist(D, bins=20, alpha=0.6, label="D-Crowded", color="red")
    plt.xlabel("Crowded Peptides")
    plt.ylabel("Run Count")
    plt.title("Crowded (Uninserted) Peptides per Run")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("plots/crowded_histogram.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    result_path = sorted(Path("experiments/lbias_60").glob("dual_raft_results_*.json"))[-1]
    
    with open(result_path) as f:
        results = json.load(f)

    plot_crowded(results)

