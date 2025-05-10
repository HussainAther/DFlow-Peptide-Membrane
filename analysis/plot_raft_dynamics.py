import json
import matplotlib.pyplot as plt
from pathlib import Path

def plot_raft_dynamics(data):
    for i, run in enumerate(data):
        plt.plot(run["L_history"], label=f"L #{i}", alpha=0.2, color="blue")
        plt.plot(run["D_history"], label=f"D #{i}", alpha=0.2, color="red")
    
    plt.xlabel("Cycle")
    plt.ylabel("Membrane Thickness")
    plt.title("Raft Thickness Over Time (L vs D)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("plots/raft_dynamics.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    result_file = "experiments/lbias_60/dual_raft_results_*.json"  # Replace with actual file
    result_path = sorted(Path("experiments/lbias_60").glob("dual_raft_results_*.json"))[-1]
    
    with open(result_path) as f:
        results = json.load(f)
    
    plot_raft_dynamics(results)

