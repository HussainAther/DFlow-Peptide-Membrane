import matplotlib.pyplot as plt
from src.coupled_thickness_walk import run_coupled_walks
import numpy as np

def simulate_absorption_times(n_runs=200, **kwargs):
    cycle_counts = []

    for _ in range(n_runs):
        hL, hD = run_coupled_walks(**kwargs)
        cycles = min(len(hL), len(hD))  # Stop when either absorbs
        cycle_counts.append(cycles)

    return cycle_counts


def plot_histogram(cycle_counts):
    plt.hist(cycle_counts, bins=20, edgecolor='black', alpha=0.8)
    plt.xlabel("Cycles before absorption")
    plt.ylabel("Frequency")
    plt.title("Distribution of Vesicle Lifespans")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    print("🧪 Running absorption time simulations...")
    cycles_list = simulate_absorption_times(
        n_runs=250,
        k_spring=0.4,
        variable_jump=True,
        start_L=12,
        start_D=12,
        min_thick=10,
        max_thick=23
    )
    print(f"✅ Completed {len(cycles_list)} simulations.")
    plot_histogram(cycles_list)

