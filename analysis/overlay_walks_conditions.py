import matplotlib.pyplot as plt
from src.coupled_thickness_walk import run_coupled_walks


def simulate_condition_series(param_name, param_values, **kwargs):
    results = []
    for value in param_values:
        print(f"Running {param_name} = {value}")
        custom_args = kwargs.copy()
        custom_args[param_name] = value
        hL, hD = run_coupled_walks(**custom_args)
        results.append((value, hL, hD))
    return results


def plot_overlay(results, label_prefix=""):
    plt.figure(figsize=(10, 6))
    for value, hL, hD in results:
        plt.plot(hL, label=f"{label_prefix}{value} (L)", alpha=0.7)
        plt.plot(hD, label=f"{label_prefix}{value} (D)", alpha=0.7, linestyle='--')
    plt.xlabel("Cycle")
    plt.ylabel("Thickness")
    plt.title("Overlay of L/D Thickness Walks Under Varying Conditions")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    glycine_values = [0.01, 0.05, 0.10]  # Change in achiral content
    results = simulate_condition_series(
        param_name="gly_prob",  # This assumes `gly_prob` is handled in your model (next step)
        param_values=glycine_values,
        k_spring=0.5,
        variable_jump=True,
        start_L=12,
        start_D=12,
        min_thick=10,
        max_thick=23
    )
    plot_overlay(results, label_prefix="gly=")

