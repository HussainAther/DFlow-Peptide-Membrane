import random
import matplotlib.pyplot as plt


def run_thickness_walk(
    start=12,
    min_thick=10,
    max_thick=23,
    p_grow=0.5,
    total_cycles=1000,
    variable_jump=False,
    verbose=False
):
    thickness = start
    history = [thickness]

    for cycle in range(total_cycles):
        anchor_len = thickness + (random.choice([-1, 0, 1]) if variable_jump else 1)

        if anchor_len > thickness:
            thickness += 1
        elif anchor_len < thickness:
            thickness -= 1
        # else: no change

        # Enforce boundaries (absorbing states)
        if thickness <= min_thick:
            if verbose:
                print(f"Ruin! Hit min thickness at cycle {cycle}")
            break
        if thickness >= max_thick:
            if verbose:
                print(f"Max growth! Hit max thickness at cycle {cycle}")
            break

        history.append(thickness)

    return history


def plot_thickness_walk(history):
    plt.plot(history, marker='o', linewidth=2)
    plt.xlabel("Cycle")
    plt.ylabel("Membrane Thickness (residues)")
    plt.title("Membrane Thickness Walk (Gambler's Ruin Analogy)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    history = run_thickness_walk(
        start=12,
        min_thick=10,
        max_thick=23,
        variable_jump=True,
        verbose=True
    )
    plot_thickness_walk(history)

