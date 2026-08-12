import random
import matplotlib.pyplot as plt


def spring_force(delta, k=1):
    """Simple harmonic spring force pulling toward zero delta."""
    return -k * delta


def run_coupled_walks(
    start_L=12,
    start_D=12,
    min_thick=10,
    max_thick=23,
    cycles=1000,
    k_spring=0.5,
    variable_jump=True,
    verbose=False
):
    thickness_L = start_L
    thickness_D = start_D
    hist_L = [thickness_L]
    hist_D = [thickness_D]

    for t in range(cycles):
        # Spring pull based on difference
        delta = thickness_L - thickness_D
        spring = spring_force(delta, k=k_spring)

        # L raft anchor
        anchor_L = thickness_L + random.choice([-1, 0, 1]) if variable_jump else thickness_L
        if anchor_L > thickness_L + spring:
            thickness_L += 1
        elif anchor_L < thickness_L + spring:
            thickness_L -= 1

        # D raft anchor
        anchor_D = thickness_D + random.choice([-1, 0, 1]) if variable_jump else thickness_D
        if anchor_D > thickness_D - spring:
            thickness_D += 1
        elif anchor_D < thickness_D - spring:
            thickness_D -= 1

        # Absorbing boundaries
        if thickness_L <= min_thick or thickness_L >= max_thick:
            if verbose:
                print(f"[L Raft] stopped at {t}")
            break
        if thickness_D <= min_thick or thickness_D >= max_thick:
            if verbose:
                print(f"[D Raft] stopped at {t}")
            break

        hist_L.append(thickness_L)
        hist_D.append(thickness_D)

    return hist_L, hist_D


def plot_coupled_walks(hist_L, hist_D):
    plt.plot(hist_L, label="L-Raft", linewidth=2)
    plt.plot(hist_D, label="D-Raft", linewidth=2, linestyle='--')
    plt.xlabel("Cycle")
    plt.ylabel("Thickness")
    plt.title("Coupled L/D Membrane Thickness Walk")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    hL, hD = run_coupled_walks(k_spring=0.4, verbose=True)
    plot_coupled_walks(hL, hD)

