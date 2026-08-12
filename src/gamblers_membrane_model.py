import random
import matplotlib.pyplot as plt

class Peptide:
    def __init__(self, anchor_len):
        self.anchor_len = anchor_len

class MembraneRaft:
    def __init__(self, start_thickness=12, min_thickness=10, max_thickness=23):
        self.thickness = start_thickness
        self.min_thickness = min_thickness
        self.max_thickness = max_thickness
        self.history = [start_thickness]

    def update_thickness(self, peptide: Peptide):
        delta = peptide.anchor_len - self.thickness

        if delta > 1:
            self.thickness += 1
        elif delta < -1:
            self.thickness -= 1

        # Clamp
        self.thickness = max(self.min_thickness, min(self.thickness, self.max_thickness))
        self.history.append(self.thickness)

    def is_absorbed(self):
        return self.thickness == self.min_thickness or self.thickness == self.max_thickness


def generate_random_peptide(min_len=6, max_len=20):
    return Peptide(anchor_len=random.randint(min_len, max_len))


def run_random_walk(max_cycles=1000):
    raft = MembraneRaft()
    cycle = 0
    while not raft.is_absorbed() and cycle < max_cycles:
        peptide = generate_random_peptide()
        raft.update_thickness(peptide)
        cycle += 1
    return raft.history


def plot_membrane_walk(history):
    plt.plot(history, linewidth=2)
    plt.xlabel("Cycle")
    plt.ylabel("Membrane Thickness (nm)")
    plt.title("Membrane Thickness Random Walk (Anchor-Length Driven)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    walk = run_random_walk()
    plot_membrane_walk(walk)

