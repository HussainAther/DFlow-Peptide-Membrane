import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving only

from src.alpha_helix_raft_simulation import RaftSimulation
from src.peptide_generator import PeptideGenerator
from src.membrane_model import MembraneModel
from src.visualization import Visualization
from src.alpha_helix_raft_simulation import RaftSimulation

# Ensure data directory exists
os.makedirs("data", exist_ok=True)

def main():
    num_cycles = 500
    num_peptides = 100
    sequence_length = 15
    membrane = MembraneModel()
    enable_catalysis = True

    # Initialize components
    generator = PeptideGenerator(num_peptides, sequence_length)
    raft_sim = RaftSimulation(grid_height=8, grid_width=12, max_rafts=20, enable_catalysis=enable_catalysis)
    raft_sim = RaftSimulation(grid_height=8, grid_width=12)
    if cycle == 0:
        raft_sim.place_peptides()
    peptides = generator.generate_peptides()

    evolution_log = []

    for cycle in range(num_cycles):
        raft_sim.simulate_step()  # Run raft movement + desorption

        # Log cycle info
        evolution_log.append({
            "Cycle": cycle,
            "Peptide_Count": len(peptides),
            "Membrane_Thickness": membrane.thickness,
            "Raft_Count": raft_sim.count_rafts(),
            "L_fraction": raft_sim.chirality_fraction('L'),
            "D_fraction": raft_sim.chirality_fraction('D'),
            "With_Catalysis": enable_catalysis
        })

        # Optional: regenerate peptide pool if needed
        if len(peptides) < 5:
            peptides += generator.generate_peptides()

    # Save logs and visualizations
    df = pd.DataFrame(evolution_log)
    df.to_csv("data/peptide_evolution.csv", index=False)
    catalysis_tag = "with_catalysis" if enable_catalysis else "no_catalysis"
    df.to_csv(f"data/peptide_evolution_{catalysis_tag}.csv", index=False)
    raft_sim.save_raft_distribution(f"data/raft_final_{catalysis_tag}.png")
    Visualization.plot_membrane_thickness(save_path="data/thickness_plot.png")
    Visualization.plot_peptide_selection(save_path="data/peptide_selection.png")


if __name__ == "__main__":
    main()

