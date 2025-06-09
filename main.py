import os
import torch
import numpy as np
import pandas as pd
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
    
    # Initialize components
    generator = PeptideGenerator(num_peptides, sequence_length)
    membrane = MembraneModel()
    raft_sim = RaftSimulation(grid_height=8, grid_width=12)
    raft_sim.simulate_step()
    peptides = generator.generate_peptides()

    evolution_log = []

    for cycle in range(num_cycles):
        # Filter peptides by stability (possibly chirality-aware)
        peptides = [p for p in peptides if membrane.check_peptide_stability(p)]

        # Run raft update with peptide placement, diffusion, desorption
        raft_sim.place_peptides(peptides)
        raft_sim.diffuse_rafts()
        raft_sim.desorb_unanchored()

        # Log cycle info
        evolution_log.append({
            "Cycle": cycle,
            "Peptide_Count": len(peptides),
            "Membrane_Thickness": membrane.thickness,
            "Raft_Count": raft_sim.count_rafts(),
            "L_fraction": raft_sim.chirality_fraction('L'),
            "D_fraction": raft_sim.chirality_fraction('D'),
        })

        # Regenerate if population is too small
        if len(peptides) < 5:
            peptides += generator.generate_peptides()
        
        # Update global membrane properties
        membrane.update_thickness(peptides)

    # Save and visualize
    df = pd.DataFrame(evolution_log)
    df.to_csv("data/peptide_evolution.csv", index=False)
    membrane.log_thickness()
    Visualization.plot_membrane_thickness()
    Visualization.plot_peptide_selection()
    raft_sim.save_raft_distribution("data/raft_final.png")

if __name__ == "__main__":
    main()

