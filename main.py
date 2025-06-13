import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Save plots without displaying

from src.alpha_helix_raft_simulation import RaftSimulation
from src.peptide_generator import PeptideGenerator
from src.membrane_model import MembraneModel
from src.visualization import Visualization

os.makedirs("data", exist_ok=True)

def main():
    num_cycles = 500
    num_peptides = 100
    sequence_length = 15
    enable_catalysis = True

    membrane = MembraneModel()
    generator = PeptideGenerator(num_peptides, sequence_length)
    raft_sim = RaftSimulation(grid_height=8, grid_width=12, max_rafts=20, enable_catalysis=enable_catalysis)
    
    raft_sim.place_peptides()
    peptides = generator.generate_peptides()

    evolution_log = []

    for cycle in range(num_cycles):
        stable_peptides = []
        rejected_count = 0

        for p in peptides:
            if membrane.check_peptide_stability(p):
                stable_peptides.append(p)
            else:
                rejected_count += 1

        peptides = stable_peptides

        raft_sim.simulate_step()

        membrane.update_thickness(peptides)

        evolution_log.append({
            "Cycle": cycle,
            "Peptide_Count": len(peptides),
            "Rejected": rejected_count,
            "Membrane_Thickness": membrane.thickness,
            "Raft_Count": raft_sim.count_rafts(),
            "L_fraction": raft_sim.chirality_fraction('L'),
            "D_fraction": raft_sim.chirality_fraction('D'),
            "With_Catalysis": enable_catalysis
        })

        if len(peptides) < 5:
            peptides += generator.generate_peptides()

    df = pd.DataFrame(evolution_log)
    catalysis_tag = "with_catalysis" if enable_catalysis else "no_catalysis"
    df.to_csv(f"data/peptide_evolution_{catalysis_tag}.csv", index=False)

    Visualization.plot_membrane_thickness(save_path="data/thickness_plot.png")
    Visualization.plot_peptide_selection(save_path="data/peptide_selection.png")
    raft_sim.save_raft_distribution(f"data/raft_final_{catalysis_tag}.png")

if __name__ == "__main__":
    main()

