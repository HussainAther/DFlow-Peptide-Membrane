import os
import torch
import numpy as np
import pandas as pd
from src.peptide_generator import PeptideGenerator
from src.membrane_model import MembraneModel
from src.visualization import Visualization

# Ensure data directory exists
os.makedirs("data", exist_ok=True)

def main():
    num_cycles = 500  # Number of selection cycles
    num_peptides = 100  # Initial peptide pool size
    sequence_length = 15  # Length of each peptide sequence
    
    # Initialize components
    generator = PeptideGenerator(num_peptides, sequence_length)
    membrane = MembraneModel()
    peptides = generator.generate_peptides()
    
    evolution_log = []
    
    for cycle in range(num_cycles):
        peptides = [p for p in peptides if membrane.check_peptide_stability(p)]
        evolution_log.append({
            "Cycle": cycle,
            "Peptide_Count": len(peptides),
            "Membrane_Thickness": membrane.thickness
        })
        
        if len(peptides) < 5:
            peptides += generator.generate_peptides()
        
        membrane.update_thickness(peptides)
    
    # Save logs
    df = pd.DataFrame(evolution_log)
    df.to_csv("data/peptide_evolution.csv", index=False)
    membrane.log_thickness()
    
    # Visualization
    Visualization.plot_membrane_thickness()
    Visualization.plot_peptide_selection()

if __name__ == "__main__":
    main()

