import argparse
import torch
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
from dflow.model import PeptideFlow  # Assuming D-Flow model exists in repo

parser = argparse.ArgumentParser()
parser.add_argument("--init", type=str, default="default", choices=["default", "space"],
                    help="Initialization mode: 'default' for random L-peptides, 'space' for prebiotic racemic amino acids")
args = parser.parse_args()

space_presets = {
    "carbonaceous_chondrite": {
        'A': 0.12, 'D': 0.08, 'E': 0.07, 'G': 0.15, 'V': 0.10,
        'L': 0.07, 'I': 0.05, 'P': 0.05, 'S': 0.08, 'T': 0.07,
        'others': 0.16  # rare or less detected ones
    },
    "irradiated_meteorite": {
        # Bias toward UV-resistant or small, stable amino acids
        'G': 0.20, 'A': 0.18, 'V': 0.14, 'D': 0.10, 'others': 0.38
    },
    "vesicle_shadowed_pool": {
        # Favors hydrophobic & amphipathic peptides
        'V': 0.18, 'L': 0.16, 'A': 0.14, 'G': 0.12, 'P': 0.10, 'others': 0.30
    }
}


class DFlowPeptideMembraneSim:
    def __init__(self, num_peptides=100, sequence_length=15):
        self.num_peptides = num_peptides
        self.sequence_length = sequence_length
        self.peptides = self.generate_peptides()
        self.membrane_thickness = np.random.uniform(2.5, 3.5)
        self.peptide_evolution_log = []
    
    def generate_peptides(self):
        """Use D-Flow to generate hybrid L-D peptides."""
        model = PeptideFlow.load_pretrained("dflow_model.pt")  # Load trained D-Flow model
        l_d_peptides = model.sample(self.num_peptides, seq_length=self.sequence_length, hybrid=True)
        return l_d_peptides
    
    def check_stability(self, peptide):
        """Determine if peptide matches membrane hydrophobic thickness."""
        hydrophobic_count = sum(1 for aa in peptide if aa in "LR")
        return abs(hydrophobic_count - self.membrane_thickness) < 0.5
    
    def iterate_selection(self, cycles=100):
        """Simulate peptide selection over multiple evolutionary cycles."""
        for cycle in range(cycles):
            self.peptides = [p for p in self.peptides if self.check_stability(p)]
            self.peptide_evolution_log.append({
                "Cycle": cycle,
                "Peptide_Count": len(self.peptides),
                "Membrane_Thickness": self.membrane_thickness
            })
            if len(self.peptides) < 5:  # Maintain genetic diversity
                self.peptides += self.generate_peptides()
            self.membrane_thickness += np.random.uniform(-0.1, 0.1)  # Simulate feedback loop
    
    def save_logs(self):
        """Save peptide evolution log and membrane thickness data."""
        df = pd.DataFrame(self.peptide_evolution_log)
        df.to_csv("data/peptide_evolution.csv", index=False)
        
        with open("data/membrane_thickness.json", "w") as f:
            json.dump({"Membrane_Thickness": self.membrane_thickness}, f)
    
    def visualize_results(self):
        """Plot membrane thickness over time."""
        plt.figure(figsize=(8, 5))
        plt.hist([sum(1 for aa in p if aa in "LR") for p in self.peptides], bins=10, color='blue', alpha=0.7)
        plt.axvline(self.membrane_thickness, color='red', linestyle='dashed', label='Membrane Thickness')
        plt.xlabel("Peptide Hydrophobicity Score")
        plt.ylabel("Frequency")
        plt.title("Peptide Selection Over Time")
        plt.legend()
        plt.show()

if __name__ == "__main__":
    sim = DFlowPeptideMembraneSim()
    sim.iterate_selection(cycles=500)
    sim.save_logs()
    sim.visualize_results()

