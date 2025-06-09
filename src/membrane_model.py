import numpy as np
import json
import os

class MembraneModel:
    def __init__(self, initial_thickness=3.0, thickness_variation=0.1):
        """
        Initializes the membrane model with a starting thickness and variation range.
        """
        self.thickness = initial_thickness
        self.thickness_variation = thickness_variation
        self.history = []

    def check_peptide_stability(self, peptide):
        """
        Determines if a peptide contributes to membrane stability.
        Stability is determined by hydrophobic matching.
        """
        hydrophobic_count = sum(1 for aa in peptide if aa in "LR")
        return abs(hydrophobic_count - self.thickness) < 0.5  # Threshold for stability

    def log_thickness(self):
        with open("data/membrane_thickness.json", "w") as f:
            f.write(f'{{"final_thickness": {self.thickness}}}')

    def update_thickness(self, peptides):
        """
        Updates the membrane thickness based on peptide interactions.
        Positive feedback loop: stabilizing peptides reinforce the current thickness.
        """
        stabilizing_peptides = [p for p in peptides if self.check_stability(p)]
        if len(stabilizing_peptides) > len(peptides) / 2:  # Majority rule
            self.thickness += np.random.uniform(-self.thickness_variation, self.thickness_variation)
        self.history.append(self.thickness)

    def save_thickness_data(self, file_path="data/membrane_thickness.json"):
        """
        Saves the membrane thickness evolution to a JSON file.
        """
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "w") as f:
            json.dump({"Membrane_Thickness": self.thickness, "History": self.history}, f, indent=4)

    def get_current_thickness(self):
        """
        Returns the current membrane thickness.
        """
        return self.thickness

