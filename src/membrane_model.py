import numpy as np
import json
import os


# src/membrane_model.py

class MembraneModel:
    def __init__(self, initial_thickness=3.0, thickness_variation=0.1):
        self.thickness = initial_thickness
        self.thickness_variation = thickness_variation
        self.history = []

    def check_peptide_stability(self, peptide):
        seq = peptide.get("sequence", "")
        if len(seq) < 10:
            return False

        charge = peptide.get("charge", 0)
        if abs(charge) > 4:
            return False  # too polar

        chirality = peptide.get("chirality", "L")
        if chirality not in ["L", "D"]:
            return False

        # Hydrophobicity match: if membrane is thick (L-rich), reject hydrophilic
        hydrophobicity = peptide.get("hydrophobicity", 0)
        if self.thickness > 3.0 and hydrophobicity < 0:
            return False
        elif self.thickness < 3.0 and hydrophobicity > 0:
            return False

        return True

    def update_thickness(self, peptides):
        if not peptides:
            self.thickness = 3.0
            return

        L_count = sum(1 for p in peptides if p.get("chirality", "L") == "L")
        D_count = sum(1 for p in peptides if p.get("chirality", "L") == "D")
        total = L_count + D_count
        ratio = (L_count - D_count) / total if total > 0 else 0

        # Stronger feedback for visual change
        self.thickness = 3.0 + 0.4 * ratio  # between 2.6 and 3.4
        self.thickness = max(2.5, min(self.thickness, 3.5))
        self.history.append(self.thickness)

    def get_current_thickness(self):
        return self.thickness

    def log_thickness(self, filename="data/membrane_thickness.json"):
        with open(filename, "w") as f:
            json.dump({"thickness": self.thickness}, f)

