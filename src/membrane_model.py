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
        Filters peptide based on basic properties like length, charge, and chirality.
        Peptide is a dict with keys: 'sequence', 'chirality', 'charge' (optional).
        """
        seq = peptide.get("sequence", "")
        if len(seq) < 10 or len(seq) > 25:
            return False  # filter extreme lengths

        charge = peptide.get("charge", 0)
        if abs(charge) > 3:
            return False  # extreme polarity rejected

        chirality = peptide.get("chirality", "L")
        if chirality not in ["L", "D"]:
            return False

        return True

    def update_thickness(self, peptides):
        """
        Updates membrane thickness based on peptide chirality distribution.
        L-rich peptides stabilize (increase) thickness, D-rich destabilize.
        """
        if not peptides:
            self.thickness = 3.0
            return

        L_count = sum(1 for p in peptides if p.get("chirality", "L") == "L")
        D_count = sum(1 for p in peptides if p.get("chirality", "L") == "D")
        total = L_count + D_count

        if total == 0:
            self.thickness = 3.0
            return

        ratio = (L_count - D_count) / total
        variation = np.random.normal(0, self.thickness_variation / 2)
        self.thickness = 3.0 + 0.2 * ratio + variation
        self.thickness = max(2.7, min(3.3, self.thickness))  # clamp within biological range
        self.history.append(self.thickness)

    def get_current_thickness(self):
        return self.thickness

    def log_thickness(self, filename="data/membrane_thickness.json"):
        """
        Save thickness history to a JSON file for further inspection.
        """
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, "w") as f:
            json.dump(self.history, f, indent=2)

