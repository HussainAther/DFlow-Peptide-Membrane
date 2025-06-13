import numpy as np
import json
import os

class MembraneModel:
    def __init__(self, initial_thickness=3.0, thickness_variation=0.1):
        self.thickness = initial_thickness
        self.thickness_variation = thickness_variation
        self.history = []

    def check_peptide_stability(self, peptide):
        """
        Basic filtering: must not be too short, too polar, or too hydrophilic.
        Peptide is a dict with keys: 'sequence', 'chirality', 'charge' (optional)
        """
        seq = peptide.get("sequence", "")
        if len(seq) < 10:
            return False

        charge = peptide.get("charge", 0)
        if abs(charge) > 3:
            return False  # Too polar

        hydrophobic_aas = set("AILMFWVY")  # Typical hydrophobic residues
        hydrophobic_count = sum(seq.count(aa) for aa in hydrophobic_aas)
        if hydrophobic_count / len(seq) < 0.3:
            return False  # Too hydrophilic

        chirality = peptide.get("chirality", "L")
        if chirality not in ["L", "D"]:
            return False

        return True

    def update_thickness(self, peptides):
        if not peptides:
            self.thickness = 3.0
            return

        L_count = sum(1 for p in peptides if p.get("chirality") == "L")
        D_count = sum(1 for p in peptides if p.get("chirality") == "D")
        total = L_count + D_count

        ratio = (L_count - D_count) / total if total > 0 else 0
        self.thickness = 3.0 + 0.2 * ratio  # Between 2.8 and 3.2

    def get_current_thickness(self):
        return self.thickness

    def log_thickness(self, filename="data/membrane_thickness.json"):
        with open(filename, "w") as f:
            json.dump({"thickness": self.thickness}, f)

