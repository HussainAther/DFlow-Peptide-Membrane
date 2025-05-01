# src/vesicle_geometry.py

import math

class Vesicle:
    def __init__(self, initial_radius_nm=25, peptide_area_nm2=1.5, flattened=False):
        self.radius_nm = initial_radius_nm  # initial radius (e.g. 25 nm ~ small vesicle)
        self.peptide_area = peptide_area_nm2  # estimated area per peptide
        self.flattened = flattened
        self.inserted_peptides = 0

    def surface_area(self):
        if self.flattened:
            return 2 * math.pi * (self.radius_nm ** 2)  # approx. double surface
        return 4 * math.pi * (self.radius_nm ** 2)

    def volume(self):
        return (4/3) * math.pi * (self.radius_nm ** 3)

    def capacity(self):
        """Max number of peptides that can insert given surface area."""
        return int(self.surface_area() / self.peptide_area)

    def insert_peptide(self):
        self.inserted_peptides += 1
        if self.inserted_peptides >= self.capacity():
            self.expand()

    def expand(self, growth_step_nm=5):
        """Increase vesicle size by step."""
        self.radius_nm += growth_step_nm
        self.inserted_peptides = 0  # reset peptide count
        print(f"🫧 Vesicle expanded → new radius: {self.radius_nm:.1f} nm")

    def status(self):
        return {
            'Radius_nm': self.radius_nm,
            'Volume_nm3': self.volume(),
            'SurfaceArea_nm2': self.surface_area(),
            'PeptideCount': self.inserted_peptides,
            'Capacity': self.capacity()
        }

