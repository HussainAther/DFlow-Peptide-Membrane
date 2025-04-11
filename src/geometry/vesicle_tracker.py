import math

class VesicleTracker:
    def __init__(self, initial_radius_nm=50):
        self.radius = initial_radius_nm  # in nanometers
        self.membrane_thickness_nm = 12
        self.area = self.calculate_surface_area()
        self.volume = self.calculate_volume()

    def calculate_surface_area(self):
        return 4 * math.pi * (self.radius ** 2)  # nm²

    def calculate_volume(self):
        return (4 / 3) * math.pi * (self.radius ** 3)  # nm³

    def update_thickness(self, new_thickness):
        self.membrane_thickness_nm = new_thickness

    def update_radius_due_to_insertion(self, num_peptides, peptide_area_nm2=1.5):
        added_area = num_peptides * peptide_area_nm2
        current_area = self.calculate_surface_area()
        new_area = current_area + added_area
        self.radius = math.sqrt(new_area / (4 * math.pi))
        self.area = new_area
        self.volume = self.calculate_volume()

