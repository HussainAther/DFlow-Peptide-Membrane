import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def calculate_snapshot_lateral_rdf(struct_path, peptide_resname, r_max=2.0, bin_width=0.02):
    """
    Calculates the 2D (XY-plane) Radial Distribution Function for a static coordinate structure.
    struct_path: Path to GROMACS structure file (.gro)
    peptide_resname: Resname of the target peptide (e.g., 'LALA' or 'LDALA')
    """
    # Load universe using a single frame
    u = mda.Universe(struct_path)
    
    # Select backbone beads
    bb_selection = u.select_atoms(f"resname {peptide_resname} and name BB")
    num_particles = len(bb_selection)
    
    if num_particles == 0:
        print(f"No backbone beads found for resname {peptide_resname}")
        return None, None

    # Set up distance bins
    bins = np.arange(0, r_max + bin_width, bin_width)
    rdf_counts = np.zeros(len(bins) - 1)
    
    box_dimensions = u.dimensions  # [Lx, Ly, Lz, alpha, beta, gamma]
    area = box_dimensions[0] * box_dimensions[1]
    density = num_particles / area
    
    # Extract XY positions
    positions = bb_selection.positions[:, :2]
    
    # Compute pairwise distances
    for i in range(num_particles):
        diff = positions - positions[i]
        
        # Apply 2D Periodic Boundary Conditions (PBC)
        diff[:, 0] -= box_dimensions[0] * np.round(diff[:, 0] / box_dimensions[0])
        diff[:, 1] -= box_dimensions[1] * np.round(diff[:, 1] / box_dimensions[1])
        
        distances = np.linalg.norm(diff, axis=1)
        distances = distances[distances > 0.0]
        
        counts, _ = np.histogram(distances, bins=bins)
        rdf_counts += counts

    # Normalize the histogram
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    bin_centers_nm = bin_centers / 10.0 # Convert Angstroms to nm
    
    for i in range(len(bin_centers)):
        ring_area = np.pi * (bins[i+1]**2 - bins[i]**2)
        expected_particles = num_particles * density * ring_area
        rdf_counts[i] /= expected_particles

    print(f"Successfully evaluated packing coordination for {peptide_resname}.")
    return bin_centers_nm, rdf_counts

# --- Execution Block ---
# Parse the static configuration layout directly from your production frame
bin_centers_l, rdf_l = calculate_snapshot_lateral_rdf("nvt_highres.gro", "LALA")
bin_centers_ld, rdf_ld = calculate_snapshot_lateral_rdf("nvt_highres.gro", "LDALA")

# Generate the comparison plot
plt.figure(figsize=(6, 4))
if bin_centers_l is not None:
    plt.plot(bin_centers_l, rdf_l, label='Pure L (Pleated)', color='blue')
if bin_centers_ld is not None:
    plt.plot(bin_centers_ld, rdf_ld, label='L/D Mix (Rippled)', color='red')

plt.xlabel('Distance r (nm)')
plt.ylabel('g(r)')
plt.title('Lateral Inter-Peptide Packing Coordination')
plt.legend()
plt.tight_layout()
plt.savefig('lateral_clustering_rdf.png', dpi=300)
print("Analysis complete. Comparison plot saved as lateral_clustering_rdf.png")
