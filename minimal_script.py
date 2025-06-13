import random
import matplotlib.pyplot as plt

# Amino acids with simple hydrophobicity scale
amino_acids = {
    'A': 1.8, 'V': 4.2, 'L': 3.8, 'I': 4.5, 'M': 1.9,
    'F': 2.8, 'Y': -1.3, 'W': -0.9, 'S': -0.8, 'T': -0.7,
    'N': -3.5, 'Q': -3.5, 'D': -3.5, 'E': -3.5, 'K': -3.9,
    'R': -4.5, 'H': -3.2, 'G': -0.4, 'P': -1.6, 'C': 2.5
}

def generate_peptide(length=15):
    seq = ''.join(random.choices(list(amino_acids.keys()), k=length))
    chirality = random.choice(['L', 'D'])
    hydro = sum(amino_acids[a] for a in seq) / length
    return {'sequence': seq, 'chirality': chirality, 'hydrophobicity': hydro}

def is_selected(peptide):
    return peptide['chirality'] == 'L' and 1.0 <= peptide['hydrophobicity'] <= 3.0

# Simulation
num_cycles = 50
pool_size = 100
history_L = []
history_D = []

peptides = [generate_peptide() for _ in range(pool_size)]

for cycle in range(num_cycles):
    # Apply selection
    peptides = [p for p in peptides if is_selected(p)]

    # Regenerate new randoms to maintain population size
    while len(peptides) < pool_size:
        peptides.append(generate_peptide())

    # Log counts
    L_count = sum(1 for p in peptides if p['chirality'] == 'L')
    D_count = pool_size - L_count
    history_L.append(L_count)
    history_D.append(D_count)

# Plotting
plt.plot(history_L, label='L peptides', color='blue')
plt.plot(history_D, label='D peptides', color='red')
plt.title("Chirality Selection Over Time")
plt.xlabel("Cycle")
plt.ylabel("Count")
plt.legend()
plt.tight_layout()
plt.savefig("simple_peptide_selection.png", dpi=150)
print("Plot saved as simple_peptide_selection.png")

