import random

# Parameters
PEPTIDE_LENGTH = 80
NUM_PEPTIDES = 10000
INITIAL_MEMBRANE_THICKNESS = 12  # in arbitrary units (e.g., Å)
MAX_MEMBRANE_THICKNESS = 23
GLYCINE_FRACTION = 0.01
ENANTIOMERIC_EXCESS = 0.50  # 0.5 means racemic
MIN_BLOCK_TO_EXPAND = 1.05  # how much longer than current thickness

# Chirality fractions
L_fraction = (1 - GLYCINE_FRACTION) * ENANTIOMERIC_EXCESS
D_fraction = (1 - GLYCINE_FRACTION) * (1 - ENANTIOMERIC_EXCESS)
G_fraction = GLYCINE_FRACTION

CHIRAL_POOL = ['L', 'D', '0']
CHIRAL_WEIGHTS = [L_fraction, D_fraction, G_fraction]

# Tracking
membrane_thickness = INITIAL_MEMBRANE_THICKNESS
membrane_growth_log = []
accepted_peptides = []

# --- Utility ---
def generate_ld0_peptide(length=PEPTIDE_LENGTH):
    return [random.choices(CHIRAL_POOL, weights=CHIRAL_WEIGHTS, k=1)[0] for _ in range(length)]

def longest_chiral_block(peptide, target_chirality):
    count = 0
    max_count = 0
    for aa in peptide:
        if aa == target_chirality:
            count += 1
            max_count = max(max_count, count)
        else:
            count = 0
    return max_count

# --- Main Simulation Loop ---
for i in range(NUM_PEPTIDES):
    peptide = generate_ld0_peptide()
    l_block = longest_chiral_block(peptide, 'L')
    d_block = longest_chiral_block(peptide, 'D')
    max_block = max(l_block, d_block)

    if max_block >= membrane_thickness * MIN_BLOCK_TO_EXPAND:
        membrane_thickness = min(membrane_thickness + 1, MAX_MEMBRANE_THICKNESS)
        accepted_peptides.append({
            'index': i,
            'peptide': ''.join(peptide),
            'block': max_block,
            'updated_thickness': membrane_thickness
        })

    membrane_growth_log.append(membrane_thickness)

print(f"[✔] Simulation complete. Final membrane thickness: {membrane_thickness}")
print(f"Accepted peptides: {len(accepted_peptides)} / {NUM_PEPTIDES}")

