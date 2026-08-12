import random
import csv
import json
from pathlib import Path

# ==================== CONFIGURATION ====================

OUTPUT_DIR = Path("logs")
OUTPUT_DIR.mkdir(exist_ok=True)

NUM_PEPTIDES = 10000
GLYCINE_FRACTION = 0.01
ENANTIOMERIC_EXCESS = 0.5  # 0.5 = racemic
INITIAL_MEMBRANE_THICKNESS = 12
MAX_MEMBRANE_THICKNESS = 23
BLOCK_GROWTH_THRESHOLD = 1.05
VESICLE_DIAMETER = 100  # nm
MAX_PEPTIDE_LENGTH = 120
CYTOPLASM_CAPACITY = 1000  # number of peptides
MAX_PEPTIDES_PER_AREA = 0.5  # peptides per nm² membrane area

CHIRAL_POOL = ['L', 'D', '0']
L_fraction = (1 - GLYCINE_FRACTION) * ENANTIOMERIC_EXCESS
D_fraction = (1 - GLYCINE_FRACTION) * (1 - ENANTIOMERIC_EXCESS)
G_fraction = GLYCINE_FRACTION
CHIRAL_WEIGHTS = [L_fraction, D_fraction, G_fraction]

# ==================== STATE VARIABLES ====================

membrane_thickness = INITIAL_MEMBRANE_THICKNESS
inserted_peptides = 0
cytoplasmic_peptides = 0
discarded_peptides = 0
cytoplasm_full = False

membrane_growth_log = []
peptide_logs = []
vesicle_state_log = []

# Membrane area (approximate sphere surface)
membrane_area = 4 * 3.1415 * (VESICLE_DIAMETER / 2) ** 2
max_peptides_membrane = int(membrane_area * MAX_PEPTIDES_PER_AREA)

# ==================== HELPER FUNCTIONS ====================

def longest_chiral_block(seq, target):
    count = 0
    max_count = 0
    for ch in seq:
        if ch == target:
            count += 1
            max_count = max(max_count, count)
        else:
            count = 0
    return max_count

def generate_peptide():
    return []

# ==================== MAIN SIMULATION ====================

for i in range(NUM_PEPTIDES):
    peptide = []
    inserted = False
    max_block = 0

    for _ in range(MAX_PEPTIDE_LENGTH):
        aa = random.choices(CHIRAL_POOL, weights=CHIRAL_WEIGHTS, k=1)[0]
        peptide.append(aa)

        l_block = longest_chiral_block(peptide, 'L')
        d_block = longest_chiral_block(peptide, 'D')
        max_block = max(l_block, d_block)

        if max_block >= membrane_thickness * BLOCK_GROWTH_THRESHOLD:
            # Peptide is accepted into membrane
            if inserted_peptides < max_peptides_membrane:
                inserted = True
                inserted_peptides += 1
                if membrane_thickness < MAX_MEMBRANE_THICKNESS:
                    membrane_thickness += 1
            break

    # If not inserted
    if not inserted:
        if len(peptide) >= MAX_PEPTIDE_LENGTH:
            if not cytoplasm_full:
                cytoplasmic_peptides += 1
                if cytoplasmic_peptides >= CYTOPLASM_CAPACITY:
                    cytoplasm_full = True
                fate = "cytoplasm"
            else:
                discarded_peptides += 1
                fate = "discarded"
        else:
            # incomplete peptide (shouldn't happen here)
            discarded_peptides += 1
            fate = "discarded"
    else:
        fate = "membrane"

    # Log this peptide
    peptide_logs.append({
        "index": i,
        "sequence": '-'.join(peptide),
        "length": len(peptide),
        "max_block": max_block,
        "inserted": inserted,
        "location": fate
    })

    membrane_growth_log.append({
        "step": i,
        "thickness_nm": membrane_thickness
    })

    vesicle_state_log.append({
        "step": i,
        "membrane_thickness": membrane_thickness,
        "peptides_membrane": inserted_peptides,
        "peptides_cytoplasm": cytoplasmic_peptides,
        "peptides_discarded": discarded_peptides,
        "cytoplasm_full": cytoplasm_full
    })

# ==================== EXPORT LOGS ====================

def write_csv(filename, fieldnames, rows):
    with open(OUTPUT_DIR / filename, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

write_csv("membrane_growth_log.csv", ["step", "thickness_nm"], membrane_growth_log)
write_csv("peptide_log.csv", ["index", "sequence", "length", "max_block", "inserted", "location"], peptide_logs)
write_csv("vesicle_state_log.csv", ["step", "membrane_thickness", "peptides_membrane",
                                    "peptides_cytoplasm", "peptides_discarded", "cytoplasm_full"], vesicle_state_log)

# Final state JSON
final_state = {
    "final_membrane_thickness": membrane_thickness,
    "total_peptides_tested": NUM_PEPTIDES,
    "inserted": inserted_peptides,
    "cytoplasm": cytoplasmic_peptides,
    "discarded": discarded_peptides,
    "cytoplasm_full": cytoplasm_full
}

with open(OUTPUT_DIR / "final_state.json", "w") as f:
    json.dump(final_state, f, indent=4)

print(f"[✔] Simulation complete.")
print(f"Membrane thickness: {membrane_thickness}")
print(f"Inserted: {inserted_peptides}, Cytoplasm: {cytoplasmic_peptides}, Discarded: {discarded_peptides}")

