import random
import csv
import json
from pathlib import Path

# ==================== CONFIGURATION ====================

OUTPUT_DIR = Path("logs_v2")
OUTPUT_DIR.mkdir(exist_ok=True)

NUM_PEPTIDES = 10000
GLYCINE_FRACTION = 0.01
ENANTIOMERIC_EXCESS = 0.5
INITIAL_MEMBRANE_THICKNESS = 12
MAX_MEMBRANE_THICKNESS = 23
HYDROPHOBIC_N_MARGIN = 2  # required block = thickness + n
MISFITTED_INSERTION_PROB = 0.5
PARTIAL_INSERTION_PROB = 0.25
VESICLE_DIAMETER = 100  # nm
MAX_PEPTIDE_LENGTH = 120
CYTOPLASM_CAPACITY = float('inf')  # effectively unlimited for now
MAX_PEPTIDES_PER_AREA = 0.5  # peptides per nm²

CHIRAL_POOL = ['L', 'D', '0']
L_fraction = (1 - GLYCINE_FRACTION) * ENANTIOMERIC_EXCESS
D_fraction = (1 - GLYCINE_FRACTION) * (1 - ENANTIOMERIC_EXCESS)
G_fraction = GLYCINE_FRACTION
CHIRAL_WEIGHTS = [L_fraction, D_fraction, G_fraction]

AA_HYDROPHOBICITY = {
    'A': 1.8,  'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8,  'K': -3.9, 'M': 1.9,  'F': 2.8,  'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# ==================== STATE VARIABLES ====================

membrane_thickness = INITIAL_MEMBRANE_THICKNESS
inserted_peptides = 0
cytoplasmic_peptides = 0
discarded_peptides = 0
cytoplasm_full = False

raft_L_count = 0
raft_D_count = 0

membrane_growth_log = []
peptide_logs = []
vesicle_state_log = []

membrane_area = 4 * 3.1415 * (VESICLE_DIAMETER / 2) ** 2
max_peptides_membrane = int(membrane_area * MAX_PEPTIDES_PER_AREA)

# ==================== HELPERS ====================

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

def compute_hydrophobicity_score(peptide):
    total = 0
    for aa in peptide:
        base = aa.replace('L-', '').replace('D-', '').replace('0-', '')
        if base in AA_HYDROPHOBICITY:
            total += AA_HYDROPHOBICITY[base]
    return round(total / len(peptide), 3)

# ==================== MAIN LOOP ====================

for i in range(NUM_PEPTIDES):
    peptide = []
    inserted = False
    partial_inserted = False
    max_block = 0
    hydrophobicity = 0
    dominant_chirality = None

    for _ in range(MAX_PEPTIDE_LENGTH):
        chiral = random.choices(CHIRAL_POOL, weights=CHIRAL_WEIGHTS, k=1)[0]
        aa = random.choice(list(AA_HYDROPHOBICITY.keys()))
        peptide.append(f"{chiral}-{aa}")

        l_block = longest_chiral_block(peptide, 'L')
        d_block = longest_chiral_block(peptide, 'D')
        max_block = max(l_block, d_block)
        hydrophobicity = compute_hydrophobicity_score(peptide)

        required_block = membrane_thickness + HYDROPHOBIC_N_MARGIN

        if max_block >= required_block:
            if inserted_peptides < max_peptides_membrane:
                inserted = True
                inserted_peptides += 1
                if l_block > d_block:
                    raft_L_count += 1
                    dominant_chirality = 'L'
                else:
                    raft_D_count += 1
                    dominant_chirality = 'D'
                if membrane_thickness < MAX_MEMBRANE_THICKNESS:
                    membrane_thickness += 1
            break

    if not inserted:
        if max_block >= membrane_thickness:
            if random.random() < PARTIAL_INSERTION_PROB:
                partial_inserted = True
                inserted_peptides += 1
                dominant_chirality = 'L' if l_block >= d_block else 'D'
                if dominant_chirality == 'L':
                    raft_L_count += 1
                else:
                    raft_D_count += 1
            else:
                # Misfitted attempt: chance of cytoplasm
                if random.random() < MISFITTED_INSERTION_PROB:
                    cytoplasmic_peptides += 1
                    fate = "cytoplasm"
                else:
                    discarded_peptides += 1
                    fate = "discarded"
        else:
            discarded_peptides += 1
            fate = "discarded"

    if inserted:
        fate = "membrane"
    elif partial_inserted:
        fate = "partial"

    # Logging peptide
    peptide_logs.append({
        "index": i,
        "sequence": '-'.join(peptide),
        "length": len(peptide),
        "max_block": max_block,
        "inserted": inserted,
        "partial": partial_inserted,
        "location": fate,
        "hydrophobicity": hydrophobicity,
        "dominant_chirality": dominant_chirality
    })

    # Log membrane thickness
    membrane_growth_log.append({
        "step": i,
        "thickness_nm": membrane_thickness
    })

    # Vesicle state tracking
    vesicle_state_log.append({
        "step": i,
        "membrane_thickness": membrane_thickness,
        "peptides_membrane": inserted_peptides,
        "peptides_cytoplasm": cytoplasmic_peptides,
        "peptides_discarded": discarded_peptides,
        "raft_L": raft_L_count,
        "raft_D": raft_D_count
    })

# ==================== EXPORT ====================

def write_csv(name, fields, rows):
    with open(OUTPUT_DIR / name, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

write_csv("membrane_growth_log.csv", ["step", "thickness_nm"], membrane_growth_log)
write_csv("peptide_log.csv", ["index", "sequence", "length", "max_block", "inserted", "partial",
                              "location", "hydrophobicity", "dominant_chirality"], peptide_logs)
write_csv("vesicle_state_log.csv", ["step", "membrane_thickness", "peptides_membrane",
                                    "peptides_cytoplasm", "peptides_discarded",
                                    "raft_L", "raft_D"], vesicle_state_log)

with open(OUTPUT_DIR / "final_state.json", "w") as f:
    json.dump({
        "final_membrane_thickness": membrane_thickness,
        "inserted": inserted_peptides,
        "cytoplasm": cytoplasmic_peptides,
        "discarded": discarded_peptides,
        "raft_L": raft_L_count,
        "raft_D": raft_D_count
    }, f, indent=4)

print(f"[✔] Simulation V2 complete. Membrane thickness: {membrane_thickness} nm")
print(f"Raft L: {raft_L_count}, Raft D: {raft_D_count}")

