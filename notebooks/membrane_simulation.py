import random
import csv
from pathlib import Path

# -------------------- CONFIGURATION --------------------
GLY_FRACTION = 0.01          # Probability of glycine (achiral)
ANCHOR_LEN = 5               # Length of LLLLL or DDDDD to anchor
NUM_CYCLES = 1000            # How many peptides to generate
START_MEMBRANE_THICKNESS = 12
MISMATCH_THRESHOLD = 5       # Number of mismatched peptides before growth
MAX_MEMBRANE_THICKNESS = 25

# -------------------- INITIALIZATION --------------------
membrane_thickness = START_MEMBRANE_THICKNESS
L_raft = []
D_raft = []
log_data = []

def generate_peptide():
    probabilities = [GLY_FRACTION, (1 - GLY_FRACTION) / 2, (1 - GLY_FRACTION) / 2]
    residues = ['0', 'D', 'L']
    seq = ''
    while True:
        seq += random.choices(residues, probabilities)[0]
        if 'D' * ANCHOR_LEN in seq or 'L' * ANCHOR_LEN in seq:
            break
    return seq

def count_mismatches(raft):
    return sum(1 for p in raft if terminal_block_length(p) > membrane_thickness)

def terminal_block_length(peptide):
    last_residue = peptide[-1]
    count = 0
    for res in reversed(peptide):
        if res == last_residue:
            count += 1
        else:
            break
    return count

# -------------------- SIMULATION LOOP --------------------
for cycle in range(1, NUM_CYCLES + 1):
    peptide = generate_peptide()
    anchor_type = 'L' if peptide.endswith('L' * ANCHOR_LEN) else 'D'

    if anchor_type == 'L':
        L_raft.append(peptide)
        mismatches = count_mismatches(L_raft)
    else:
        D_raft.append(peptide)
        mismatches = count_mismatches(D_raft)

    # Check for membrane growth
    if mismatches >= MISMATCH_THRESHOLD and membrane_thickness < MAX_MEMBRANE_THICKNESS:
        membrane_thickness += 1
        print(f"[Cycle {cycle}] Membrane grew to thickness {membrane_thickness}")

    # Log state
    log_data.append({
        'Cycle': cycle,
        'MembraneThickness': membrane_thickness,
        'L_Raft_Size': len(L_raft),
        'D_Raft_Size': len(D_raft),
        'LastPeptideLength': len(peptide),
        'AnchorType': anchor_type
    })

# -------------------- EXPORT LOG --------------------
output_dir = Path("experiments/phase1_run")
output_dir.mkdir(parents=True, exist_ok=True)

with open(output_dir / "membrane_log.csv", 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=log_data[0].keys())
    writer.writeheader()
    writer.writerows(log_data)

print(f"\n✅ Simulation complete. Log saved to {output_dir/'membrane_log.csv'}")

