import random
import csv
import argparse
from pathlib import Path
from src.uv_filter import is_uv_degraded
from src.vesicle_geometry import Vesicle

# -------------------- ARGPARSE --------------------
parser = argparse.ArgumentParser(description="DL-Peptide Membrane Simulation with UV degradation")
parser.add_argument('gly_fraction', type=float, nargs='?', default=0.01, help="Glycine fraction (0-1)")
parser.add_argument('anchor_len', type=int, nargs='?', default=5, help="Anchor length (e.g., 5 for LLLLL/DDDDD)")
parser.add_argument('mismatch_threshold', type=int, nargs='?', default=5, help="Peptides causing mismatch before growth")
parser.add_argument('num_cycles', type=int, nargs='?', default=1000, help="Number of peptides to generate")
parser.add_argument('output_dir', type=str, nargs='?', default="experiments/manual_run", help="Directory to save logs")
parser.add_argument('--uv_prob', type=float, default=0.7, help="UV degradation probability if D-aromatics exceed limit")
parser.add_argument('--max_d_aromatics', type=int, default=1, help="Maximum D-aromatics allowed before degradation")

args = parser.parse_args()

# -------------------- CONFIG --------------------
GLY_FRACTION = args.gly_fraction
ANCHOR_LEN = args.anchor_len
MISMATCH_THRESHOLD = args.mismatch_threshold
NUM_CYCLES = args.num_cycles
UV_PROB = args.uv_prob
MAX_D_AROMATICS = args.max_d_aromatics

output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

START_MEMBRANE_THICKNESS = 12
MAX_MEMBRANE_THICKNESS = 25

# -------------------- INITIALIZATION --------------------
membrane_thickness = START_MEMBRANE_THICKNESS
L_raft = []
D_raft = []
log_data = []
uv_degraded_count = 0

vesicle = Vesicle(initial_radius_nm=25, flattened=False)
vesicle.insert_peptide()

log_data.append({
    'Cycle': cycle,
    'MembraneThickness': membrane_thickness,
    'L_Raft_Size': len(L_raft),
    'D_Raft_Size': len(D_raft),
    'L_to_D_Ratio': len(L_raft) / max(len(D_raft), 1),
    'VesicleRadius_nm': vesicle.radius_nm,
    'VesicleVolume_nm3': vesicle.volume(),
    'VesicleSurfaceArea_nm2': vesicle.surface_area(),
    'PeptidesInVesicle': vesicle.inserted_peptides
})


# -------------------- FUNCTIONS --------------------
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

    # UV degradation filter
    if is_uv_degraded(peptide, max_d_aromatics=MAX_D_AROMATICS, degradation_prob=UV_PROB):
        uv_degraded_count += 1
        continue

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
with open(output_dir / "membrane_log.csv", 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=log_data[0].keys())
    writer.writeheader()
    writer.writerows(log_data)

# -------------------- UV SUMMARY --------------------
print(f"\n☀️ UV-degraded peptides: {uv_degraded_count} / {NUM_CYCLES}")
print(f"✅ Simulation complete. Log saved to {output_dir/'membrane_log.csv'}")

