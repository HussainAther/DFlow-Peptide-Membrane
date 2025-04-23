import os
import subprocess
from itertools import product
from pathlib import Path

# Define parameter ranges
glycine_fractions = [0.01, 0.05]
anchor_lengths = [5, 6]
mismatch_thresholds = [5, 7]

# Number of peptides per run
num_cycles = 1000

run_id = 1
base_output = Path("experiments/batch_runs")
base_output.mkdir(parents=True, exist_ok=True)

for g, m, threshold in product(glycine_fractions, anchor_lengths, mismatch_thresholds):
    run_name = f"run_{run_id:03d}_g{g}_m{m}_t{threshold}"
    output_dir = base_output / run_name

    cmd = [
        "python", "membrane_simulation.py",
        str(g), str(m), str(threshold), str(num_cycles), str(output_dir)
    ]

    print(f"🚀 Starting {run_name} ...")
    subprocess.run(cmd)

    run_id += 1

print("\n✅ All batch simulations completed!")

