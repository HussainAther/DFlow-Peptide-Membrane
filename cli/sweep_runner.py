import itertools
import subprocess
from pathlib import Path
import argparse
import json

# --------- ARGUMENTS ---------
parser = argparse.ArgumentParser(description="Run a full parameter sweep grid")
parser.add_argument("--output_root", type=str, default="experiments", help="Parent folder for all sweeps")
parser.add_argument("--cycles", type=int, default=1000, help="Number of simulation cycles per run")
parser.add_argument("--mismatch_threshold", type=int, default=1, help="Hydrophobic mismatch threshold")
parser.add_argument("--flattened", action="store_true", help="Use flattened vesicle geometry")
args = parser.parse_args()

output_root = Path(args.output_root)
output_root.mkdir(parents=True, exist_ok=True)

# --------- PARAMETER GRID ---------
glycine_probs = [0.0, 0.01, 0.05, 0.1]
hydro_lengths = [4, 5, 6, 7]

grid = list(itertools.product(glycine_probs, hydro_lengths))

print(f"🔁 Running {len(grid)} total sweeps...")

# --------- RUN EACH EXPERIMENT ---------
for i, (g, h) in enumerate(grid, start=1):
    sweep_name = f"sweep_{i:03}_g{g}_h{h}"
    outdir = output_root / sweep_name
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python", "cli/run_experiment.py",
        str(g),
        str(h),
        str(args.mismatch_threshold),
        str(args.cycles),
        str(outdir)
    ]
    if args.flattened:
        cmd.append("--flattened")

    print(f"🚀 Running {sweep_name} ...")
    subprocess.run(cmd)

