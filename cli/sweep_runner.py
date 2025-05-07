import itertools
import subprocess
from pathlib import Path
import argparse
import json
import pandas as pd

# --------- ARGUMENT PARSING ---------
parser = argparse.ArgumentParser(description="Run a full parameter sweep grid and save summaries")
parser.add_argument("--output_root", type=str, default="experiments", help="Parent folder for all sweeps")
parser.add_argument("--cycles", type=int, default=1000, help="Number of simulation cycles per run")
parser.add_argument("--mismatch_threshold", type=int, default=1, help="Hydrophobic mismatch threshold")
parser.add_argument("--flattened", action="store_true", help="Use flattened vesicle geometry")
args = parser.parse_args()

output_root = Path(args.output_root)
output_root.mkdir(parents=True, exist_ok=True)

# --------- DEFINE SWEEP GRID ---------
glycine_probs = [0.0, 0.01, 0.05, 0.1]
hydro_lengths = [4, 5, 6, 7]
grid = list(itertools.product(glycine_probs, hydro_lengths))

summary_data = []

print(f"🔁 Running {len(grid)} total sweeps...\n")

# --------- EXECUTE EACH SWEEP ---------
for i, (g, h) in enumerate(grid, start=1):
    sweep_name = f"sweep_{i:03}_g{g}_h{h}"
    outdir = output_root / sweep_name
    outdir.mkdir(parents=True, exist_ok=True)

    # Build run command
    cmd = [
        "python", "cli/run_experiment.py",
        str(g),                    # glycine_prob
        str(h),                    # hydro_block_len
        str(args.mismatch_threshold),
        str(args.cycles),
        str(outdir)
    ]
    if args.flattened:
        cmd.append("--flattened")

    print(f"🚀 Running {sweep_name} ...")
    subprocess.run(cmd)

    # Try to read final membrane state
    try:
        log = pd.read_csv(outdir / "membrane_log.csv")
        final = log.iloc[-1]
        summary_data.append({
            "sweep": sweep_name,
            "glycine_prob": g,
            "hydro_block_len": h,
            "final_thickness": final.get("MembraneThickness", None),
            "final_LD_ratio": final.get("L_to_D_Ratio", None),
            "final_volume": final.get("VesicleVolume_nm3", None),
            "final_density": final.get("PeptideDensity_per_nm2", None),
            "cycles": len(log)
        })
    except Exception as e:
        print(f"⚠️ Could not parse {sweep_name}: {e}")
        continue

# --------- WRITE SUMMARY JSON ---------
summary_file = output_root / "sweep_summary.json"
with open(summary_file, "w") as f:
    json.dump(summary_data, f, indent=4)

print(f"\n✅ Sweep summary saved to {summary_file}")

