import argparse
import datetime
import json
from pathlib import Path
from src import config
from src.simulation import run_simulation
from src.logger import save_logs

# -----------------------------
# CLI Argument Parsing
# -----------------------------
parser = argparse.ArgumentParser(description="Run a peptide-membrane simulation experiment.")
parser.add_argument('--init', type=str, default='default', help='Initialization mode: default | space')
parser.add_argument('--num', type=int, default=10000, help='Number of peptides to simulate')
parser.add_argument('--margin', type=int, default=2, help='Hydrophobic block margin n')
parser.add_argument('--partial_prob', type=float, default=0.25, help='Probability to retain partial inserts')
parser.add_argument('--misfit_prob', type=float, default=0.5, help='Probability of cytoplasm retention on misfit')
parser.add_argument('--outdir', type=str, default=None, help='Custom output directory (default: timestamped)')

args = parser.parse_args()

# -----------------------------
# Set Up Output Directory
# -----------------------------
timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H%M%S")
outdir = Path(args.outdir) if args.outdir else Path(f"experiments/exp_{timestamp}")
outdir.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Override config values
# -----------------------------
config.NUM_PEPTIDES = args.num
config.HYDROPHOBIC_N_MARGIN = args.margin
config.PARTIAL_INSERTION_PROB = args.partial_prob
config.MISFITTED_INSERTION_PROB = args.misfit_prob
config.INIT_MODE = args.init

# -----------------------------
# Run Simulation
# -----------------------------
results = run_simulation(config)

# -----------------------------
# Save Outputs
# -----------------------------
save_logs(results, outdir)

# Save metadata
meta = {
    "timestamp": timestamp,
    "init_mode": args.init,
    "num_peptides": args.num,
    "hydrophobic_margin": args.margin,
    "partial_insertion_prob": args.partial_prob,
    "misfitted_insertion_prob": args.misfit_prob
}
with open(outdir / "meta.json", "w") as f:
    json.dump(meta, f, indent=4)

print(f"\n[✔] Experiment complete.")
print(f"Results saved to: {outdir.resolve()}")

