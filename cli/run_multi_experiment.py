import subprocess
import argparse
import json
import yaml
import time
from pathlib import Path
from datetime import datetime
from tqdm import tqdm  # pip install tqdm

def load_sweep_config(file_path):
    path = Path(file_path)
    if path.suffix in [".yaml", ".yml"]:
        with open(path, "r") as f:
            return yaml.safe_load(f)
    elif path.suffix == ".json":
        with open(path, "r") as f:
            return json.load(f)
    else:
        raise ValueError("Unsupported file format. Use .json or .yaml")

def run_analysis(outdir):
    try:
        subprocess.run([
            "python", "analysis/chirality_tracker.py",
            "--log_path", str(Path(outdir) / "logs" / "peptide_log.csv"),
            "--outdir", str(Path(outdir) / "analysis")
        ])
    except Exception as e:
        print(f"[⚠] Analysis failed for {outdir}: {e}")

def extract_summary(outdir):
    final_state_file = Path(outdir) / "final_state.json"
    if final_state_file.exists():
        with open(final_state_file, "r") as f:
            data = json.load(f)
        return {
            "membrane_thickness": data.get("membrane_thickness", "NA"),
            "days_elapsed": data.get("days_elapsed", "NA"),
            "raft_L": data.get("raft_L_count", "NA"),
            "raft_D": data.get("raft_D_count", "NA")
        }
    return {"membrane_thickness": "NA", "days_elapsed": "NA", "raft_L": "NA", "raft_D": "NA"}

def main():
    parser = argparse.ArgumentParser(description="Run multiple peptide simulation experiments.")
    parser.add_argument("--config", required=True, help="Path to parameter sweep file (YAML or JSON)")
    parser.add_argument("--base_outdir", default="experiments", help="Base output directory")
    args = parser.parse_args()

    base_dir = Path(args.base_outdir)
    base_dir.mkdir(exist_ok=True)

    sweep_config = load_sweep_config(args.config)
    sweep_summary = []

    for idx, exp in tqdm(enumerate(sweep_config), total=len(sweep_config), desc="Running sweep"):

        label = exp.get("label", f"exp_{idx+1}")
        timestamp = datetime.now().strftime("%Y_%m_%d_%H%M%S")
        outdir = base_dir / f"sweep_{str(idx+1).zfill(3)}_{label}"

        cmd = [
            "python", "cli/run_experiment.py",
            "--init", exp.get("init", "default"),
            "--num", str(exp.get("num", 1000)),
            "--margin", str(exp.get("margin", 2)),
            "--partial_prob", str(exp.get("partial_prob", 0.25)),
            "--misfit_prob", str(exp.get("misfit_prob", 0.5)),
            "--outdir", str(outdir)
        ]

        print(f"\n🚀 Experiment {idx+1}: {label}")
        subprocess.run(cmd)

        # Auto-run analysis
        run_analysis(outdir)

        # Extract and log summary
        summary = extract_summary(outdir)
        summary["label"] = label
        summary["folder"] = str(outdir)
        sweep_summary.append(summary)

    # Save sweep summary
    summary_file = base_dir / "sweep_summary.json"
    with open(summary_file, "w") as f:
        json.dump(sweep_summary, f, indent=4)

    print("\n✅ Sweep completed. Summary written to:", summary_file)

if __name__ == "__main__":
    main()

