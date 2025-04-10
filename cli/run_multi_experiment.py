import subprocess
import argparse
import json
import yaml
from pathlib import Path
from datetime import datetime

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

def main():
    parser = argparse.ArgumentParser(description="Run multiple peptide simulation experiments.")
    parser.add_argument("--config", required=True, help="Path to parameter sweep file (YAML or JSON)")
    parser.add_argument("--base_outdir", default="experiments", help="Base output directory")
    args = parser.parse_args()

    base_dir = Path(args.base_outdir)
    base_dir.mkdir(exist_ok=True)

    sweep_config = load_sweep_config(args.config)

    for idx, exp in enumerate(sweep_config):
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

        print(f"\n🚀 Running experiment {idx+1}: {label}")
        print("📎 Command:", " ".join(cmd))
        subprocess.run(cmd)

if __name__ == "__main__":
    main()

