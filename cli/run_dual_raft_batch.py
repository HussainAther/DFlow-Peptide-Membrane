import argparse
import json
from pathlib import Path
from src.dual_raft_simulation import run_dual_raft_simulation
import datetime

def main(args):
    all_runs = []
    for i in range(args.runs):
        result = run_dual_raft_simulation(
            max_cycles=args.cycles,
            l_bias=args.lbias,
            seed=i + args.seed
        )
        all_runs.append(result)

    out_path = Path(args.outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    with open(out_path / f"dual_raft_results_{timestamp}.json", "w") as f:
        json.dump(all_runs, f, indent=2)

    print(f"✅ Saved {args.runs} simulations to {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=int, default=100)
    parser.add_argument("--cycles", type=int, default=1000)
    parser.add_argument("--lbias", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", type=str, default="experiments/raft_batch/")
    main(parser.parse_args())

