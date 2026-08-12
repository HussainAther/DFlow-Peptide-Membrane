import csv, json
from pathlib import Path

def save_logs(results, outdir):
    logs_dir = outdir / "logs"
    logs_dir.mkdir(exist_ok=True)

    def save_csv(name, fields, rows):
        with open(logs_dir / name, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    save_csv("membrane_growth_log.csv", ["step", "thickness_nm"], results["membrane_log"])
    save_csv("peptide_log.csv", results["peptide_fields"], results["peptides"])
    save_csv("vesicle_state_log.csv", results["vesicle_fields"], results["vesicle"])

    with open(outdir / "final_state.json", "w") as f:
        json.dump(results["final_state"], f, indent=4)

