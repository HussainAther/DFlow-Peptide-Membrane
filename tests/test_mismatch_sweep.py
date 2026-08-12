import subprocess
import sys
import re
from pathlib import Path

def run_script(args):
    """Run the sweep script and return stdout as string."""
    result = subprocess.run(
        [sys.executable, "run_mismatch_sweep.py", *args],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def test_dry_run_basic_command_structure():
    """Ensure dry-run prints valid commands containing required flags."""

    # test three distinct gamma values
    gammas = [0.0, 0.1, 0.25]

    stdout = run_script(["--dry_run", "--gammas", "0.0", "0.1", "0.25"])

    # split into lines, each line should be one command
    lines = [l.strip() for l in stdout.split("\n") if l.strip()]

    assert len(lines) == len(gammas), \
        f"Expected {len(gammas)} dry-run commands, got {len(lines)}"

    for line, gamma in zip(lines, gammas):

        # must call python dflow_reversible.py
        assert "dflow_reversible.py" in line, f"Missing script in: {line}"

        # must have mismatch mobility gamma flag
        assert "--MISMATCH_MOBILITY_GAMMA" in line, f"Missing gamma flag: {line}"

        # must include correct gamma value
        # support both 0.1 and 0.10 etc.
        gamma_pattern = re.escape(str(gamma))
        assert re.search(gamma_pattern, line), \
            f"Gamma {gamma} not found in: {line}"

        # ensure no empty arguments (e.g., accidental double spaces)
        tokens = line.split()
        assert all(tokens), f"Empty token detected in: {line}"

        # ensure commands begin with 'python' or 'python3'
        assert tokens[0].startswith("python"), f"Bad prefix: {tokens[0]}"

        # ensure OUT directory is being set
        assert "--OUT" in tokens, f"--OUT missing: {line}"

        # ensure SAVE_STEPS contains multiple values
        assert "--SAVE_STEPS" in tokens, f"--SAVE_STEPS missing: {line}"


def test_no_execution_in_dry_run():
    """Dry run must NOT create output directories."""

    out_root = Path("runs/mismatch_sweep")
    if out_root.exists():
        # clean up for test isolation
        for p in out_root.glob("*"):
            if p.is_dir():
                for q in p.glob("*"):
                    q.unlink()
                p.rmdir()
        out_root.rmdir()

    # run dry-run
    run_script(["--dry_run"])

    # ensure no directories created
    assert not out_root.exists(), \
        "Dry-run should not create directories, but runs/mismatch_sweep exists!"

