import matplotlib.pyplot as plt
from collections import Counter
from src.peptide_generator import PeptideGenerator
from src.membrane_model import MembraneModel

# Setup
generator = PeptideGenerator(num_peptides=100, sequence_length=10)
membrane = MembraneModel()

# Generate peptides
peptides = generator.generate_peptides()

# Track results
stable = []
rejected = []
rejection_reasons = []

# Redefine check with reason logging
def check_stability_with_reason(peptide):
    seq = peptide.get("sequence", "")
    length = len(seq)

    # Length constraint
    if length < 10:
        return False, "Too short"
    if length > 20:
        return False, "Too long"

    # Charge (simulate realistic range)
    charge = peptide.get("charge", 0)
    if abs(charge) > 2:
        return False, "Too polar"

    # Chirality must be present
    chirality = peptide.get("chirality", "")
    if chirality not in ["L", "D"]:
        return False, "Bad chirality"

    # Hydrophobicity check (simple)
    hydrophobic_residues = set("AILMFWV")
    hydrophobic_count = sum(1 for aa in seq if aa in hydrophobic_residues)
    if hydrophobic_count / length < 0.3:
        return False, "Too hydrophilic"

    return True, "Stable"


# Evaluate
for p in peptides:
    ok, reason = check_stability_with_reason(p)
    if ok:
        stable.append(p)
    else:
        rejected.append(p)
        rejection_reasons.append(reason)

# Print summary
print(f"Total peptides: {len(peptides)}")
print(f"Stable: {len(stable)}")
print(f"Rejected: {len(rejected)}")

# Plot rejection reasons
counts = Counter(rejection_reasons)
plt.bar(counts.keys(), counts.values(), color="coral")
plt.title("Rejection Reasons for Peptides")
plt.ylabel("Count")
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig("data/stability_rejection_summary.png")
plt.close()

