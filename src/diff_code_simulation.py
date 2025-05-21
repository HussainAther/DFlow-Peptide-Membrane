import random
import json
import matplotlib.pyplot as plt
from pathlib import Path

def generate_peptide(anchor_min=4, anchor_max=24, glycine_prob=0.01):
    """Generate a peptide with random chirality until it forms an anchor."""
    peptide = []
    max_length = 80
    for _ in range(max_length):
        r = random.random()
        if r < glycine_prob:
            aa = '0'
        else:
            aa = random.choice(['L', 'D'])
        peptide.append(aa)

        # Check for anchor termination
        if len(peptide) >= anchor_min and peptide[-5:] == [aa] * 5 and aa in ['L', 'D']:
            break
    return peptide

def anchor_length(peptide):
    """Return length of terminal homochiral anchor if valid, else None."""
    if len(peptide) < 5:
        return None
    anchor = peptide[-5:]
    if all(a == 'L' for a in anchor):
        return len(anchor)
    elif all(a == 'D' for a in anchor):
        return len(anchor)
    return None

def simulate_diff_raft(num_peptides=100, out_json='analysis/diff_code_output.json', out_img='analysis/final_raft_state.png'):
    L_raft, D_raft, cytoplasm = [], [], []
    log = []

    for i in range(num_peptides):
        peptide = generate_peptide()
        alen = anchor_length(peptide)
        if alen is None:
            cytoplasm.append(peptide)
            label = 'cytoplasm'
        elif peptide[-1] == 'L':
            L_raft.append((peptide, alen))
            label = 'L'
        elif peptide[-1] == 'D':
            D_raft.append((peptide, alen))
            label = 'D'

        log.append({
            "step": i + 1,
            "raft_L_count": len(L_raft),
            "raft_D_count": len(D_raft),
            "cytoplasm_count": len(cytoplasm),
            "added_label": label,
            "anchor_len": alen if alen else 0
        })

    # Save results
    Path(out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, 'w') as f:
        json.dump(log, f, indent=2)

    # Plot end-state raft histogram
    plt.figure(figsize=(8, 4))
    plt.hist([alen for _, alen in L_raft], bins=range(4, 25), alpha=0.6, label='L Raft')
    plt.hist([alen for _, alen in D_raft], bins=range(4, 25), alpha=0.6, label='D Raft')
    plt.xlabel('Anchor Length')
    plt.ylabel('Count')
    plt.title('Final Anchor Lengths in Rafts')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_img)
    plt.close()

    print(f"✅ Simulation complete. {len(L_raft)} L, {len(D_raft)} D, {len(cytoplasm)} cytoplasm peptides.")
    print(f"📊 Results saved to {out_json} and {out_img}")

if __name__ == "__main__":
    simulate_diff_raft()

