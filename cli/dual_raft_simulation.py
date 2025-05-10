import random

class Peptide:
    def __init__(self, chirality, anchor_len):
        self.chirality = chirality  # 'L' or 'D'
        self.anchor_len = anchor_len

class MembraneRaft:
    def __init__(self, start_thickness=12, min_thick=10, max_thick=23):
        self.thickness = start_thickness
        self.min = min_thick
        self.max = max_thick
        self.history = [start_thickness]
        self.crowded = 0  # peptides not inserted

    def try_insert(self, peptide):
        delta = peptide.anchor_len - self.thickness
        inserted = False

        if abs(delta) <= 1:
            self.thickness += (1 if delta > 0 else -1 if delta < 0 else 0)
            self.thickness = max(self.min, min(self.max, self.thickness))
            inserted = True
        else:
            self.crowded += 1  # peptide rejected

        self.history.append(self.thickness)
        return inserted

    def is_done(self):
        return self.thickness == self.min or self.thickness == self.max

def generate_peptide(min_len=6, max_len=20, l_bias=0.5):
    chirality = random.choices(['L', 'D'], weights=[l_bias, 1 - l_bias])[0]
    anchor = random.randint(min_len, max_len)
    return Peptide(chirality, anchor)

def run_dual_raft_simulation(
    max_cycles=1000, l_bias=0.5, seed=None, min_len=6, max_len=20
):
    if seed is not None:
        random.seed(seed)

    L_raft = MembraneRaft()
    D_raft = MembraneRaft()

    for _ in range(max_cycles):
        p = generate_peptide(l_bias=l_bias, min_len=min_len, max_len=max_len)
        if p.chirality == 'L':
            L_raft.try_insert(p)
        else:
            D_raft.try_insert(p)

        if L_raft.is_done() and D_raft.is_done():
            break

    return {
        "L_history": L_raft.history,
        "D_history": D_raft.history,
        "L_crowded": L_raft.crowded,
        "D_crowded": D_raft.crowded,
        "cycles": len(L_raft.history)
    }

