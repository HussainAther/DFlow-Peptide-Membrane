from .config import *
# from .peptide import ...
# from .membrane import ...
# from .vesicle import ...

def run_simulation(config):
    # Placeholder — insert simulation logic here
    membrane_log = []
    peptides = []
    vesicle = []
    final_state = {}

    return {
        "membrane_log": membrane_log,
        "peptides": peptides,
        "peptide_fields": ["index", "sequence", "length", "max_block", "inserted", "partial", "location", "hydrophobicity", "dominant_chirality"],
        "vesicle": vesicle,
        "vesicle_fields": ["step", "membrane_thickness", "peptides_membrane", "peptides_cytoplasm", "peptides_discarded", "raft_L", "raft_D"],
        "final_state": final_state
    }

