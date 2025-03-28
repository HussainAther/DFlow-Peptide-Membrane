import random

# Space-derived presets (amino acid distributions from meteorite samples)
SPACE_PRESETS = {
    "carbonaceous_chondrite": {
        "G": 0.14, "A": 0.13, "D": 0.12, "E": 0.08,
        "V": 0.07, "L": 0.06, "I": 0.05, "P": 0.04,
        "S": 0.08, "T": 0.05, "R": 0.04, "K": 0.03,
        "F": 0.03, "Y": 0.02, "C": 0.02, "M": 0.01,
        "H": 0.02, "Q": 0.02, "N": 0.02, "W": 0.01
    },
    "irradiated_meteorite": {
        # UV-degraded bias — more stable side chains survive
        "G": 0.18, "A": 0.16, "D":

