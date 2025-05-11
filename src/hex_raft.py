import random
import numpy as np
from collections import deque

class Peptide:
    def __init__(self, chirality, anchor_len):
        self.chirality = chirality
        self.anchor_len = anchor_len

class HexRaft:
    def __init__(self, chirality):
        self.chirality = chirality
        self.peptides = []

    def add_peptide(self, peptide):
        if peptide.chirality == self.chirality:
            self.peptides.append(peptide)

    def sort_and_place(self):
        self.peptides.sort(key=lambda p: p.anchor_len, reverse=True)

        # Place in spiral hex grid: (q, r) axial coords
        spiral_coords = self.generate_hex_spiral(len(self.peptides))
        placed = {coord: pep for coord, pep in zip(spiral_coords, self.peptides)}
        return placed

    def get_perimeter_avg(self):
        placed = self.sort_and_place()
        max_ring = max(max(abs(q), abs(r), abs(-q-r)) for q, r in placed)

        perimeter_peps = [pep.anchor_len for (q, r), pep in placed.items()
                          if max(abs(q), abs(r), abs(-q-r)) == max_ring]

        return np.mean(perimeter_peps) if perimeter_peps else 0

    def generate_hex_spiral(self, n):
        coords = []
        q = r = 0
        directions = [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]

        ring = 0
        while len(coords) < n:
            if ring == 0:
                coords.append((0, 0))
                ring += 1
                continue

            q, r = -ring, 0
            for d in directions:
                for _ in range(ring):
                    if len(coords) >= n:
                        break
                    coords.append((q, r))
                    q += d[0]
                    r += d[1]
            ring += 1
        return coords

