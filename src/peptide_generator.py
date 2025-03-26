import torch
import random
import numpy as np
from dflow.model import PeptideFlow  # Assuming D-Flow model exists in repo

class PeptideGenerator:
    def __init__(self, num_peptides=100, sequence_length=15):
        self.num_peptides = num_peptides
        self.sequence_length = sequence_length
        self.model = PeptideFlow.load_pretrained("dflow_model.pt")  # Load trained D-Flow model
    
    def generate_peptides(self, hybrid=True):
        """Generate hybrid L-D peptides using D-Flow."""
        peptides = self.model.sample(self.num_peptides, seq_length=self.sequence_length, hybrid=hybrid)
        return peptides
    
    def mutate_peptide(self, peptide, mutation_rate=0.1):
        """Introduce mutations into a peptide sequence based on a mutation rate."""
        peptide_list = list(peptide)
        for i in range(len(peptide_list)):
            if np.random.rand() < mutation_rate:
                peptide_list[i] = np.random.choice(["L", "R", "0"])  # Random replacement
        return "".join(peptide_list)

    def generate_space_peptide(length=80, preset="carbonaceous_chondrite"):
        aa_probs = space_presets[preset]
        aa_pool = list(aa_probs.keys())
        weights = [aa_probs[aa] for aa in aa_pool]

        peptide = []
        for _ in range(length):
            aa = random.choices(aa_pool, weights=weights, k=1)[0]
            chirality = random.choice(["L", "D"])  # racemic
            peptide.append(f"{chirality}-{aa}")
        return " ".join(peptide)

    
    def recombine_peptides(self, peptide1, peptide2):
        """Perform crossover between two peptide sequences to simulate recombination."""
        if len(peptide1) != len(peptide2):
            raise ValueError("Peptides must be of the same length for recombination.")
        crossover_point = np.random.randint(1, len(peptide1) - 1)
        new_peptide = peptide1[:crossover_point] + peptide2[crossover_point:]
        return new_peptide

if __name__ == "__main__":
    generator = PeptideGenerator()
    peptides = generator.generate_peptides()
    mutated_peptide = generator.mutate_peptide(peptides[0])
    recombined_peptide = generator.recombine_peptides(peptides[0], peptides[1])
    
    print("Generated Peptides:", peptides[:5])
    print("Mutated Peptide:", mutated_peptide)
    print("Recombined Peptide:", recombined_peptide)

