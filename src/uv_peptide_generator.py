import random
import argparse

class UVPeptideGenerator:
    def __init__(self, length=80, max_d_aromatics=2, preset="carbonaceous_chondrite"):
        self.length = length
        self.max_d_aromatics = max_d_aromatics
        self.preset = preset
        self.presets = {
            "carbonaceous_chondrite": {
                'A': 0.12, 'D': 0.08, 'E': 0.07, 'G': 0.15, 'V': 0.10,
                'L': 0.07, 'I': 0.05, 'P': 0.05, 'S': 0.08, 'T': 0.07,
                'F': 0.03, 'Y': 0.02, 'W': 0.02, 'others': 0.09
            },
            # Add more presets as needed
        }
        self.d_aromatics = {"D-F", "D-Y", "D-W"}

    def count_d_aromatics(self, peptide):
        return sum(1 for aa in peptide if aa in self.d_aromatics)

    def generate_single_peptide(self):
        aa_probs = self.presets[self.preset]
        aa_pool = list(aa_probs.keys())
        weights = [aa_probs[aa] for aa in aa_pool]

        while True:
            peptide = []
            for _ in range(self.length):
                aa = random.choices(aa_pool, weights=weights, k=1)[0]
                chirality = random.choice(["L", "D"])
                peptide.append(f"{chirality}-{aa}")
            if self.count_d_aromatics(peptide) <= self.max_d_aromatics:
                return peptide

    def generate_library(self, n=1000, output_file="uv_filtered_peptides.txt"):
        with open(output_file, "w") as f:
            for i in range(n):
                peptide = self.generate_single_peptide()
                f.write(f">peptide_{i+1}\n{' '.join(peptide)}\n")
        print(f"[✔] Generated {n} UV-filtered peptides → {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Generate UV-filtered LD peptide libraries.")
    parser.add_argument("--num", type=int, default=1000, help="Number of peptides to generate")
    parser.add_argument("--length", type=int, default=80, help="Length of each peptide")
    parser.add_argument("--preset", type=str, default="carbonaceous_chondrite", help="Amino acid composition preset")
    parser.add_argument("--max_d_aromatics", type=int, default=2, help="Maximum allowed D-aromatic residues")
    parser.add_argument("--output", type=str, default="uv_filtered_peptides.txt", help="Output file name")

    args = parser.parse_args()

    generator = UVPeptideGenerator(
        length=args.length,
        max_d_aromatics=args.max_d_aromatics,
        preset=args.preset
    )
    generator.generate_library(n=args.num, output_file=args.output)

if __name__ == "__main__":
    main()

