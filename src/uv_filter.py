# src/uv_filter.py

D_AROMATIC_RESIDUES = ['F', 'Y', 'W']

def count_d_aromatics(peptide):
    return sum(1 for i in range(len(peptide)-1)
               if peptide[i] == 'D' and peptide[i+1] in D_AROMATIC_RESIDUES)

def is_uv_degraded(peptide, max_d_aromatics=1, degradation_prob=0.5):
    from random import random
    num_d_aromatics = count_d_aromatics(peptide)
    if num_d_aromatics > max_d_aromatics:
        return random() < degradation_prob
    return False

