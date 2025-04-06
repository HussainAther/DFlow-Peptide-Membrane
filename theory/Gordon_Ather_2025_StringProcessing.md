# A String-Processing Coarse-Grained Approach to the Origin of Peptides During the Origin of Life

**Richard Gordon**  
Gulf Specimen Marine Laboratory & Aquarium  
Panacea, FL 32346, USA  
DickGordonCan@protonmail.com  

**Syed Hussain Ather**

---

We have previously presented a qualitative hypothesis for the origin of membrane peptides during the origin of life based on a positive feedback mechanism involving hydrophobic mismatching [1,2]. Here we attempt a plausible quantitation of this hypothesis. For other models, see [6].

---

## 1. Assumptions and Framework

1. **Archaea First Hypothesis**:  
   We assume that early vesicle membranes were highly permeable, especially compared to modern Bacteria [5].  
   → Amino acids and other monomers could freely diffuse into vesicles in a prebiotic lake setting [3,4].

2. **Flat Membranes**:  
   Vesicles are assumed to remain flattened (not spherical) during simulation [7,8].

3. **Peptide Rafts Only**:  
   No lipid rafts are assumed to form abiotically [9]. Only peptide rafts (L-rafts and D-rafts) are modeled.

4. **Daily Polymerization**:  
   Peptides form one amino acid per diurnal cycle. When a hydrophobic sequence of 4+ like-chiral residues appears, it can form an \(\alpha\)-helix and insert into the membrane.

---

## 2. Chirality and Peptide Dynamics

5. **LD0 Model**:  
   Three monomer types are used:
   - `L`: Levorotary amino acids  
   - `D`: Dextrorotary amino acids  
   - `0`: Achiral (e.g., glycine)

   A sequence of 4+ consecutive `L` or `D` forms a hydrophobic \(\alpha\)-helix, which can insert into the membrane **if** its length matches current membrane thickness + 1 residue (mismatch rule).  
   Peptides without such motifs exit or flip/flop to the inner membrane surface.

6. **Termination**:  
   Once transferred, the peptide terminates and a new one begins polymerizing.

---

## 3. Membrane Feedback and Raft Growth

7. **Raft Tracking**:  
   We track two peptide rafts:
   - **L-Rafts**: Peptides with L-block helices  
   - **D-Rafts**: Peptides with D-block helices  
   Each entry logs the dominant chirality and helix length.

8. **Membrane Growth Trigger**:  
   When the number of inserted peptides with mismatched length reaches the **perimeter threshold**, membrane thickness is increased to match. Amphiphile supply is assumed not limiting.

9. **Simplification**:  
   No tilted \(\alpha\)-helices are modeled — insertions are assumed perpendicular to membrane plane.

10. **Termination Condition**:  
   The simulation ends when membrane thickness reaches 4 nm. Output includes total simulation time (days), raft stats, and chirality distributions.

---

## References

[1] Gordon, R. & Gordon, N.K. (2023). *How to make a transmembrane domain at the origin of life.* In: **Conflicting Models for the Origin of Life**, Eds. Smoukov, Seckbach & Gordon, Wiley-Scrivener: 131–174.

[2] Gordon, R. & Gordon, N.K. (2024). *How to make a transmembrane domain at the origin of life.* In: **Origin of Life via Archaea**, Eds. Gordon & Seckbach, Wiley-Scrivener: 229–284.

[3] Gordon, R. et al. (2023). *The fish ladder toy model...* In: **Conflicting Models for the Origin of Life**, Wiley-Scrivener: 425–458.

[4] Gordon, R. et al. (2024). *The fish ladder toy model...* In: **Origin of Life via Archaea**, Wiley-Scrivener: 185–228.

[5] Lapinska, U. et al. (2023). *Archaeal membranes more permeable than bacterial ones*. PLoS Biol 21(4): e3002048.

[6] Mulkidjanian, A.Y. et al. (2009). *Co-evolution of membranes and proteins*. Trends Biochem Sci 34(4): 206–215.

[7] Osipenko, D.S. et al. (2016). *Redistribution in curved membranes*. Biol Membr 33(3): 176–186.

[8] Osipenko, D.S. et al. (2016). *Redistribution in lipid domains*. Biochem Moscow Supplement A 10(4): 259–268.

[9] Wang, C. et al. (2015). *Push and pull in lipid raft formation*. J Am Chem Soc 137(2): 664–666.

---

