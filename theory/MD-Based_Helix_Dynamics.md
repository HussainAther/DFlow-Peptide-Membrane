# HYDROPHOBIC MISMATCH PAPERS REVIEW (MD-Based Helix Dynamics) – Varshith Madishetty

---

### **de Jesus & Allen (2013) — “Hydrophobic mismatch modulates transmembrane helix tilt, bending, and stretching in atomistic simulations”**
**Method/Model:** WALP peptides in bilayers; 1.4 μs atomistic MD with PMF decomposition (tilt, bending, stretching).  
**DFlow Hook:** Defines ΔE ∝ (ΔL)² energy relationship; provides atomistic baseline for mismatch energetics and Boltzmann weighting in DFlow.

---

### **Kim & Im (2010) — “Revisiting Hydrophobic Mismatch with Free Energy Simulation Studies of Transmembrane Helix Tilt and Rotation”**
**Method/Model:** WALPn peptides in DMPC/POPC; umbrella-sampling MD with tilt–rotation PMFs.  
**DFlow Hook:** Establishes tilt–rotation landscape; validates ΔE vs ΔL form used in DFlow’s β parameterization.

---

### **Morris et al. (2011) – “Exploring Hydrophobic Mismatch using MD Simulations of Gramicidin A in Lipid Bilayers”**
**Method/Model:** Gramicidin A channels in DLPC, DMPC, DOPC, POPC bilayers; all-atom MD examining mismatch-induced curvature.  
**DFlow Hook:** Shows local bilayer deformation and finite decay length; informs raft-size scaling and mismatch relaxation in DFlow.

---

### **Grau et al. (2019) – “The Role of Hydrophobic Mismatch on Transmembrane Helix Dimerization in Living Cells”**
**Method/Model:** Glycophorin A TM helices; cell assays (ToxRED/BiFC) + atomistic MD on mismatch-driven dimerization.  
**DFlow Hook:** Demonstrates cooperative mismatch adaptation; supports DFlow’s reversible R₁± raft association/dissociation rules.

---

### **Wallace et al. (2006) – “Effect of Hydrophobic Mismatch on Phase Behavior of Lipid Membranes”**
**Method/Model:** Binary lipid bilayer; Monte Carlo “mattress model” with J/κ ratio (surface tension vs compressibility).  
**DFlow Hook:** Connects continuum elasticity to DFlow’s ΔE–J/κ scans; framework for raft coarsening and stability modeling.

---

