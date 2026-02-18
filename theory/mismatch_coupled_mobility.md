# Mismatch-coupled mobility and raft diffusion in DFlow

In the reversible DFlow implementation, lateral mobility of peptide rafts is
treated as a hydrodynamic process on a two-dimensional viscous membrane. For a
raft of effective radius \(a\) embedded in a membrane with surface viscosity
\(\eta_{2D}\), the baseline mobility follows a Saffman–Delbrück–type scaling

\[
\mu_{0} = \frac{1}{4\pi \eta_{2D} a}.
\]

Hydrophobic mismatch is defined locally at each triad as

\[
\Delta h = \left| L_{p} - h_{\text{leaflet}} \right|,
\]

where \(L_{p}\) is the hydrophobic length of the peptide anchor and
\(h_{\text{leaflet}}\) is the local leaflet thickness inferred from the
amphiphile configuration. Motivated by Andersen–Koeppe and by mismatch-dependent
diffusion measurements (Ni, Fabian, Ramadurai and co-workers), we introduce a
quadratic mismatch penalty for the effective mobility,

\[
\mu_{\mathrm{eff}} = \frac{\mu_{0}}{1 + \gamma_{\text{mis}} \, \Delta h^{2}},
\]

where \(\gamma_{\text{mis}}\) is a dimensionless mismatch–mobility coupling
parameter (exposed in the code as `MISMATCH_MOBILITY_GAMMA`).

This rescaled mobility multiplies the attempt rates for the R2\(^+\)/R2\(^-\)
fission–fusion pair and the A1 amphiphile-swap event, so that rafts embedded in
strongly mismatched environments diffuse and remodel more slowly, while
well-matched domains retain higher mobilities. Importantly, the Metropolis
acceptance probability \(\exp(-\beta \Delta E)\) is left unchanged, so detailed
balance and microscopic reversibility are preserved even for large values of
\(\gamma_{\text{mis}}\).

For visualization, DFlow can render frames in a “mismatch” color mode, where
each triad is colored by its instantaneous mismatch magnitude \(\Delta h\) using
a continuous colormap (typically *viridis*) with optional per-frame
normalization. This provides a direct view of mismatch hotspots and helps
diagnose the symmetry of R2\(^+\) (fission at strained edges) and R2\(^-\)
(fusion relieving mismatch) under mismatch-scaled mobility.
