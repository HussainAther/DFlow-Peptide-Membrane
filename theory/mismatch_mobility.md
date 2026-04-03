# Mismatch-Dependent Mobility and Raft Kinetics

In the **reversible DFlow model**, we allow the effective in-plane mobility of raft motions to depend on the local hydrophobic mismatch between embedded amphiphilic peptides and the surrounding bilayer.

At the continuum level, the lateral mobility of an inclusion of effective radius $a$ in a 2D viscous sheet with surface viscosity $\eta_{2D}$ is taken to follow a Stokes-like scaling:

$$
\mu_0 = \frac{1}{4\pi \eta_{2D} a},
$$

so that any kinetic rate $k$ controlled by lateral diffusion (e.g., rigid-body raft drift) scales as

$$
k \propto \mu_0.
$$

## Hydrophobic Mismatch Field

To connect this to the hydrophobic mismatch literature, we define a local mismatch field:

$$
\Delta h = \left| L_p - h_{\mathrm{leaflet}} \right|,
$$

where:
- $L_p$ is the hydrophobic length of a peptide (in practice $L_p \approx 0.15\,\text{nm}$ per residue in the current implementation),
- $h_{\mathrm{leaflet}}$ is the local leaflet thickness inferred from the carbon-chain composition of nearby lipids.

This mirrors the Andersen–Koeppe picture in which hydrophobic mismatch generates an elastic deformation of the bilayer and an associated energetic penalty that is approximately **quadratic in $\Delta h$** for small deformations.

## Mismatch-Dependent Effective Mobility

In the code, this is captured by a mismatch-dependent correction to the bare mobility:

$$
\mu_{\mathrm{eff}}(\Delta h) = \frac{\mu_0}{1 + \gamma \, \Delta h^2},
$$

where $\gamma = \texttt{MISMATCH MOBILITY GAMMA}$ is a dimensionless coupling parameter exposed in the top-level configuration files.

The quadratic dependence ensures symmetry between positive and negative mismatch ($+\Delta h$ vs. $-\Delta h$) and recovers the mismatch-free limit when $\gamma \to 0$:

$$
\gamma = 0 \quad \Rightarrow \quad \mu_{\mathrm{eff}} = \mu_0.
$$

## Operational Use in Kinetics

Operationally, we use $\mu_{\mathrm{eff}}$ as a multiplicative factor on diffusion-controlled event rates:

- **R4 rigid-body raft drift**:
  $$
  k_{\mathrm{R4}} \propto \mu_{\mathrm{eff}}(\Delta h),
  $$
  where $\Delta h$ may be computed per raft (e.g., average boundary mismatch) or from a global boundary-averaged field.

- **R2 fission/fusion and A1 swap moves** (optional, reversible branch):
  The same $\mu_{\mathrm{eff}}$ factor is used to rescale mismatch-sensitive association/dissociation attempt rates, so that strongly mismatched rafts drift and remodel **more slowly** than well-matched ones — consistent with higher effective friction.

Because the reversible R2 moves are Metropolis-accepted with respect to the total energy $E$ (including both bond and mismatch terms), this mobility prefactor changes **attempt frequencies** but **preserves detailed balance** as long as the proposal kernel remains symmetric in configuration space.

Validation:
- The sweep in `dflow_reversible.py` (energy vs. largest raft size) and
- the diffusion-histogram comparison in `analysis/raft_diffusion_interplay.py`

confirm that introducing $\gamma > 0$ leaves the **equilibrium size distribution** and **R4 step-length statistics** essentially unchanged, aside from a uniform slowdown consistent with reduced mobility.

## Relation to Hydrophobic Mismatch Literature

The functional form above is deliberately minimal but grounded in the standard hydrophobic mismatch framework:

- **Andersen & Koeppe** review how hydrophobic length mismatch between transmembrane segments and the host bilayer leads to elastic deformations and quadratic energetic penalties in $\Delta h$, motivating our use of a $\Delta h^2$ coupling in the mobility.
- **Ni, Ramadurai, Fabian and co-workers** quantify how peptide and protein diffusion in model membranes depends on bilayer thickness and mismatch, supporting the qualitative picture that strong mismatch increases effective friction and slows lateral motion.

In our implementation we **do not attempt to fit absolute mobilities**; instead, $\gamma$ plays the role of a **tunable phenomenological parameter** that can be swept in the reversible/diurnal drivers to explore regimes from:

- “mismatch-blind” ($\gamma=0$)
- to “strongly slowdown-coupled” ($\gamma \gg 1$).

These references provide the physical justification for treating hydrophobic mismatch as a slowly varying field that modulates the effective mobility of raft-scale motions, while the stochastic DFlow kinetics enforce microscopic reversibility through explicit Metropolis acceptance of R2 moves.