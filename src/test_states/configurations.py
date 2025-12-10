from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, Set, Any, List

Coord = Tuple[int, int]


@dataclass
class ConfigState:
    """
    Canonical configuration on the hexagonal triad lattice.

    This is a lightweight container that mirrors the internal DFlow
    membrane data model closely enough for visualization and testing.

    Fields:
        label: config ID, e.g. "C0", "C1", ...
        description: human-readable text for the manuscript/figures.
        triads: set of axial coordinates (q, r) in the rhombic patch.
        amph: mapping (q, r) -> amphiphile dict
              { "carbon": int, "mono_di": 1|2, "pep": bool }
        peptides: mapping pid -> peptide dict
              {
                "cent": (q, r),
                "orient": int,
                "inside": bool,
                "chir": "L"|"D",
                "length": int,
                "Lp_nm": float,
                "helical": bool,
              }
        rafts: mapping raft_label -> set of peptide IDs
    """
    label: str
    description: str
    triads: Set[Coord]
    amph: Dict[Coord, Dict[str, Any]]
    peptides: Dict[int, Dict[str, Any]]
    rafts: Dict[str, Set[int]]


# ---- helpers ----------------------------------------------------------


def _make_rhombic_triads(N0: int) -> Set[Coord]:
    triads: Set[Coord] = set()
    for q in range(-N0, N0 + 1):
        for r in range(-N0, N0 + 1):
            if abs(q) <= N0 and abs(r) <= N0 and abs(q + r) <= N0:
                triads.add((q, r))
    return triads


def _checker_amph(triads: Set[Coord]) -> Dict[Coord, Dict[str, Any]]:
    amph: Dict[Coord, Dict[str, Any]] = {}
    for (q, r) in triads:
        mono_di = 1 if (q + r) % 2 == 0 else 2
        amph[(q, r)] = {
            "carbon": 18,
            "mono_di": mono_di,
            "pep": False,
        }
    return amph


def _add_peptides(
    amph: Dict[Coord, Dict[str, Any]],
    coords: List[Coord],
    chir: str,
    start_pid: int = 0,
    length: int = 10,
    Lp_nm: float = 3.0,
) -> Tuple[Dict[int, Dict[str, Any]], int]:
    peptides: Dict[int, Dict[str, Any]] = {}
    pid = start_pid
    for coord in coords:
        peptides[pid] = {
            "cent": coord,
            "orient": 0,
            "inside": True,
            "chir": chir,
            "length": length,
            "Lp_nm": Lp_nm,
            "helical": True,
        }
        amph[coord]["pep"] = True
        pid += 1
    return peptides, pid


# ---- canonical configurations -----------------------------------------


def make_C0_bare_membrane(N0: int = 3) -> ConfigState:
    """
    C0: bare amphiphile membrane with no peptides.

    Used as the "origin" configuration in the atlas.
    """
    triads = _make_rhombic_triads(N0)
    amph = _checker_amph(triads)
    peptides: Dict[int, Dict[str, Any]] = {}
    rafts: Dict[str, Set[int]] = {}

    return ConfigState(
        label="C0",
        description="Bare amphiphile membrane (no peptides).",
        triads=triads,
        amph=amph,
        peptides=peptides,
        rafts=rafts,
    )


def make_C1_single_peptide(N0: int = 3) -> ConfigState:
    """
    C1: bare membrane with a single isolated peptide.

    Used to illustrate the first non-trivial peptide state.
    """
    triads = _make_rhombic_triads(N0)
    amph = _checker_amph(triads)

    center = (0, 0)
    if center not in triads:
        # fallback: pick any triad
        center = next(iter(triads))

    peptides, _ = _add_peptides(amph, [center], chir="L", start_pid=0, length=8, Lp_nm=2.5)
    rafts: Dict[str, Set[int]] = {"L": {0}}

    return ConfigState(
        label="C1",
        description="Single peptide on a rhombic membrane patch.",
        triads=triads,
        amph=amph,
        peptides=peptides,
        rafts=rafts,
    )


def make_C2_mini_raft(N0: int = 3) -> ConfigState:
    """
    C2: Minimal 2-peptide raft.

    Used as the smallest raft configuration for R1 and P2 tests.
    """
    triads = _make_rhombic_triads(N0)
    amph = _checker_amph(triads)

    coords = [(0, 0), (1, 0)]
    coords = [c for c in coords if c in triads]
    peptides, next_pid = _add_peptides(amph, coords, chir="L", start_pid=0)
    rafts: Dict[str, Set[int]] = {"L": set(peptides.keys())}

    return ConfigState(
        label="C2",
        description="Minimal 2-peptide raft on the membrane.",
        triads=triads,
        amph=amph,
        peptides=peptides,
        rafts=rafts,
    )


def make_C3_compact_L_raft(N0: int = 4) -> ConfigState:
    """
    C3: Compact L-raft (2x2 block) on the membrane.

    Good for R1 association/dissociation and R2 fusion examples.
    """
    triads = _make_rhombic_triads(N0)
    amph = _checker_amph(triads)

    raft_coords = [(-1, 2), (0, 2), (-1, 3), (0, 3)]
    raft_coords = [c for c in raft_coords if c in triads]

    peptides, _ = _add_peptides(amph, raft_coords, chir="L", start_pid=0)
    rafts: Dict[str, Set[int]] = {"L": set(peptides.keys())}

    return ConfigState(
        label="C3",
        description="Compact 2x2 L-raft embedded in the membrane.",
        triads=triads,
        amph=amph,
        peptides=peptides,
        rafts=rafts,
    )


def make_C4_elongated_D_raft(N0: int = 4) -> ConfigState:
    """
    C4: Elongated D-raft with a neck.

    Used for R1 bonds along the bar, R2 fission/fusion, and R4 drift.
    """
    triads = _make_rhombic_triads(N0)
    amph = _checker_amph(triads)

    D_raft_sites = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1)]
    D_raft_sites = [c for c in D_raft_sites if c in triads]

    peptides, _ = _add_peptides(amph, D_raft_sites, chir="D", start_pid=0)
    rafts: Dict[str, Set[int]] = {"D": set(peptides.keys())}

    return ConfigState(
        label="C4",
        description="Elongated D-raft with a neck (for R2 and R4 examples).",
        triads=triads,
        amph=amph,
        peptides=peptides,
        rafts=rafts,
    )


def make_C5_two_rafts(N0: int = 4) -> ConfigState:
    """
    C5: Two separate rafts (L and D) on the same membrane.

    Good for multi-raft configurations and R1/R2 between rafts.
    """
    triads = _make_rhombic_triads(N0)
    amph = _checker_amph(triads)

    D_coords = [(0, 0), (1, 0), (0, 1)]
    L_coords = [(-2, 2), (-1, 2), (-2, 3), (-1, 3)]

    D_coords = [c for c in D_coords if c in triads]
    L_coords = [c for c in L_coords if c in triads]

    peptides: Dict[int, Dict[str, Any]] = {}
    rafts: Dict[str, Set[int]] = {"D": set(), "L": set()}

    pid = 0
    D_peps, pid = _add_peptides(amph, D_coords, chir="D", start_pid=pid)
    L_peps, pid = _add_peptides(amph, L_coords, chir="L", start_pid=pid)

    peptides.update(D_peps)
    peptides.update(L_peps)
    rafts["D"].update(D_peps.keys())
    rafts["L"].update(L_peps.keys())

    return ConfigState(
        label="C5",
        description="Two separate rafts (D and L) on the same membrane.",
        triads=triads,
        amph=amph,
        peptides=peptides,
        rafts=rafts,
    )


# Registry of configuration constructors
CONFIG_BUILDERS = {
    "C0": make_C0_bare_membrane,
    "C1": make_C1_single_peptide,
    "C2": make_C2_mini_raft,
    "C3": make_C3_compact_L_raft,
    "C4": make_C4_elongated_D_raft,
    "C5": make_C5_two_rafts,
    # TODO: extend to C6–C11 as needed
}
