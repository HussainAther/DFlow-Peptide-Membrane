#!/usr/bin/env python3
"""
Generate canonical membrane configurations (C0–C5) and the configuration
transition graph for the Origin-of-Life manuscript.

Outputs:
    figures/configs/C0.png ... C5.png
    figures/configs/config_graph.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, Tuple, List

import matplotlib.pyplot as plt
from matplotlib import patches

from src.test_states.configurations import CONFIG_BUILDERS, ConfigState  # adjust import if needed

Coord = Tuple[int, int]

HEX_R = 1.0
SQRT3 = math.sqrt(3.0)


# --- hex helpers -------------------------------------------------------


def axial_to_xy(q: int, r: int, R: float = HEX_R) -> Tuple[float, float]:
    """
    Convert axial hex coords (q, r) to 2D positions for plotting.
    Pointy-top layout.
    """
    x = R * 1.5 * q
    y = R * (SQRT3 * (r + 0.5 * q))
    return x, y


def hex_vertices(x: float, y: float, R: float = HEX_R) -> List[Tuple[float, float]]:
    return [
        (
            x + R * math.cos(math.radians(60 * k + 30)),
            y + R * math.sin(math.radians(60 * k + 30)),
        )
        for k in range(6)
    ]


# --- rendering ---------------------------------------------------------


def render_config(cfg: ConfigState, out_path: Path) -> None:
    """
    Render a single configuration as a simple hex lattice figure.

    This is intentionally minimalist: background hexes, with
    peptides drawn as colored markers (L vs D in different colors).
    """
    triads = cfg.triads
    amph = cfg.amph
    peptides = cfg.peptides

    fig, ax = plt.subplots(figsize=(4, 4))

    # Draw background hexes
    for (q, r) in triads:
        x, y = axial_to_xy(q, r)
        face = "#F3F6FC"
        edge = "#CBD5E8"

        # optional: darker edge if peptide present
        if amph[(q, r)]["pep"]:
            face = "#E2EDFF"
            edge = "#8898C5"

        poly = patches.Polygon(
            hex_vertices(x, y),
            closed=True,
            linewidth=0.6,
            edgecolor=edge,
            facecolor=face,
            zorder=1,
        )
        ax.add_patch(poly)

    # Draw peptides as circles on top
    for pid, pep in peptides.items():
        q, r = pep["cent"]
        x, y = axial_to_xy(q, r)

        chir = pep.get("chir", "L")
        if chir == "L":
            color = "#5B2C6F"  # purple
        elif chir == "D":
            color = "#B22222"  # red
        else:
            color = "#333333"

        ax.scatter(x, y, s=60, color=color, edgecolor="black", zorder=3)

    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title(f"{cfg.label}: {cfg.description}", fontsize=9)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)


# --- config graph ------------------------------------------------------


def build_config_graph() -> Dict:
    """
    Define a simple configuration graph for the manuscript.

    Nodes are C0–C5.
    Edges are labeled by the dominant event class that links them.
    This is a conceptual graph, not a full Markov chain.
    """
    nodes = []
    edges = []

    # node list
    for label, builder in CONFIG_BUILDERS.items():
        cfg: ConfigState = builder()
        nodes.append(
            {
                "id": cfg.label,
                "description": cfg.description,
            }
        )

    # conceptual transitions (you can refine this)
    # C0 --(P1/P2)--> C1 --> C2 --> C3 --> C4 --> C5, etc.
    def add_edge(src: str, dst: str, event: str):
        edges.append({"source": src, "target": dst, "event": event})

    add_edge("C0", "C1", "P1/P2")
    add_edge("C1", "C2", "P2")
    add_edge("C2", "C3", "R1")
    add_edge("C3", "C4", "R2")
    add_edge("C4", "C5", "R2/R1")
    # reversible arrows (conceptual)
    add_edge("C2", "C1", "P2⁻")
    add_edge("C3", "C2", "R1⁻")
    add_edge("C4", "C3", "R2⁻")
    add_edge("C5", "C4", "R1⁻/R2⁻")

    return {"nodes": nodes, "edges": edges}


# --- main --------------------------------------------------------------


def main(out_root: str = "figures/configs") -> None:
    out_root_path = Path(out_root)
    out_root_path.mkdir(parents=True, exist_ok=True)

    # 1. Render each configuration
    for label, builder in CONFIG_BUILDERS.items():
        cfg = builder()
        out_path = out_root_path / f"{label}.png"
        print(f"[gen_configurations] Rendering {label} -> {out_path}")
        render_config(cfg, out_path)

    # 2. Save configuration graph JSON
    graph = build_config_graph()
    graph_path = out_root_path / "config_graph.json"
    with graph_path.open("w", encoding="utf-8") as f:
        json.dump(graph, f, indent=2)
    print(f"[gen_configurations] Wrote config graph to {graph_path}")


if __name__ == "__main__":
    main()
