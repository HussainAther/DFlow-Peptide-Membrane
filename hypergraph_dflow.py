#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hypergraph_dflow.py

DFlow hypergraph / incidence-graph utilities with R4_drift,
with overlay labels moved OUTSIDE the rafts (Option C).

Changes in this version:
- Anchor labels (D1, D2, D3, L1) placed outside raft region
  with thin leader lines pointing to their anchor hexes.
- Labels are small, no filled background boxes.
- R4_drift arrow label kept small and out of the way.
"""

import math
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as patches

try:
    import networkx as nx
except ImportError as e:
    raise SystemExit(
        "This script requires networkx. Install it with:\n"
        "    pip install networkx\n"
        f"Original error: {e}"
    )

# ---------------------------------------------------------------------------
# 1. Hypergraph / incidence graph specification
# ---------------------------------------------------------------------------

ENTITY_NODES: List[str] = [
    "D1_anchor",
    "D2_anchor",
    "D3_anchor",
    "L1_anchor",
]

EVENT_TO_ENTITIES: Dict[str, List[str]] = {
    "A1_forward": ["D1_anchor", "D2_anchor", "D3_anchor"],
    "A1_reverse": ["D1_anchor", "D2_anchor", "D3_anchor"],
    "P2_insert": ["D1_anchor", "D2_anchor", "D3_anchor"],
    "P2_desorb": ["D1_anchor", "D2_anchor", "D3_anchor"],
    "R1_assoc": ["D2_anchor", "L1_anchor"],
    "R1_dissoc": ["D2_anchor", "L1_anchor"],
    "R2_fusion": ["D1_anchor", "L1_anchor"],
    "P4_flip_LtoD": ["L1_anchor"],
    "P4_flip_DtoL": ["D1_anchor"],
    "R4_drift": ["D1_anchor", "D2_anchor", "D3_anchor", "L1_anchor"],
}


def build_incidence_graph() -> "nx.Graph":
    """Build a bipartite incidence graph (entities ↔ events)."""
    G = nx.Graph()

    for ent in ENTITY_NODES:
        G.add_node(ent, kind="entity", bipartite=0)

    for ev, ents in EVENT_TO_ENTITIES.items():
        G.add_node(ev, kind="event", bipartite=1)
        for ent in ents:
            if ent not in G:
                raise ValueError(f"Unknown entity '{ent}' referenced by event '{ev}'")
            G.add_edge(ev, ent)

    return G


def draw_incidence_graph(G: "nx.Graph", outpath: str = "hypergraph_incidence.png") -> None:
    """Draw the incidence graph with events on top and entities at the bottom."""
    entities = [n for n, d in G.nodes(data=True) if d.get("kind") == "entity"]
    events = [n for n, d in G.nodes(data=True) if d.get("kind") == "event"]

    entities_sorted = sorted(entities)
    events_sorted = sorted(events)

    pos = {}
    n_ent = len(entities_sorted)
    n_evt = len(events_sorted)

    for i, ent in enumerate(entities_sorted):
        x = i if n_ent <= 1 else i / (n_ent - 1)
        pos[ent] = (x, 0.0)

    for j, ev in enumerate(events_sorted):
        x = j if n_evt <= 1 else j / (n_evt - 1)
        pos[ev] = (x, 1.0)

    fig = plt.figure(figsize=(10, 4))
    ax = plt.gca()

    scale = fig.get_figwidth() / 10.0
    ent_node_size = 260 * scale**2
    evt_node_size = 420 * scale**2
    label_fontsize = 7 * scale

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="gray", width=1.0, alpha=0.7)

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=entities_sorted,
        node_color="#4C72B0",
        node_size=ent_node_size,
        ax=ax,
    )
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=events_sorted,
        node_color="#FFB570",
        node_size=evt_node_size,
        ax=ax,
    )

    nx.draw_networkx_labels(G, pos, ax=ax, font_size=label_fontsize)

    ax.set_axis_off()
    plt.title("DFlow incidence graph (entities \u2194 events)", fontsize=10 * scale)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# 2. Rhombic lattice overlay with event anchors (Option C)
# ---------------------------------------------------------------------------

HEX_R = 1.0  # hexagon radius


def axial_to_xy(rc: Tuple[int, int], R: float = HEX_R) -> Tuple[float, float]:
    """Axial (q,r) -> Cartesian (x,y) for pointy-top hex coordinates."""
    q, r = rc
    SQRT3 = math.sqrt(3.0)
    x = R * 1.5 * q
    y = R * (SQRT3 * (r + 0.5 * q))
    return x, y


def hex_vertices(x: float, y: float, R: float = HEX_R):
    """Vertices for a pointy-top hex centered at (x,y)."""
    verts = []
    for k in range(6):
        angle = math.radians(60 * k + 30)
        vx = x + R * math.cos(angle)
        vy = y + R * math.sin(angle)
        verts.append((vx, vy))
    return verts


def rhombic_cells(extent=(4, 4)):
    """Return list of axial coords (q,r) inside a rhombic window."""
    qq, rr = extent
    cells = []
    for q in range(-qq, qq + 1):
        for r in range(-rr, rr + 1):
            if abs(q) <= qq and abs(r) <= rr and abs(q + r) <= qq + rr:
                cells.append((q, r))
    return cells


ENTITY_ANCHORS: Dict[str, Tuple[int, int]] = {
    "D1_anchor": (0, 0),
    "D2_anchor": (1, 0),
    "D3_anchor": (0, 1),
    "L1_anchor": (0, 3),
}


def draw_overlay(outpath: str = "hypergraph_overlay.png") -> None:
    """
    Draw a rhombic membrane patch with two rafts (D and L) and overlay
    event-anchor markers.

    Anchor labels are placed *outside* the raft region with thin leader
    lines pointing back to the anchors (Option C).
    """
    cells = rhombic_cells(extent=(4, 4))

    D_raft = {(0, 0), (1, 0), (0, 1), (-1, 1)}
    L_raft = {(0, 3), (1, 3), (0, 4), (-1, 4)}

    fig = plt.figure(figsize=(6, 6))
    ax = plt.gca()
    scale = fig.get_figwidth() / 6.0

    # Background membrane
    for (q, r) in cells:
        x, y = axial_to_xy((q, r))
        verts = hex_vertices(x, y)
        poly = patches.Polygon(
            verts, closed=True, linewidth=0.4, edgecolor="#C0D6F1", facecolor="#E4EFFB"
        )
        ax.add_patch(poly)

    # Rafts
    for (q, r) in D_raft:
        x, y = axial_to_xy((q, r))
        verts = hex_vertices(x, y)
        poly = patches.Polygon(
            verts, closed=True, linewidth=0.7, edgecolor="#B22222", facecolor="#E74C3C"
        )
        ax.add_patch(poly)

    for (q, r) in L_raft:
        x, y = axial_to_xy((q, r))
        verts = hex_vertices(x, y)
        poly = patches.Polygon(
            verts, closed=True, linewidth=0.7, edgecolor="#5B2C6F", facecolor="#AF7AC5"
        )
        ax.add_patch(poly)

    # Anchor markers (small circles)
    for ent, (q, r) in ENTITY_ANCHORS.items():
        x, y = axial_to_xy((q, r))
        ax.scatter(
            x,
            y,
            s=35 * scale**2,
            color="#FFB570",
            edgecolor="black",
            zorder=5,
        )

    # Label positions outside rafts, with leader lines back to anchors.
    # These offsets were eyeballed to keep labels away from the raft hexes.
    label_offsets = {
        "D1_anchor": (-2.3, -0.7),
        "D2_anchor": (2.2, -0.4),
        "D3_anchor": (2.2, 1.2),
        "L1_anchor": (2.2, 3.8),
    }

    for ent, (q, r) in ENTITY_ANCHORS.items():
        x_anchor, y_anchor = axial_to_xy((q, r))
        dx, dy = label_offsets[ent]
        x_label = x_anchor + dx
        y_label = y_anchor + dy

        # Leader line
        ax.plot(
            [x_anchor, x_label],
            [y_anchor, y_label],
            color="black",
            linewidth=0.5,
            zorder=4,
        )

        # Text (no box)
        short = ent.split("_")[0]  # "D1_anchor" -> "D1"
        ax.text(
            x_label,
            y_label,
            short,
            ha="center",
            va="center",
            fontsize=6 * scale,
            zorder=6,
        )

    # R4_drift arrow & label, small and off to the side
    cx = sum(axial_to_xy(p)[0] for p in D_raft) / len(D_raft)
    cy = sum(axial_to_xy(p)[1] for p in D_raft) / len(D_raft)
    tx, ty = cx + 2.5, cy + 0.2

    ax.annotate(
        "R4_drift",
        xy=(cx, cy),
        xytext=(tx, ty),
        arrowprops=dict(arrowstyle="->", lw=1.0),
        fontsize=6 * scale,
    )

    ax.set_aspect("equal")
    ax.set_axis_off()
    plt.title(
        "Rhombic membrane lattice with rafts + event anchors",
        fontsize=9 * scale,
    )
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# 3. Main
# ---------------------------------------------------------------------------

def main():
    G = build_incidence_graph()
    draw_incidence_graph(G, outpath="hypergraph_incidence.png")
    draw_overlay(outpath="hypergraph_overlay.png")
    print("Wrote hypergraph_incidence.png and hypergraph_overlay.png")


if __name__ == "__main__":
    main()
