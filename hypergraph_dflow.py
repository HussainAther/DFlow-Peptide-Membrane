#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hypergraph_dflow.py

Rebuild of the DFlow "hypergraph" / incidence-graph utilities,
with explicit inclusion of the R4_drift event.

This script does two things:

1. Builds an incidence graph between "entities" (blue, bottom) and
   "events" (orange, top) and saves it as:
       hypergraph_incidence.png

2. Overlays a small set of event anchors on a rhombic hexagonal
   membrane patch with two rafts (D and L) and saves:
       hypergraph_overlay.png

Dependencies:
    numpy
    matplotlib
    networkx   # pip install networkx

You can drop this file in the repo root and run:

    python hypergraph_dflow.py

to regenerate the two figures.
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
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

# The design here is deliberately simple: we hard-code a small set of
# "entities" and "events" that are meant to be archetypal, not exhaustive.
#
# Entities are roughly:
#   D*_anchor  -> domains in a D-raft (central raft)
#   L*_anchor  -> domains in a neighboring L-raft
#
# Events are the DFlow event types we care about for the conceptual
# picture. We include R4_drift explicitly and connect it to all raft
# anchors, since it moves the entire raft domain as a rigid body.

ENTITY_NODES: List[str] = [
    "D1_anchor",
    "D2_anchor",
    "D3_anchor",
    "L1_anchor",
]

EVENT_TO_ENTITIES: Dict[str, List[str]] = {
    # Amphiphile thickness / swap events
    "A1_forward": ["D1_anchor", "D2_anchor", "D3_anchor"],
    "A1_reverse": ["D1_anchor", "D2_anchor", "D3_anchor"],

    # Peptide insertion / desorption at the D raft
    "P2_insert": ["D1_anchor", "D2_anchor", "D3_anchor"],
    "P2_desorb": ["D1_anchor", "D2_anchor", "D3_anchor"],

    # Raft association / dissociation between neighboring domains
    "R1_assoc": ["D2_anchor", "L1_anchor"],
    "R1_dissoc": ["D2_anchor", "L1_anchor"],

    # Optional fusion (kept as an event type for completeness)
    "R2_fusion": ["D1_anchor", "L1_anchor"],

    # Peptide flips between L and D metadata states
    "P4_flip_LtoD": ["L1_anchor"],
    "P4_flip_DtoL": ["D1_anchor"],

    # R4: size-dependent 2D Stokes drift of the entire raft domain
    # We connect it to all raft anchors, since the drift is a rigid-body
    # translation of the raft as a whole.
    "R4_drift": ["D1_anchor", "D2_anchor", "D3_anchor", "L1_anchor"],
}


def build_incidence_graph() -> "nx.Graph":
    """
    Build a bipartite graph: entities (bipartite=0) and events (bipartite=1).

    Returns
    -------
    G : networkx.Graph
        Graph with node attributes:
            - kind: "entity" or "event"
            - bipartite: 0 or 1
    """
    G = nx.Graph()

    # Add entity nodes
    for ent in ENTITY_NODES:
        G.add_node(ent, kind="entity", bipartite=0)

    # Add event nodes and incidence edges
    for ev, ents in EVENT_TO_ENTITIES.items():
        G.add_node(ev, kind="event", bipartite=1)
        for ent in ents:
            if ent not in G:
                raise ValueError(f"Unknown entity '{ent}' referenced by event '{ev}'")
            G.add_edge(ev, ent)

    return G


def draw_incidence_graph(G: "nx.Graph", outpath: str = "hypergraph_incidence.png") -> None:
    """
    Draw the incidence graph with events on top and entities at the bottom.

    The layout is deliberately simple and deterministic: entities are placed
    along a line at y=0, events along another at y=1, and we rescale x
    coordinates so edges have a nice spread.
    """
    # Separate bipartite sets
    entities = [n for n, d in G.nodes(data=True) if d.get("kind") == "entity"]
    events = [n for n, d in G.nodes(data=True) if d.get("kind") == "event"]

    # Deterministic ordering
    entities_sorted = sorted(entities)
    events_sorted = sorted(events)

    pos: Dict[str, Tuple[float, float]] = {}

    # Place entities on bottom line
    n_ent = len(entities_sorted)
    for i, ent in enumerate(entities_sorted):
        x = i if n_ent <= 1 else i / (n_ent - 1)
        pos[ent] = (x, 0.0)

    # Place events on top line
    n_evt = len(events_sorted)
    for j, ev in enumerate(events_sorted):
        x = j if n_evt <= 1 else j / (n_evt - 1)
        pos[ev] = (x, 1.0)

    plt.figure(figsize=(10, 4))
    ax = plt.gca()

    # Draw edges
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="gray", width=1.0, alpha=0.7)

    # Draw nodes
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=entities_sorted,
        node_color="#4C72B0",  # blue
        node_size=400,
        ax=ax,
    )
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=events_sorted,
        node_color="#FFB570",  # orange
        node_size=600,
        ax=ax,
    )

    # Labels
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8)

    ax.set_axis_off()
    plt.title("DFlow incidence graph (entities \u2194 events)", fontsize=12)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# 2. Rhombic lattice overlay with event anchors
# ---------------------------------------------------------------------------

HEX_R = 1.0  # hexagon radius for plotting (visual only)


def axial_to_xy(rc: Tuple[int, int], R: float = HEX_R) -> Tuple[float, float]:
    """
    Axial (q,r) -> cartesian (x,y) for pointy-top hex coordinates.

    x = R * (3/2) * q
    y = R * (sqrt(3) * r + sqrt(3)/2 * q)
    """
    q, r = rc
    SQRT3 = math.sqrt(3.0)
    x = R * 1.5 * q
    y = R * (SQRT3 * (r + 0.5 * q))
    return x, y


def hex_vertices(x: float, y: float, R: float = HEX_R):
    """Vertices for a pointy-top hex centered at (x,y)."""
    SQRT3 = math.sqrt(3.0)
    # Start at 30 degrees for pointy-top
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


# Anchor coordinates for entities (purely illustrative)
ENTITY_ANCHORS: Dict[str, Tuple[int, int]] = {
    "D1_anchor": (0, 0),
    "D2_anchor": (1, 0),
    "D3_anchor": (0, 1),
    "L1_anchor": (0, 3),  # second raft, above the D raft
}


def draw_overlay(outpath: str = "hypergraph_overlay.png") -> None:
    """
    Draw a rhombic membrane patch with two rafts (D and L) and overlay
    a few event-anchor markers, including R4_drift.
    """
    cells = rhombic_cells(extent=(4, 4))

    # Define two rafts as compact clusters around their anchors
    D_raft = {(0, 0), (1, 0), (0, 1), (-1, 1)}
    L_raft = {(0, 3), (1, 3), (0, 4), (-1, 4)}

    plt.figure(figsize=(8, 4))
    ax = plt.gca()

    # Draw background membrane
    for (q, r) in cells:
        x, y = axial_to_xy((q, r))
        verts = hex_vertices(x, y)
        poly = patches.Polygon(
            verts, closed=True, linewidth=0.5, edgecolor="#C0D6F1", facecolor="#E4EFFB"
        )
        ax.add_patch(poly)

    # Draw rafts
    for (q, r) in D_raft:
        x, y = axial_to_xy((q, r))
        verts = hex_vertices(x, y)
        poly = patches.Polygon(
            verts, closed=True, linewidth=0.8, edgecolor="#B22222", facecolor="#E74C3C"
        )
        ax.add_patch(poly)

    for (q, r) in L_raft:
        x, y = axial_to_xy((q, r))
        verts = hex_vertices(x, y)
        poly = patches.Polygon(
            verts, closed=True, linewidth=0.8, edgecolor="#5B2C6F", facecolor="#AF7AC5"
        )
        ax.add_patch(poly)

    # Event anchors: place small circles at entity anchors, one per entity
    for ent, (q, r) in ENTITY_ANCHORS.items():
        x, y = axial_to_xy((q, r))
        ax.scatter(x, y, s=60, color="#FFB570", edgecolor="black", zorder=5)
        ax.text(
            x,
            y + 0.3,
            ent,
            ha="center",
            va="bottom",
            fontsize=7,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=0.5),
            zorder=6,
        )

    # R4_drift arrow: from D raft centroid to a drifted position
    dq, dr = 1, 0  # drift direction (axial)
    cx = sum(axial_to_xy(p)[0] for p in D_raft) / len(D_raft)
    cy = sum(axial_to_xy(p)[1] for p in D_raft) / len(D_raft)
    cx2, cy2 = axial_to_xy((1, 0))  # a representative target cell

    ax.annotate(
        "R4_drift",
        xy=(cx, cy),
        xytext=(cx2 + 0.5, cy2 - 0.5),
        arrowprops=dict(arrowstyle="->", lw=1.2),
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=0.5),
    )

    ax.set_aspect("equal")
    ax.set_axis_off()
    plt.title("Rhombic membrane lattice with rafts + event anchors", fontsize=12)
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
