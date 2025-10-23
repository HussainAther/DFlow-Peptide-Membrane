#!/usr/bin/env python3
import json, math, random, argparse, copy
from pathlib import Path
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon, Rectangle
from collections import defaultdict, deque

# ---- import/copy small bits from your sim ----
SQRT3 = math.sqrt(3.0)
AX_NEI = [(1,0),(1,-1),(0,-1),(-1,0),(-1,1),(0,1)]
TRIAD_SHAPES = [
    [(0,0),(1,0),(0,1)],[(0,0),(0,1),(-1,1)],[(0,0),(-1,1),(-1,0)],
    [(0,0),(-1,0),(0,-1)],[(0,0),(0,-1),(1,-1)],[(0,0),(1,-1),(1,0)]
]
def axial_to_xy(q,r,s): return (s*1.5*q, s*SQRT3*(r+q/2))
def triad_cells(cent, orient):
    aq, ar = cent; offs = TRIAD_SHAPES[orient % 6]
    return [(aq+dq, ar+dr) for (dq,dr) in offs]
def _axial_in_rhombus(q, r, n): return (abs(q) <= n) and (abs(r) <= n) and (abs(q+r) <= n)
def _peptide_length_to_nm(n_res:int) -> float: return 0.15 * n_res

# ---- utility functions for the gallery ----
def get_palette(n: int, name: str = "tab20"):
    cmap = plt.get_cmap(name)
    if hasattr(cmap, "colors"):
        return list(cmap.colors)[:n]
    return [cmap(i / max(1, n - 1)) for i in range(n)]

def load_state(state_json):
    """Loads state, converting string keys back to tuples."""
    with open(state_json) as f: S = json.load(f)
    # Convert string keys "q,r" back to tuple keys (q,r) for the logic functions
    amph_fixed = {}
    for k_str, v in S["amph"].items():
        q, r = map(int, k_str.split(','))
        amph_fixed[(q,r)] = v
    S["amph"] = amph_fixed
    # Convert peptide centers from list to tuple
    for p in S["peptides"].values():
        p["cent"] = tuple(p["cent"])
    return S

def draw_zoom(S, changed_before, changed_after, outpath, title=""):
    n = S["n"]; s = S["hex_radius"]
    amph = S["amph"]; peps = S["peptides"]

    # Determine zoom center: max_pep_q, max_pep_r or (0,0)
    pep_cents = [p["cent"] for p in peps.values()]
    if pep_cents:
        center_q = sum(q for q,_ in pep_cents) // len(pep_cents)
        center_r = sum(r for _,r in pep_cents) // len(pep_cents)
    else:
        center_q, center_r = 0, 0

    # Determine zoom bounds around center
    qmin,qmax = center_q - 5, center_q + 5
    rmin,rmax = center_r - 5, center_r + 5

    def _draw(ax, amph, peps, mark_cells, edge_color, mark_style=''):
        minC,maxC = 10,20
        # Raft colors
        raft_ids = sorted({p["raft"] for p in peps.values()}) if peps else []
        palette = get_palette(max(1, len(raft_ids)), "tab20")
        rid2c = {rid: palette[i] for i,rid in enumerate(raft_ids)}

        # Draw amph cells
        for (q,r), info in amph.items():
            if q<qmin or q>qmax or r<rmin or r>rmax: continue
            x,y = axial_to_xy(q,r,s)
            t = (info["carbon"]-minC)/max(1,(maxC-minC)); t = min(max(t,0.0),1.0)
            base = 0.92 - 0.25*t; color = (base,base,base)
            patch = RegularPolygon((x,y), 6, s, facecolor=color, edgecolor='k', linewidth=0.5)
            ax.add_patch(patch)

        # Draw peptides
        for pid,p in peps.items():
            for (q,r) in triad_cells(p["cent"], p["orient"]):
                if q<qmin or q>qmax or r<rmin or r>rmax: continue
                x,y = axial_to_xy(q,r,s)
                lw = 1.5 if p.get("inside",False) else 0.6
                patch = RegularPolygon((x,y), 6, s,
                                       facecolor=rid2c.get(p["raft"], (0.7,0.2,0.2)),
                                       edgecolor='black', linewidth=lw)
                ax.add_patch(patch)

        # Highlight changed cells
        for (q,r) in mark_cells:
            if q<qmin or q>qmax or r<rmin or r>rmax: continue
            x,y = axial_to_xy(q,r,s)
            patch = RegularPolygon((x,y), 6, s,
                                   fill=False, edgecolor=edge_color, linewidth=2.0)
            ax.add_patch(patch)
            if mark_style == 'X':
                 ax.text(x, y, 'X', ha='center', va='center', fontsize=16, color=edge_color, fontweight='bold')
            if mark_style == '+':
                 ax.text(x, y, '+', ha='center', va='center', fontsize=16, color=edge_color, fontweight='bold')
            if mark_style == '-':
                 ax.text(x, y, '-', ha='center', va='center', fontsize=16, color=edge_color, fontweight='bold')


        # limits
        x_c, y_c = axial_to_xy(center_q, center_r, s)
        ax.set_aspect('equal'); ax.axis('off')
        zoom_range = 7.0 * s # fixed visual range
        ax.set_xlim(x_c - zoom_range, x_c + zoom_range)
        ax.set_ylim(y_c - zoom_range, y_c + zoom_range)


    fig,axs = plt.subplots(1,2, figsize=(8,4))
    style_b = '-' if title.endswith('desorb') else ''
    style_a = '+' if title.endswith('insert') else ''
    _draw(axs[0], S["amph"], S["peptides"], changed_before, "red", mark_style=style_b)
    _draw(axs[1], S["amph_after"], S["peptides_after"], changed_after, "limegreen", mark_style=style_a)

    axs[0].set_title("Before"); axs[1].set_title("After")
    fig.suptitle(title)
    fig.tight_layout()
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=200, bbox_inches="tight"); plt.close(fig)

# ---- apply single events deterministically ----
def pick_demo_amph(S):
    # deterministic: pick amph near (0,0)
    for q in range(-3, 4):
        for r in range(-3, 4):
            if (q,r) in S["amph"] and not S["amph"][(q,r)].get("pep", False):
                return (q,r)
    return None

def pick_demo_peptide(S):
    # deterministic: pick smallest pid
    if not S["peptides"]: return None
    return min(S["peptides"].keys(), key=int)

def find_triad_space(S):
    # Find a spot for insertion, prioritizing next to an existing peptide
    pep_cells = set()
    for p in S["peptides"].values():
        pep_cells.update(triad_cells(p["cent"], p["orient"]))

    free_cells = [c for c, info in S["amph"].items() if not info.get("pep", False)]
    if not free_cells: return None, None, None

    # Try to find a free cell next to an existing peptide
    candidates = []
    for c in free_cells:
        for dq, dr in AX_NEI:
            neighbor = (c[0] + dq, c[1] + dr)
            if neighbor in pep_cells:
                candidates.append(c)
                break

    # If no adjacency, use any free cell
    anchor = random.choice(candidates) if candidates else random.choice(free_cells)

    # Try 6 orientations for the chosen anchor
    for orient in range(6):
        cells = triad_cells(anchor, orient)
        if all(c in S["amph"] and not S["amph"][c].get("pep", False) for c in cells):
            return anchor, orient, cells
    return None, None, None

# --- EVENT FUNCTIONS ---

def event_A1_swap(S):
    # pick center amph (not necessarily peptide center) and a neighbor
    a = pick_demo_amph(S)
    if a is None: return [], [], None
    q,r = a
    for dq,dr in AX_NEI:
        b = (q+dq, r+dr)
        if b in S["amph"] and not S["amph"][b].get("pep", False): # don't swap into a peptide cell
            S2 = copy.deepcopy(S)
            S2["amph"][a]["carbon"], S2["amph"][b]["carbon"] = S2["amph"][b]["carbon"], S2["amph"][a]["carbon"]
            # Swap changes the state of both 'a' and 'b' cells
            return [a,b], [a,b], S2
    return [], [], None

def event_A3_plus_thicken(S):
    a = pick_demo_amph(S)
    if a is None: return [], [], None
    S2 = copy.deepcopy(S)
    S2["amph"][a]["carbon"] = min(S2["amph"][a]["carbon"]+1, 20) # Max C=20
    return [a], [a], S2

def event_A3_minus_thin(S):
    a = pick_demo_amph(S)
    if a is None: return [], [], None
    S2 = copy.deepcopy(S)
    S2["amph"][a]["carbon"] = max(S2["amph"][a]["carbon"]-1, 10) # Min C=10
    return [a], [a], S2

def event_P2_plus_insert(S):
    anchor, orient, cells = find_triad_space(S)
    if not anchor: return [], [], None
    S2 = copy.deepcopy(S)
    pid_new = max(S["peptides"].keys()) + 1 if S["peptides"] else 1
    rid_new = max(p["raft"] for p in S["peptides"].values()) + 1 if S["peptides"] else 1
    nres = 10 # fixed length for demo
    S2["peptides"][pid_new] = {"cent": anchor, "orient": orient, "raft": rid_new, "inside": True, "nres": nres}
    for c in cells:
        S2["amph"][c]["pep"] = True
    return [], cells, S2

def event_P2_minus_desorb(S):
    pid = pick_demo_peptide(S)
    if pid is None: return [], [], None
    p = S["peptides"][pid]
    cells = triad_cells(p["cent"], p["orient"])
    S2 = copy.deepcopy(S)
    del S2["peptides"][pid]
    for c in cells:
        if c in S2["amph"]: S2["amph"][c]["pep"] = False
    return cells, [], S2 # cells were occupied, now they are free

def event_P3_step(S):
    pid = pick_demo_peptide(S)
    if pid is None: return [], [], None
    p = S["peptides"][pid]
    q,r = p["cent"]
    old_cells = triad_cells(p["cent"], p["orient"])

    # Find a free neighbor position
    for dq,dr in AX_NEI:
        new_cent = (q+dq, r+dr)
        new_cells = triad_cells(new_cent, p["orient"])
        # Check if the new cells are free (and exist) and not the old cells
        is_free = True
        for c in new_cells:
            if c not in S["amph"] or (c in S["amph"] and S["amph"][c].get("pep", False) and c not in old_cells):
                is_free = False; break
        if is_free:
            S2 = copy.deepcopy(S)
            S2["peptides"][pid]["cent"] = new_cent
            # Update 'pep' status: old cells are free, new cells are occupied
            for c in old_cells:
                if c in S2["amph"]: S2["amph"][c]["pep"] = False
            for c in new_cells:
                if c in S2["amph"]: S2["amph"][c]["pep"] = True
            
            changed_cells = list(set(old_cells) | set(new_cells))
            return old_cells, new_cells, S2
    return [], [], None

def event_P4_flip(S):
    pid = pick_demo_peptide(S)
    if pid is None: return [], [], None
    p = S["peptides"][pid]
    S2 = copy.deepcopy(S)
    S2["peptides"][pid]["inside"] = (not p.get("inside", False))
    # Changed cells = those triad cells, marked by different color/linewidth after flip
    cells = triad_cells(p["cent"], p["orient"])
    return cells, cells, S2

def _get_connected_components_pids(peptides, amph, n):
    ids = list(peptides.keys())
    if not ids: return []
    idx = {pid:i for i,pid in enumerate(ids)}
    adj = [[] for _ in ids]
    cell2pids = defaultdict(list)
    for pid, p in peptides.items():
        for c in triad_cells(p["cent"], p["orient"]):
            cell2pids[c].append(pid)
    for pid, p in peptides.items():
        seen = set()
        for (q,r) in triad_cells(p["cent"], p["orient"]):
            for dq,dr in AX_NEI:
                nb = (q+dq, r+dr)
                if nb in cell2pids:
                    for pid2 in cell2pids[nb]:
                        if pid2 != pid and pid2 not in seen:
                            adj[idx[pid]].append(idx[pid2])
                            adj[idx[pid2]].append(idx[pid])
                            seen.add(pid2)
    comps = []
    vis = [False]*len(ids)
    for i in range(len(ids)):
        if vis[i]: continue
        q = deque([i]); vis[i]=True; group=[ids[i]]
        while q:
            u = q.popleft()
            for v in adj[u]:
                if not vis[v]:
                    vis[v]=True; q.append(v); group.append(ids[v])
        comps.append(group)
    return comps

def event_R1_plus_assoc(S):
    # Find two rafts that are currently touching but have different raft IDs
    components = _get_connected_components_pids(S["peptides"], S["amph"], S["n"])
    if len(components) < 2: return [], [], None

    cell2pid = {}
    for pid, p in S["peptides"].items():
        for c in triad_cells(p["cent"], p["orient"]):
            cell2pid[c] = pid

    touching_rafts = set()
    for pid1, p1 in S["peptides"].items():
        r1 = p1["raft"]
        for c in triad_cells(p1["cent"], p1["orient"]):
            for dq,dr in AX_NEI:
                nb = (c[0]+dq, c[1]+dr)
                if nb in cell2pid:
                    pid2 = cell2pid[nb]
                    r2 = S["peptides"][pid2]["raft"]
                    if r1 != r2:
                        touching_rafts.add(tuple(sorted((r1, r2))))

    if not touching_rafts: return [], [], None

    rA, rB = random.choice(list(touching_rafts))
    target_raft, source_raft = (rA, rB) if rA < rB else (rB, rA)

    S2 = copy.deepcopy(S)
    pids_moved = []
    for pid, p in S2["peptides"].items():
        if p["raft"] == source_raft:
            p["raft"] = target_raft
            pids_moved.append(pid)

    # Highlight all cells of the peptides involved in the merge
    changed_cells = []
    for pid, p in S["peptides"].items():
        if p["raft"] in (rA, rB):
            changed_cells.extend(triad_cells(p["cent"], p["orient"]))

    return changed_cells, changed_cells, S2

def event_R1_minus_dissoc(S):
    components = _get_connected_components_pids(S["peptides"], S["amph"], S["n"])
    comps = [c for c in components if len(c) >= 2]
    if not comps: return [], [], None

    # Pick the largest component to split deterministically
    comp_to_split = max(comps, key=len)
    
    # Split: move roughly half to a new raft ID
    new_raft = max(p["raft"] for p in S["peptides"].values()) + 1
    random.seed(123) # deterministic split
    to_move = set(random.sample(comp_to_split, k=len(comp_to_split)//2))

    S2 = copy.deepcopy(S)
    for pid in to_move:
        S2["peptides"][pid]["raft"] = new_raft

    # Highlight all cells of the peptides in the component that was split
    changed_cells = []
    for pid, p in S["peptides"].items():
        if pid in comp_to_split:
            changed_cells.extend(triad_cells(p["cent"], p["orient"]))

    return changed_cells, changed_cells, S2

# --- REGISTER EVENTS ---
EVENTS = {
    "A1_swap": event_A1_swap,
    "A3_plus_thicken": event_A3_plus_thicken,
    "A3_minus_thin": event_A3_minus_thin,
    "P2_plus_insert": event_P2_plus_insert,
    "P2_minus_desorb": event_P2_minus_desorb,
    "P3_step": event_P3_step,
    "P4_flip": event_P4_flip,
    "R1_plus_assoc": event_R1_plus_assoc,
    "R1_minus_dissoc": event_R1_minus_dissoc,
}

def main():
    ap = argparse.ArgumentParser("Generate before/after case figures from a saved state JSON")
    ap.add_argument("--STATE_JSON", required=True, help="exported state json with amph, peptides, n, hex_radius")
    ap.add_argument("--CASE", required=True, choices=list(EVENTS.keys()))
    ap.add_argument("--OUT", required=True)
    args = ap.parse_args()

    # Seed for deterministic event selection (e.g., which peptide to pick)
    random.seed(42)

    S = load_state(args.STATE_JSON)
    
    # Run event function on a copy of the state
    before_changed, after_changed, S2 = EVENTS[args.CASE](copy.deepcopy(S))
    
    if S2 is None:
        print(f"No-op for {args.CASE} (no suitable target found in state)."); return

    # Attach after-state for drawing
    S["amph_after"] = S2["amph"]
    S["peptides_after"] = S2["peptides"]
    
    draw_zoom(S, before_changed, after_changed, args.OUT, title=args.CASE)
    print(f"Wrote {args.OUT}")

if __name__ == "__main__":
    main()
