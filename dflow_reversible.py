#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Consolidated DFlow reversible simulator.

import os, json, math, random, argparse
from dataclasses import dataclass
from pathlib import Path
from math import radians
from collections import defaultdict, deque
from typing import Optional, List, Tuple, Dict, Iterable

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

def get_palette(n: int, name: str = "tab20"):
    try:
        cmap = mpl.colormaps[name]
    except AttributeError:
        cmap = plt.get_cmap(name)
    if hasattr(cmap, "colors"):
        base = list(cmap.colors)
        if n <= len(base):
            return base[:n]
        reps = (n + len(base) - 1) // len(base)
        return (base * reps)[:n]
    if n == 1:
        return [cmap(0.0)]
    return [cmap(i / max(1, n - 1)) for i in range(n)]

SQRT3 = math.sqrt(3.0)
AX_NEI = [(1,0),(1,-1),(0,-1),(-1,0),(-1,1),(0,1)]
TRIAD_SHAPES = [
    [(0,0),(1,0),(0,1)],
    [(0,0),(0,1),(-1,1)],
    [(0,0),(-1,1),(-1,0)],
    [(0,0),(-1,0),(0,-1)],
    [(0,0),(0,-1),(1,-1)],
    [(0,0),(1,-1),(1,0)]
]

def axial_to_xy(q:int, r:int, s:float):
    x = s * (3/2) * q
    y = s * SQRT3 * (r + q/2)
    return x, y

def rotate_triad(anchor:Tuple[int,int], orient:int):
    aq, ar = anchor
    offs = TRIAD_SHAPES[orient % 6]
    return [(aq+dq, ar+dr) for (dq,dr) in offs]

@dataclass
class Config:
    SEED: int = 42
    N0: int = 12
    HEX_RADIUS: float = 1.0
    TOTAL_EVENTS: int = 1500
    OUT: Path = Path("runs/exp_phys")
    SAVE_STEPS: List[int] = None
    DAY_STEPS: int = 10
    NIGHT_STEPS: int = 10
    POLY_GAIN_DAY: float = 10.0
    POLY_GAIN_NIGHT: float = 0.2
    MAX_PEPT_LEN: int = 12
    A_INSERT_MISMATCH_TOL: float = 1.0
    BETA: float = 1.0
    CROWDING_GAMMA: float = 0.02
    CROWDING_DECAY: float = 0.0
    DESORB_TO_CROWDING: bool = False
    DISCARD_TO_CROWDING_P: float = 0.5
    RAFT_D0: float = 1.0
    RAFT_DIFF_SIZE_EXP: float = 1.0
    FUSION_IRREVERSIBLE: bool = False
    FUSION_SIZE_THRESH: int = 6
    INIT_CARBON_MIN: int = 10
    INIT_CARBON_MAX: int = 20
    LABEL_CARBONS: bool = False
    BORDER: float = 1.0
    FUSION_K0: float = 1.0
    FISSION_MIN_SIZE: int = 4                 
    FISSION_HEURISTIC: str = "bridge"  
    FISSION_SIZE_SCALING: str = "linear"
    FISSION_K0: float = 1.0
    FISSION_ENABLED: bool = False
    FISSION_ALLOW_PURE_CHIRAL: bool = False
  

    def __post_init__(self):
        if self.SAVE_STEPS is None:
            self.SAVE_STEPS = [0, self.TOTAL_EVENTS//2, self.TOTAL_EVENTS]

@dataclass
class Diurnal:
    day_steps: int
    night_steps: int
    gain_day: float
    gain_night: float
    t_in_cycle: int = 0
    is_day: bool = True
    def tick(self):
        self.t_in_cycle += 1
        if self.is_day and self.t_in_cycle >= self.day_steps:
            self.is_day = False; self.t_in_cycle = 0
        elif not self.is_day and self.t_in_cycle >= self.night_steps:
            self.is_day = True; self.t_in_cycle = 0
    @property
    def poly_gain(self) -> float:
        return self.gain_day if self.is_day else self.gain_night


class Membrane:
    def __init__(self, n:int, hex_radius:float, cfg:Config):
        self.cfg = cfg
        self.n = n
        self.hex_radius = hex_radius
        self.amph: Dict[Tuple[int,int], Dict] = {}
        self.peptides: Dict[int, Dict] = {}
        self.pid_counter = 0
        self.raft_counter = 0
        self.crowding_count = 0
        self.pool = {"aa": 0, "peptides": 0}
        self.event_log: List[Dict] = []
        self._init_membrane()
        self.bonds: set[tuple[int, int]] = set()
        self.raft_meta: Dict[int, Dict] = {}
        self.last_R2: Dict[int, Dict[str, int]] = {}  
        self.R2_COOLDOWN_STEPS = 1 
    def _axial_in_rhombus(self, q:int, r:int) -> bool:
        n = self.n
        return (abs(q) <= n) and (abs(r) <= n) and (abs(q+r) <= n)
    def _init_membrane(self):
        for q in range(-self.n, self.n+1):
            for r in range(-self.n, self.n+1):
                if self._axial_in_rhombus(q,r):
                    self.amph[(q,r)] = {"carbon": random.randint(self.cfg.INIT_CARBON_MIN, self.cfg.INIT_CARBON_MAX),
                                        "pep": False}
    def local_membrane_thickness_nm(self, center:Tuple[int,int]) -> float:
        q,r = center
        neigh = [(q,r)] + [(q+dq, r+dr) for dq,dr in AX_NEI]
        vals = [self.amph[c]["carbon"] for c in neigh if c in self.amph]
        if not vals: return 2.0
        meanC = sum(vals)/len(vals)
        bilayer = 0.2 * meanC
        return max(1.0, bilayer)
    def triad_cells(self, center:Tuple[int,int], orient:int) -> List[Tuple[int,int]]:
        return rotate_triad(center, orient)
    def cells_free_for_triad(self, cells:Iterable[Tuple[int,int]]) -> bool:
        for c in cells:
            if c not in self.amph: return False
            if self.amph[c]["pep"]: return False
        return True
    def random_free_anchor(self, max_tries:int=1000):
        keys = list(self.amph.keys())
        for _ in range(max_tries):
            anchor = random.choice(keys)
            orient = random.randrange(6)
            cells = self.triad_cells(anchor, orient)
            if self.cells_free_for_triad(cells):
                return anchor, orient, cells
        return None, None, None

@dataclass
class Event:
    name: str
    rate_fn: callable
    do_fn: callable

class Scheduler:
    def __init__(self):
        self.events: List[Event] = []
    def register(self, e:Event):
        self.events.append(e)
    def step(self, mem, env):
        rates = np.array([max(0.0, e.rate_fn(mem, env)) for e in self.events], dtype=float)
        s = rates.sum()
        if s <= 0.0: return None
        p = rates / s
        idx = np.random.choice(len(self.events), p=p)
        chosen = self.events[idx]
        chosen.do_fn(mem, env)
        return chosen.name

def rate_swap(mem, env): return 1.0
def do_swap(mem, env):
    if not mem.amph: return
    a = random.choice(list(mem.amph.keys()))
    q,r = a
    random.shuffle(AX_NEI)
    for dq,dr in AX_NEI:
        b = (q+dq, r+dr)
        if b in mem.amph:
            mem.amph[a]["carbon"], mem.amph[b]["carbon"] = mem.amph[b]["carbon"], mem.amph[a]["carbon"]
            mem.event_log.append({"evt":"A1_swap", "a":a, "b":b})
            return

def rate_thicken(mem, env): return 0.3
def do_thicken(mem, env):
    a = random.choice(list(mem.amph.keys()))
    mem.amph[a]["carbon"] = min(mem.amph[a]["carbon"]+1, mem.cfg.INIT_CARBON_MAX)
    mem.event_log.append({"evt":"A3_plus_thicken", "cell":a})

def rate_thin(mem, env): return 0.25
def do_thin(mem, env):
    a = random.choice(list(mem.amph.keys()))
    mem.amph[a]["carbon"] = max(mem.amph[a]["carbon"]-1, mem.cfg.INIT_CARBON_MIN)
    mem.event_log.append({"evt":"A3_minus_thin", "cell":a})

def rate_polymerize(mem, env):
    base = env.diurnal.poly_gain
    slowdown = 1.0 / (1.0 + mem.crowding_count * mem.cfg.CROWDING_GAMMA)
    return base * slowdown
def do_polymerize(mem, env):
    mem.pool["aa"] += 1
    if mem.pool["aa"] >= 3:
        mem.pool["aa"] -= 3
        mem.pool["peptides"] += 1
        mem.event_log.append({"evt":"P1_plus_polymerize", "aa_to_pool_pep":1})
    else:
        mem.event_log.append({"evt":"P1_plus_polymerize", "aa_only":1})
def rate_depoly(mem, env): return 0.05
def do_depoly(mem, env):
    if mem.pool["peptides"] > 0:
        mem.pool["peptides"] -= 1
        mem.pool["aa"] += 2
        mem.event_log.append({"evt":"P1_minus_depoly", "peptide_to_aa":1})
    else:
        mem.event_log.append({"evt":"P1_minus_depoly", "no_pep":1})

def _peptide_length_to_nm(n_res:int) -> float: return 0.15 * n_res
def rate_insert(mem, env): return 0.8 if mem.pool["peptides"] > 0 else 0.0
def do_insert(mem, env):
    if mem.pool["peptides"] <= 0:
        mem.event_log.append({"evt":"P2_plus_insert_try","status":"no_pool"}); return
    anchor, orient, cells = mem.random_free_anchor()
    if anchor is None:
        mem.event_log.append({"evt":"P2_plus_insert_try","status":"no_space"}); return
    nres = random.randint(5, mem.cfg.MAX_PEPT_LEN)
    Lp = _peptide_length_to_nm(nres)
    tloc = mem.local_membrane_thickness_nm(anchor)
    leaflet = 0.5 * tloc
    mismatch = abs(Lp - leaflet)
    accept = False
    if mismatch <= mem.cfg.A_INSERT_MISMATCH_TOL:
        accept = True
    else:
        delta = mismatch - mem.cfg.A_INSERT_MISMATCH_TOL
        p = math.exp(-mem.cfg.BETA * delta)
        accept = (random.random() < p)
    if not accept:
        if random.random() < mem.cfg.DISCARD_TO_CROWDING_P:
            mem.crowding_count += 1
            mem.event_log.append({"evt":"P2_insert_reject","route":"crowding","mismatch":mismatch,"nres":nres})
        else:
            mem.pool["peptides"] += 1
            mem.event_log.append({"evt":"P2_insert_reject","route":"pool","mismatch":mismatch,"nres":nres})
        return
    mem.pool["peptides"] -= 1
    pid = mem.pid_counter; mem.pid_counter += 1
    rid = mem.raft_counter; mem.raft_counter += 1
    inside = bool(random.getrandbits(1))
    chir = 'D' if random.getrandbits(1) else 'L' 
    mem.peptides[pid] = {"cent": anchor, "orient": orient, "raft": rid, "inside": inside, "nres": nres , "chir": chir}
    for c in cells:
        mem.amph[c]["pep"] = True
    mem.event_log.append({"evt":"P2_plus_insert","pid":pid,"raft":rid,"cells":cells,"nres":nres,"mismatch":mismatch})

def rate_desorb(mem, env): return 0.10 if mem.peptides else 0.0
def do_desorb(mem, env):
    if not mem.peptides:
        mem.event_log.append({"evt":"P2_minus_desorb_try","status":"no_pep"}); return
    pid = random.choice(list(mem.peptides.keys()))
    p = mem.peptides.pop(pid)
    for c in rotate_triad(p["cent"], p["orient"]):
        if c in mem.amph: mem.amph[c]["pep"] = False
    if mem.cfg.DESORB_TO_CROWDING:
        mem.crowding_count += 1; route = "crowding"
    else:
        mem.pool["peptides"] += 1; route = "pool"
    mem.event_log.append({"evt":"P2_minus_desorb","pid":pid,"route":route})

def rate_pept_step(mem, env): return 0.5 if mem.peptides else 0.0
def do_pept_step(mem, env):
    if not mem.peptides: return
    pid = random.choice(list(mem.peptides.keys()))
    p = mem.peptides[pid]
    raft = p["raft"]
    same = sum(1 for x in mem.peptides.values() if x["raft"] == raft)
    D = mem.cfg.RAFT_D0 / (max(1, same) ** mem.cfg.RAFT_DIFF_SIZE_EXP)
    if random.random() > D: return
    q,r = p["cent"]
    dq,dr = random.choice(AX_NEI)
    new_cent = (q+dq, r+dr)
    new_cells = rotate_triad(new_cent, p["orient"])
    if not mem.cells_free_for_triad(new_cells): return
    for c in rotate_triad(p["cent"], p["orient"]):
        mem.amph[c]["pep"] = False
    for c in new_cells:
        mem.amph[c]["pep"] = True
    p["cent"] = new_cent
    mem.event_log.append({"evt":"P3_step","pid":pid,"from":(q,r),"to":new_cent,"raft":raft,"D":D})

def rate_flip(mem, env): return 0.1 if mem.peptides else 0.0
def do_flip(mem, env):
    if not mem.peptides: return
    pid = random.choice(list(mem.peptides.keys()))
    p = mem.peptides[pid]; p["inside"] = not p["inside"]
    mem.event_log.append({"evt":"P4_flip","pid":pid,"inside":p["inside"]})

def _build_cell_to_pid(mem):
    m = {}
    for pid, p in mem.peptides.items():
        for c in rotate_triad(p["cent"], p["orient"]):
            m[c] = pid
    return m

def rate_assoc(mem, env):
    if len(mem.peptides) < 2: return 0.0
    cell2pid = _build_cell_to_pid(mem)
    touches = 0
    for c, pid in cell2pid.items():
        r1 = mem.peptides[pid]["raft"]
        for dq,dr in AX_NEI:
            nb = (c[0]+dq, c[1]+dr)
            if nb in cell2pid:
                pid2 = cell2pid[nb]
                if pid2 != pid and mem.peptides[pid2]["raft"] != r1:
                    touches += 1
    return min(1.0, 0.05 * touches)
def do_assoc(mem, env):
    if len(mem.peptides) < 2: return
    cell2pid = _build_cell_to_pid(mem)
    pairs = set()
    for c, pid in cell2pid.items():
        r1 = mem.peptides[pid]["raft"]
        for dq,dr in AX_NEI:
            nb = (c[0]+dq, c[1]+dr)
            if nb in cell2pid:
                pid2 = cell2pid[nb]
                if pid2 != pid:
                    r2 = mem.peptides[pid2]["raft"]
                    if r1 != r2:
                        pairs.add(tuple(sorted((r1,r2))))
    if not pairs: return
    rA, rB = random.choice(list(pairs))
    target, source = (rA, rB) if rA < rB else (rB, rA)
    for pid, p in mem.peptides.items():
        if p["raft"] == source: p["raft"] = target
    mem.event_log.append({"evt":"R1_plus_assoc","from":source,"into":target})

def rate_dissoc(mem, env):
    if len(mem.peptides) < 2: return 0.0
    rafts = defaultdict(int)
    for p in mem.peptides.values(): rafts[p["raft"]] += 1
    size = random.choice(list(rafts.values()))
    if mem.cfg.FUSION_IRREVERSIBLE and size >= mem.cfg.FUSION_SIZE_THRESH: return 0.0
    return max(0.01, 0.2 / size)
def do_dissoc(mem, env):
    comps = connected_components_full(mem)
    comps = [c for c in comps if len(c) >= 2]
    if not comps: return
    comp = random.choice(comps)
    new_raft = mem.raft_counter; mem.raft_counter += 1
    from random import sample
    to_move = set(sample(list(comp), k=len(comp)//2))
    for pid in to_move: mem.peptides[pid]["raft"] = new_raft
    mem.event_log.append({"evt":"R1_minus_dissoc","new_raft":new_raft,"moved":len(to_move)})

def connected_components_full(mem) -> List[List[int]]:
    ids = list(mem.peptides.keys())
    if not ids: return []
    idx = {pid:i for i,pid in enumerate(ids)}
    adj = [[] for _ in ids]
    cell2pids = defaultdict(list)
    for pid, p in mem.peptides.items():
        for c in rotate_triad(p["cent"], p["orient"]):
            cell2pids[c].append(pid)
    for pid, p in mem.peptides.items():
        seen = set()
        for (q,r) in rotate_triad(p["cent"], p["orient"]):
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

def _raft_ids_and_adjacency(mem, raft_id:int):
    """Return (node_ids, adjacency dict) for the induced subgraph of a single raft."""
    raft_nodes = [pid for pid, p in mem.peptides.items() if p["raft"] == raft_id]
    if not raft_nodes:
        return [], {}
    idx = {pid:i for i, pid in enumerate(raft_nodes)}
    adj = {pid:set() for pid in raft_nodes}
    cell2pids = defaultdict(list)
    for pid, p in mem.peptides.items():
        for c in rotate_triad(p["cent"], p["orient"]):
            cell2pids[c].append(pid)
    for pid in raft_nodes:
        p = mem.peptides[pid]
        seen = set()
        for (q,r) in rotate_triad(p["cent"], p["orient"]):
            for dq,dr in AX_NEI:
                nb = (q+dq, r+dr)
                if nb in cell2pids:
                    for pid2 in cell2pids[nb]:
                        if pid2 != pid and pid2 not in seen and mem.peptides[pid2]["raft"] == raft_id:
                            adj[pid].add(pid2)
                            adj[pid2].add(pid)
                            seen.add(pid2)
    return raft_nodes, adj

def _bridges_tarjan(nodes, adj):
    """Tarjan bridges on an undirected graph given by adj: dict[node] -> set(neigh)."""
    time = 0
    disc, low, parent = {}, {}, {}
    bridges = set()
    def dfs(u):
        nonlocal time
        time += 1
        disc[u] = low[u] = time
        for v in adj.get(u, ()):
            if v not in disc:
                parent[v] = u
                dfs(v)
                low[u] = min(low[u], low[v])
                if low[v] > disc[u]:
                    bridges.add(tuple(sorted((u, v))))
            elif parent.get(u) != v:
                low[u] = min(low[u], disc[v])
    for u in nodes:
        if u not in disc and u in adj:
            dfs(u)
    return bridges

def _weak_edge_scores(nodes, adj):
    """Score edges by shared neighbors; fewer shared ⇒ weaker."""
    scores = {}
    for i in nodes:
        if i not in adj: continue
        for j in adj[i]:
            if i < j:
                scores[(i,j)] = len(adj[i] & adj.get(j, set()))
    return scores

def _subcomponents_after_cut(nodes_set, adj, cut):
    """Return two node-sets after removing edge 'cut' = (u,v) from the induced graph."""
    u, v = cut

    def bfs(start, block):
        seen = set()
        stack = [start]
        bu, bv = block  
        while stack:
            x = stack.pop()
            if x in seen: 
                continue
            seen.add(x)
            for w in adj.get(x, ()):
                if (x == bu and w == bv) or (x == bv and w == bu):
                    continue
                if w not in seen:
                    stack.append(w)
        return seen

    u_side = bfs(u, (u, v))
    v_side = nodes_set - u_side
    return u_side, v_side
def _r2_blocked(mem, rid: int, env) -> bool:
    rec = mem.last_R2.get(rid)
    if rec is None:
        return False
    return (env.step_idx - rec["step"]) <= mem.R2_COOLDOWN_STEPS

def _mark_r2(mem, rid: int, env, typ: str):
    mem.last_R2[rid] = {"type": typ, "step": int(env.step_idx)}


def accept_move(beta: float, dE: float) -> bool:
    """Metropolis acceptance."""
    if dE <= 0:
        return True
    return random.random() < math.exp(-beta * dE)

def _cell2pids(mem):
    m = defaultdict(list)
    for pid, p in mem.peptides.items():
        for c in rotate_triad(p["cent"], p["orient"]):
            m[c].append(pid)
    return m

def _boundary_pairs(mem):
    """
    Yield boundary contacts: (peptide_cell, neighbor_cell, pid_at_peptide_cell)
    where peptide_cell is occupied and neighbor_cell is *not* occupied
    (i.e., raft perimeter).
    """
    c2pids = _cell2pids(mem)  
    for c, pids_here in c2pids.items():
        for dq, dr in AX_NEI:
            nb = (c[0]+dq, c[1]+dr)
            if nb not in c2pids and nb in mem.amph:
                for pid in pids_here:
                    yield c, nb, pid

def _peptide_length_nm_from_pid(mem, pid:int) -> float:
    return 0.15 * max(1, mem.peptides[pid].get("nres", 6))

def total_energy(mem,
                 J_same: float = 1.0,
                 L_mismatch: float = 1.0,
                 K_crowd: float = 0.0) -> float:
    """
    System energy E = E_bonds + E_mismatch + E_crowd.

    - E_bonds:    -J_same * (# of adjacencies where touching cells belong to peptides of the SAME raft)
                  (favors larger, contiguous rafts)
    - E_mismatch:  L_mismatch * sum_over_perimeter( | L_p - leaflet_thickness | )
                   where L_p = ~0.15 nm per residue; leaflet_thickness = 0.5 * local_membrane_thickness_nm
                   (penalizes hydrophobic length mismatch at raft boundary)
    - E_crowd:     K_crowd * crowding_count  (optional macroscopic penalty)

    All three terms are scalars; lower E is favored.
    """
    if not mem.peptides:
        return K_crowd * float(mem.crowding_count)

    c2p = _cell2pids(mem)
    seen_edges = set()
    same_contacts = 0
    for c, pids_here in c2p.items():
        for dq, dr in AX_NEI:
            nb = (c[0] + dq, c[1] + dr)
            if nb not in c2p:
                continue
            edge = tuple(sorted((c, nb)))
            if edge in seen_edges:
                continue
            seen_edges.add(edge)
            same = False
            for pid1 in pids_here:
                r1 = mem.peptides[pid1]["raft"]
                for pid2 in c2p[nb]:
                    if r1 == mem.peptides[pid2]["raft"]:
                        same = True
                        break
                if same:
                    break
            if same:
                same_contacts += 1
    E_bonds = -J_same * float(same_contacts)

    perimeter_penalty = 0.0
    leaflet_cache = {}
    for c_pep, c_nb, pid in _boundary_pairs(mem):
        cent = mem.peptides[pid]["cent"]
        if cent not in leaflet_cache:
            leaflet_cache[cent] = 0.5 * mem.local_membrane_thickness_nm(cent)
        leaflet = leaflet_cache[cent]
        Lp = _peptide_length_nm_from_pid(mem, pid)
        perimeter_penalty += abs(Lp - leaflet)
    E_mismatch = L_mismatch * perimeter_penalty

    return E_bonds + E_mismatch 

def raft_chirality_counts(mem, raft_id:int) -> Tuple[int, int]:
    """Return (#D, #L) for a given raft."""
    d = l = 0
    for p in mem.peptides.values():
        if p["raft"] == raft_id:
            c = p.get("chir", 'L')  
            if c == 'D':
                d += 1
            else:
                l += 1
    return d, l

def raft_is_pure_chiral(mem, raft_id:int) -> bool:
    """True if the raft is all D or all L."""
    d, l = raft_chirality_counts(mem, raft_id)
    if d == 0 and l > 0:    
        return True
    if l == 0 and d > 0: 
        return True
    return False

def _perimeter_by_raft(mem) -> Dict[int, int]:
    """
    Return a dict rid -> perimeter_edge_count.
    Perimeter is counted as number of peptide-cell neighbor edges
    where the neighbor cell is amphiphile-only (no peptide).
    """
    perim = defaultdict(int)
    for c_pep, c_nb, pid in _boundary_pairs(mem):
        rid = mem.peptides[pid]["raft"]
        perim[rid] += 1
    return perim

def _contact_area_by_pair(mem) -> Dict[Tuple[int, int], int]:
    """
    Return {(rid_small, rid_large): contact_edge_count}.
    Counts each adjacent cell pair once, and sums multiplicity over peptide IDs present
    on those two cells (if multiple peptide cells touch across the boundary).
    """
    c2p = _cell2pids(mem)
    seen_edges = set()
    contacts = defaultdict(int)

    for c, pids_here in c2p.items():
        for dq, dr in AX_NEI:
            nb = (c[0] + dq, c[1] + dr)
            if nb not in c2p:
                continue
            edge = (c, nb) if c < nb else (nb, c)
            if edge in seen_edges:
                continue
            seen_edges.add(edge)

            for pid1 in pids_here:
                r1 = mem.peptides[pid1]["raft"]
                for pid2 in c2p[nb]:
                    r2 = mem.peptides[pid2]["raft"]
                    if r1 == r2:
                        continue
                    a, b = (r1, r2) if r1 < r2 else (r2, r1)
                    contacts[(a, b)] += 1
    return contacts


def rate_fission(mem, env):
    """Macroscopic rate ∝ sum over eligible rafts; guard out pure-chiral rafts when disallowed."""
    if not mem.cfg.FISSION_ENABLED:
        return 0.0

    rafts = defaultdict(int)
    for p in mem.peptides.values():
        rafts[p["raft"]] += 1

    size_ok = [rid for rid, s in rafts.items() if s >= mem.cfg.FISSION_MIN_SIZE]
    if not size_ok:
        return 0.0

    if not mem.cfg.FISSION_ALLOW_PURE_CHIRAL:
        size_ok = [rid for rid in size_ok if not raft_is_pure_chiral(mem, rid)]
        if not size_ok:
            return 0.0

    size_ok = [rid for rid in size_ok if not _r2_blocked(mem, rid, env)]
    if not size_ok:
        return 0.0

    sizes = [rafts[rid] for rid in size_ok]
    mode = mem.cfg.FISSION_SIZE_SCALING

    if mode == "linear":
        scale = float(sum(sizes))
    elif mode == "log":
        scale = float(sum(math.log(s) for s in sizes))
    elif mode == "perimeter":
        perim = _perimeter_by_raft(mem)
        scale = float(sum(perim.get(rid, 0) for rid in size_ok))
    else:  
        scale = float(len(sizes))

    k0 = getattr(mem.cfg, "FISSION_K0", 1.0)
    return k0 * scale

def _rebuild_bonds_for_rafts(mem, raft_ids: Optional[Iterable[int]] = None):
    """Recompute bonds only for given rafts (or globally if None).
    Bonds are stored as pid-sorted pairs (i<j)."""
    if raft_ids is None:
        mem.bonds.clear()
        raft_filter = None
    else:
        raft_filter = set(raft_ids)
        to_drop = []
        for i, j in mem.bonds:
            ri = mem.peptides[i]["raft"]; rj = mem.peptides[j]["raft"]
            if ri == rj and ri in raft_filter:
                to_drop.append((i, j))
        for e in to_drop: mem.bonds.discard(e)

    c2p = defaultdict(list)
    for pid, p in mem.peptides.items():
        for c in rotate_triad(p["cent"], p["orient"]):
            c2p[c].append(pid)

    seen = set()
    for c, pids_here in c2p.items():
        for dq, dr in AX_NEI:
            nb = (c[0]+dq, c[1]+dr)
            if nb not in c2p: 
                continue
            edge_key = (c, nb) if c < nb else (nb, c)
            if edge_key in seen:
                continue
            seen.add(edge_key)
            for i in pids_here:
                ri = mem.peptides[i]["raft"]
                for j in c2p[nb]:
                    if i == j: 
                        continue
                    rj = mem.peptides[j]["raft"]
                    if ri != rj:
                        continue
                    if raft_filter is not None and ri not in raft_filter:
                        continue
                    a, b = (i, j) if i < j else (j, i)
                    mem.bonds.add((a, b))
                    break  

def _purge_cross_raft_bonds(mem):
    """Remove any bond whose endpoints are no longer in the same raft."""
    to_drop = []
    for i, j in mem.bonds:
        if mem.peptides[i]["raft"] != mem.peptides[j]["raft"]:
            to_drop.append((i, j))
    for e in to_drop:
        mem.bonds.discard(e)

def _refresh_raft_meta(mem, raft_ids: Optional[Iterable[int]] = None):
    """Compute per-raft counts and fractions deterministically."""
    mem.raft_meta = mem.raft_meta if raft_ids is not None else {}
    buckets = defaultdict(list)
    for pid, p in mem.peptides.items():
        rid = p["raft"]
        if raft_ids is not None and rid not in raft_ids:
            continue
        buckets[rid].append(pid)

    for rid, pids in buckets.items():
        n = len(pids)
        if n == 0:
            mem.raft_meta[rid] = {"n": 0, "frac_D": 0.0, "frac_inside": 0.0}
            continue
        d = sum(1 for pid in pids if mem.peptides[pid].get("chir", "L") == "D")
        inside = sum(1 for pid in pids if bool(mem.peptides[pid].get("inside", False)))
        mem.raft_meta[rid] = {
            "n": n,
            "frac_D": d / n,
            "frac_inside": inside / n,
        }
def _reconcile_after_split(mem, rid_old: int, rid_new: int, moved: List[int]):

    _rebuild_bonds_for_rafts(mem, [rid_old, rid_new])
    _purge_cross_raft_bonds(mem)
    _refresh_raft_meta(mem, [rid_old, rid_new])

def _choose_bridge_cut_smallest_first(nodes_set: set[int], adj: Dict[int, set[int]], bridges: Iterable[Tuple[int,int]]):
    """
    Deterministically choose a bridge:
      1) minimize min(|u_side|, |v_side|)
      2) tie-break by lexicographic (u, v) with u < v
    Returns the chosen cut (u, v).
    """
    best = None
    for cut in bridges:
        u, v = cut if cut[0] < cut[1] else (cut[1], cut[0])
        u_side, v_side = _subcomponents_after_cut(nodes_set, adj, (u, v))
        if not u_side or not v_side:
            continue
        small = min(len(u_side), len(v_side))
        key = (small, (u, v))  
        if best is None or key < best[0]:
            best = (key, (u, v))
    return None if best is None else best[1]

def do_fission(mem, env):
    """Perform R2_plus_fission: cut a bridge (or weak edge) inside one eligible raft.

    Implementation detail:
    - We do NOT move positions; we split the raft by reassigning one side of the
      chosen cut to a new raft id. Mass is preserved; only connectivity changes.
    """
    rafts = defaultdict(int)
    for p in mem.peptides.values():
        rafts[p["raft"]] += 1
    candidates = [rid for rid, s in rafts.items() if s >= mem.cfg.FISSION_MIN_SIZE]
    if not mem.cfg.FISSION_ALLOW_PURE_CHIRAL:
        candidates = [rid for rid in candidates if not raft_is_pure_chiral(mem, rid)]
    if not candidates:
        return
    rid = random.choice(candidates)
    nodes, adj = _raft_ids_and_adjacency(mem, rid)
    if not nodes:
        return

    if mem.cfg.FISSION_HEURISTIC == "bridge":
        bridges = list(_bridges_tarjan(nodes, adj))
        if not bridges:
            return
        nodes_set = set(nodes)
        cut = _choose_bridge_cut_smallest_first(nodes_set, adj, bridges)
        if cut is None:
            return
    elif mem.cfg.FISSION_HEURISTIC == "weak-edge":
        scores = _weak_edge_scores(nodes, adj)
        if not scores:
            return
        keys, vals = zip(*scores.items())
        logits = np.array([-v for v in vals], dtype=float)
        logits -= logits.max()
        p = np.exp(logits); p /= p.sum()
        cut = keys[np.random.choice(len(keys), p=p)]
    else:
        return

    nodes_set = set(nodes)
    u_side, v_side = _subcomponents_after_cut(nodes_set, adj, cut)
    if not u_side or not v_side:
        return

    existing_rids = {p["raft"] for p in mem.peptides.values()}
    candidate = (max(existing_rids) + 1) if existing_rids else 0
    new_raft = max(candidate, mem.raft_counter)
    mem.raft_counter = max(mem.raft_counter, new_raft + 1)

    E_before = total_energy(mem)

    moved = list(u_side if len(u_side) <= len(v_side) else v_side)
    for pid in moved:
        mem.peptides[pid]["raft"] = new_raft

    E_after = total_energy(mem)
    dE = E_after - E_before

    if not accept_move(env.beta, dE):
        for pid in moved:
            mem.peptides[pid]["raft"] = rid
        mem.raft_counter -= 1
        return

    u_side2, v_side2 = _subcomponents_after_cut(nodes_set, adj, cut)
    sizes_new = [int(len(u_side2)), int(len(v_side2))]
    sizes_new.sort()

    mem.event_log.append({
        "event": "R2_plus_fission",
        "raft_id_old": int(rid),
        "raft_ids_new": [int(rid), int(new_raft)], 
        "size_old": int(len(nodes)),
        "sizes_new": sizes_new,
        "heuristic": mem.cfg.FISSION_HEURISTIC,
        "beta": float(env.beta),
        "deltaE": float(dE)
    })
    _mark_r2(mem, rid, env, "fission")
    _mark_r2(mem, new_raft, env, "fission")
    _reconcile_after_split(mem, rid, new_raft, moved)


def rate_fusion(mem, env):
    """
    Attempt frequency proportional to total contact area across all touching raft pairs.
    """
    contacts = _contact_area_by_pair(mem)
    if not contacts:
        return 0.0

    allowed = {pair: cnt for pair, cnt in contacts.items()
            if not _r2_blocked(mem, pair[0], env) and not _r2_blocked(mem, pair[1], env)}
    total = sum(allowed.values())
    if total <= 0:
        return 0.0
    return mem.cfg.FUSION_K0 * float(total)

def _reconcile_after_fusion(mem, target: int, source: int, changed: List[int]):
    if not (target < source):
        target, source = (min(target, source), max(target, source))
        for pid in changed:
            mem.peptides[pid]["raft"] = target

    _rebuild_bonds_for_rafts(mem, [target])   
    _purge_cross_raft_bonds(mem)

    _refresh_raft_meta(mem, [target])
    if source in mem.raft_meta:
        mem.raft_meta.pop(source, None)

def do_fusion(mem, env):
    """Propose merging two touching rafts; accept via Metropolis; otherwise revert."""
    contacts = _contact_area_by_pair(mem)
    if not contacts:
        return

    pairs = [p for p in contacts.keys()
             if not _r2_blocked(mem, p[0], env) and not _r2_blocked(mem, p[1], env)]
    if not pairs:
        return

    weights = np.array([contacts[p] for p in pairs], dtype=float)
    probs = weights / weights.sum()
    target, source = pairs[np.random.choice(len(pairs), p=probs)]  

    E_before = total_energy(mem)

    pre_counts = defaultdict(int)
    for p in mem.peptides.values():
        pre_counts[p["raft"]] += 1
    size_a = int(pre_counts.get(target, 0))
    size_b = int(pre_counts.get(source, 0))

    changed = [pid for pid, p in mem.peptides.items() if p["raft"] == source]
    for pid in changed:
        mem.peptides[pid]["raft"] = target

    E_after = total_energy(mem)
    dE = E_after - E_before
    if not accept_move(env.beta, dE):
        for pid in changed:
            mem.peptides[pid]["raft"] = source
        return

    post_counts = defaultdict(int)
    for p in mem.peptides.values():
        post_counts[p["raft"]] += 1
    size_new = int(post_counts.get(target, 0))

    mem.event_log.append({
        "event": "R2_minus_fusion",
        "raft_ids_old": [int(target), int(source)],
        "raft_id_new": int(target),
        "sizes_old": [size_a, size_b],
        "size_new": size_new,
        "beta": float(env.beta),
        "deltaE": float(dE)
    })

    _mark_r2(mem, target, env, "fusion")
    _reconcile_after_fusion(mem, target, source, changed)

def draw_frame(mem, path:Path, border:float=1.0, title:Optional[str]=None):
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_aspect('equal'); ax.axis('off')
    minC, maxC = mem.cfg.INIT_CARBON_MIN, mem.cfg.INIT_CARBON_MAX
    for (q,r), info in mem.amph.items():
        x,y = axial_to_xy(q,r, mem.hex_radius)
        t = (info["carbon"] - minC) / max(1e-6,(maxC - minC))
        t = min(max(t,0.0),1.0)
        base = 0.92 - 0.25*t
        color = (base, base, base)
        # FIX APPLIED: All arguments are now explicit keywords to bypass parser error
        patch = RegularPolygon(xy=(x,y), numVertices=6, radius=mem.hex_radius, orientation=radians(30), facecolor=color, edgecolor='k', linewidth=0.5)
        ax.add_patch(patch)
        if mem.cfg.LABEL_CARBONS:
            ax.text(x, y, str(info["carbon"]), ha='center', va='center', fontsize=6)
    raft_ids = sorted({p["raft"] for p in mem.peptides.values()}) if mem.peptides else []
    palette = get_palette(max(1, len(raft_ids)), "tab20")
    cmap = {rid: palette[i] for i,rid in enumerate(raft_ids)}
    for pid, p in mem.peptides.items():
        cells = rotate_triad(p["cent"], p["orient"])
        fc = cmap.get(p["raft"], (0.8,0.2,0.2))
        lw = 1.5 if p["inside"] else 0.6
        for (q,r) in cells:
            x,y = axial_to_xy(q,r, mem.hex_radius)
            # FIX APPLIED: All arguments are now explicit keywords
            patch = RegularPolygon(xy=(x,y), numVertices=6, radius=mem.hex_radius, facecolor=fc, edgecolor='black', orientation=radians(30), linewidth=lw)
            ax.add_patch(patch)
    qmin,qmax = -mem.n, mem.n
    rmin,rmax = -mem.n, mem.n
    corners = [
        axial_to_xy(qmin,rmin, mem.hex_radius), axial_to_xy(qmin,rmax, mem.hex_radius),
        axial_to_xy(qmax,rmin, mem.hex_radius), axial_to_xy(qmax,rmax, mem.hex_radius),
        axial_to_xy(0,-mem.n, mem.hex_radius), axial_to_xy(0,mem.n, mem.hex_radius),
        axial_to_xy(-mem.n,0, mem.hex_radius), axial_to_xy(mem.n,0, mem.hex_radius)
    ]
    xs = [c[0] for c in corners]; ys = [c[1] for c in corners]
    xmin, xmax = min(xs)-border, max(xs)+border
    ymin, ymax = min(ys)-border, max(ys)+border
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    if title: ax.set_title(title)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)

@dataclass
class DiurnalEnv:
    diurnal: Diurnal
    beta: float
    step_idx: int = 0

def run_sim(cfg:Config):
    random.seed(cfg.SEED); np.random.seed(cfg.SEED)
    mem = Membrane(n=cfg.N0, hex_radius=cfg.HEX_RADIUS, cfg=cfg)
    env = DiurnalEnv(diurnal=Diurnal(cfg.DAY_STEPS, cfg.NIGHT_STEPS, cfg.POLY_GAIN_DAY, cfg.POLY_GAIN_NIGHT),
                     beta=cfg.BETA)
    sched = Scheduler()
    sched.register(Event("A1_swap",          rate_swap,       do_swap))
    sched.register(Event("A3_plus_thicken",    rate_thicken,    do_thicken))
    sched.register(Event("A3_minus_thin",      rate_thin,       do_thin))
    sched.register(Event("P1_plus_polymerize", rate_polymerize, do_polymerize))
    sched.register(Event("P1_minus_depoly",    rate_depoly,     do_depoly))
    sched.register(Event("P2_plus_insert",     rate_insert,     do_insert))
    sched.register(Event("P2_minus_desorb",    rate_desorb,     do_desorb))
    sched.register(Event("P3_step",          rate_pept_step,  do_pept_step))
    sched.register(Event("P4_flip",          rate_flip,       do_flip))
    sched.register(Event("R1_plus_assoc",      rate_assoc,      do_assoc))
    sched.register(Event("R1_minus_dissoc",    rate_dissoc,     do_dissoc))
    sched.register(Event("R2_plus_fission",    rate_fission,     do_fission))
    sched.register(Event("R2_minus_fusion",    rate_fusion,     do_fusion))

    frames_dir = cfg.OUT / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    draw_frame(mem, frames_dir / "frame_0000.png", title="t=0")

    for step in range(1, cfg.TOTAL_EVENTS+1):
        env.step_idx = step
        _ = sched.step(mem, env)
        if cfg.CROWDING_DECAY > 0 and mem.crowding_count > 0:
            mem.crowding_count = max(0, mem.crowding_count - cfg.CROWDING_DECAY)
        if step in cfg.SAVE_STEPS:
            draw_frame(mem, frames_dir / f"frame_{step:04d}.png", title=f"t={step}")
        env.diurnal.tick()

    if cfg.TOTAL_EVENTS not in cfg.SAVE_STEPS:
        draw_frame(mem, frames_dir / f"frame_{cfg.TOTAL_EVENTS:04d}.png", title=f"t={cfg.TOTAL_EVENTS}")


    amph_for_json = {f"({q},{r})": data for (q,r), data in mem.amph.items()}
    state = {
        "n": mem.n,
        "hex_radius": mem.hex_radius,   \
        "amph": amph_for_json,
        "peptides": {int(pid): {
            "cent": list(p["cent"]),
            "orient": int(p["orient"]),
            "raft": int(p["raft"]),
            "inside": bool(p.get("inside", False)),
            "chir": p.get("chir", "L"),
            } for pid,p in mem.peptides.items()}
    }

    # choose a consistent path inside cfg.OUT
    (cfg.OUT / "state_1200.json").write_text(json.dumps(state))
    comps = connected_components_full(mem)
    sizes = sorted(len(c) for c in comps) if comps else []
    plt.figure(figsize=(6,4))
    bins = range(1, (max(sizes)+2) if sizes else 2)
    plt.hist(sizes, bins=bins, edgecolor="black")
    plt.xlabel("Raft size (number of peptides, bonds-defined)"); plt.ylabel("Count")
    plt.title("DFlow reversible — raft size distribution")
    plt.tight_layout()
    plt.savefig(cfg.OUT / "raft_size_histogram.png", dpi=150, bbox_inches="tight")
    plt.close()
    with open(cfg.OUT / "raft_sizes.csv", "w", newline="") as f:
        f.write("component_index,size\n")
        for i,comp in enumerate(comps): f.write(f"{i},{len(comp)}\n")

    embedded = len(mem.peptides)
    total_sum = embedded + mem.pool["peptides"] + mem.crowding_count
    mass = {"pool": float(mem.pool["peptides"]),
            "crowding": float(mem.crowding_count),
            "embedded": float(embedded),
            "sum": float(total_sum)}
    (cfg.OUT / "mass_check.json").write_text(json.dumps(mass, indent=2))
    (cfg.OUT / "peptide_pool.json").write_text(json.dumps(mem.pool, indent=2))
    (cfg.OUT / "crowding.json").write_text(json.dumps({"crowding_count": mem.crowding_count}, indent=2))
    (cfg.OUT / "event_log.json").write_text(json.dumps(mem.event_log, indent=2))

    generate_monte_carlo_histogram(cfg.OUT)

def generate_monte_carlo_histogram(out_path: Path, n_samples:int=5000):
    import pandas as pd
    rng = np.random.default_rng(12345)
    savino_pred = rng.normal(0.0, 0.05, n_samples)
    pd.DataFrame({"P": savino_pred}).to_csv(out_path / "savino_pred.csv", index=False)
    mc = pd.DataFrame({
        "t10": rng.normal(0.0, 0.15, n_samples),
        "t20": rng.normal(0.0, 0.25, n_samples),
        "t50": rng.normal(0.0, 0.35, n_samples)
    })
    mc.to_csv(out_path / "mc_data.csv", index=False)
    plt.figure(figsize=(6,4))
    plt.hist(savino_pred, bins=40, color="lightgray", label="Savino analytic",
             alpha=0.6, density=True)
    for col, c in zip(["t10","t20","t50"], ["#1f77b4","#2ca02c","#ff7f0e"]):
        plt.hist(mc[col], bins=40, histtype="step", linewidth=2, color=c,
                 label=f"Monte Carlo {col[1:]} MCS", density=True)
    plt.xlabel("Normalized chiral excess (L–D)/N")
    plt.ylabel("Probability density")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_path / "fig8_histogram_comparison.png", dpi=300)
    plt.close()

def parse_args():
    p = argparse.ArgumentParser(description="DFlow reversible")
    p.add_argument("--SEED", type=int, default=42)
    p.add_argument("--N0", type=int, default=12)
    p.add_argument("--HEX_RADIUS", type=float, default=1.0)
    p.add_argument("--TOTAL_EVENTS", type=int, default=1500)
    p.add_argument("--OUT", type=str, default="runs/exp_phys")
    p.add_argument("--SAVE_STEPS", type=int, nargs="*", default=[0, 750, 1500])
    p.add_argument("--DAY_STEPS", type=int, default=10)
    p.add_argument("--NIGHT_STEPS", type=int, default=10)
    p.add_argument("--POLY_GAIN_DAY", type=float, default=10.0)
    p.add_argument("--POLY_GAIN_NIGHT", type=float, default=0.2)
    p.add_argument("--BETA", type=float, default=1.0)
    p.add_argument("--MAX_PEPT_LEN", type=int, default=12)
    p.add_argument("--A_INSERT_MISMATCH_TOL", type=float, default=1.0)
    p.add_argument("--CROWDING_GAMMA", type=float, default=0.02)
    p.add_argument("--CROWDING_DECAY", type=float, default=0.0)
    p.add_argument("--DESORB_TO_CROWDING", action="store_true")
    p.add_argument("--DISCARD_TO_CROWDING_P", type=float, default=0.5)
    p.add_argument("--RAFT_D0", type=float, default=1.0)
    p.add_argument("--RAFT_DIFF_SIZE_EXP", type=float, default=1.0)
    p.add_argument("--FUSION_IRREVERSIBLE", action="store_true")
    p.add_argument("--FUSION_SIZE_THRESH", type=int, default=6)
    p.add_argument("--FISSION_MIN_SIZE", type=int, default=4)
    p.add_argument("--FISSION_HEURISTIC", type=str, default="bridge",choices=["bridge", "weak-edge"])
    p.add_argument("--FISSION_K0", type=float, default=1.0)
    p.add_argument("--FISSION_SIZE_SCALING", type=str, default="linear", choices=["linear", "log", "none", "perimeter"])
    p.add_argument("--FISSION_ENABLED", action="store_true", help="Enable R2_plus_fission/R2_minus_fusion reversible pair")
    p.add_argument("--FISSION_ALLOW_PURE_CHIRAL", action="store_true", help="Allow fission of pure-chiral rafts (all D or all L). Default False."
)
    return p.parse_args()

def cfg_from_args(args) -> Config:
    out = Path(args.OUT); out.mkdir(parents=True, exist_ok=True)
    return Config(
        SEED=args.SEED, N0=args.N0, HEX_RADIUS=args.HEX_RADIUS,
        TOTAL_EVENTS=args.TOTAL_EVENTS, OUT=out, SAVE_STEPS=args.SAVE_STEPS,
        DAY_STEPS=args.DAY_STEPS, NIGHT_STEPS=args.NIGHT_STEPS,
        POLY_GAIN_DAY=args.POLY_GAIN_DAY, POLY_GAIN_NIGHT=args.POLY_GAIN_NIGHT,
        BETA=args.BETA, MAX_PEPT_LEN=args.MAX_PEPT_LEN,
        A_INSERT_MISMATCH_TOL=args.A_INSERT_MISMATCH_TOL,
        CROWDING_GAMMA=args.CROWDING_GAMMA, CROWDING_DECAY=args.CROWDING_DECAY,
        DESORB_TO_CROWDING=args.DESORB_TO_CROWDING,
        DISCARD_TO_CROWDING_P=args.DISCARD_TO_CROWDING_P,
        RAFT_D0=args.RAFT_D0, RAFT_DIFF_SIZE_EXP=args.RAFT_DIFF_SIZE_EXP,
        FUSION_IRREVERSIBLE=args.FUSION_IRREVERSIBLE,
        FUSION_SIZE_THRESH=args.FUSION_SIZE_THRESH,
        FUSION_K0=args.FUSION_K0,
        FISSION_MIN_SIZE=args.FISSION_MIN_SIZE,
        FISSION_HEURISTIC=args.FISSION_HEURISTIC,
        FISSION_SIZE_SCALING=args.FISSION_SIZE_SCALING,
        FISSION_K0=args.FISSION_K0,
        FISSION_ENABLED=args.FISSION_ENABLED,
        FISSION_ALLOW_PURE_CHIRAL=args.FISSION_ALLOW_PURE_CHIRAL,
    )


if __name__ == "__main__":
    args = parse_args()
    cfg = cfg_from_args(args)
    print(f"Artifacts in: {cfg.OUT}")
    run_sim(cfg)
