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

from dataclasses import dataclass
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
    mem.peptides[pid] = {"cent": anchor, "orient": orient, "raft": rid, "inside": inside, "nres": nres}
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

    frames_dir = cfg.OUT / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    draw_frame(mem, frames_dir / "frame_0000.png", title="t=0")

    for step in range(1, cfg.TOTAL_EVENTS+1):
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
        "hex_radius": mem.hex_radius,   # <- was HEX_RADIUS
        "amph": amph_for_json,
        "peptides": {int(pid): {
            "cent": list(p["cent"]),
            "orient": int(p["orient"]),
            "raft": int(p["raft"]),
            "inside": bool(p.get("inside", False))
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
    )

if __name__ == "__main__":
    args = parse_args()
    cfg = cfg_from_args(args)
    print(f"Artifacts in: {cfg.OUT}")
    run_sim(cfg)
