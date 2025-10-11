#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DFlow (physically-grounded, reversible) — with mismatch window + soft acceptance,
peptide length cap (~12 aa), crowding effect on polymerization, size-dependent raft diffusion,
and optional quasi-irreversible fusion for large rafts (off by default).

Usage:
  python dflow_reversible.py --N0 12 --TOTAL_EVENTS 1500 --OUT runs/exp_phys --SAVE_STEPS 0 750 1500
"""
import os, json, math, random, argparse
from pathlib import Path
from collections import defaultdict, deque
from typing import Optional, List

import numpy as np
import matplotlib as mpl
mpl.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

# -----------------------------
# Physical knobs / parameters
# -----------------------------
# Peptide constraints & mismatch gating
MAX_PEPT_LEN = 12            # Deamer/Damer note; hard cap
EPS_MISMATCH_NM = 0.5        # hard window for |Lp - d_site| <= eps (strict gate)
BETA = 1.0                   # inverse temperature for soft acceptance (Boltzmann)
SOFT_GATE_FRACTION = 0.6     # portion of EPS used as "free" zone; beyond it costs energy

# Crowding (rejected peptides -> cytoplasmic crowding that slows polymerization)
CROWDING_GAMMA = 0.02        # slowdown coefficient: higher => more slowdown
CROWDING_DECAY = 0.0         # per-step decay of crowding (0 = no decay)

# Route desorbed peptides: False=return to pool (default), True=send to crowding
DESORB_TO_CROWDING = False

# Raft diffusion (2D Stokes-like proxy): D(size) ~ D0 / size^exp
RAFT_D0 = 1.0
RAFT_DIFF_SIZE_EXP = 1.0

# Optional quasi-irreversible fusion for large rafts (OFF by default)
FUSION_IRREVERSIBLE = False
FUSION_SIZE_THRESH = 6       # rafts >= this size keep internal bonds (no dissoc)

# -----------------------------
# Geometry (flat-top hex, axial coords)
# -----------------------------
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

def rotate_triad(anchor, orient:int):
    aq, ar = anchor
    offs = TRIAD_SHAPES[orient % 6]
    return [(aq+dq, ar+dr) for (dq,dr) in offs]

def get_palette(n:int, name:str="tab20"):
    try:    cmap = mpl.colormaps[name]
    except: cmap = plt.get_cmap(name)
    if hasattr(cmap, "colors"):
        base = list(cmap.colors)
        if n <= len(base): return base[:n]
        reps = (n + len(base) - 1)//len(base)
        return (base*reps)[:n]
    if n == 1: return [cmap(0.0)]
    return [cmap(i/(n-1)) for i in range(n)]

# -----------------------------
# Diurnal driver (bias polymerization)
# -----------------------------
class Diurnal:
    def __init__(self, day_steps:int=10, night_steps:int=10,
                 poly_gain_day:float=10.0, poly_gain_night:float=0.2):
        self.day_steps = max(1, day_steps)
        self.night_steps = max(1, night_steps)
        self.poly_gain_day = float(poly_gain_day)
        self.poly_gain_night = float(poly_gain_night)
        self.t = 0

    @property
    def is_day(self)->bool:
        return (self.t % (self.day_steps + self.night_steps)) < self.day_steps

    def polymerization_gain(self)->float:
        return self.poly_gain_day if self.is_day else self.poly_gain_night

    def tick(self):
        self.t += 1

# -----------------------------
# Membrane & peptide pool state
# -----------------------------
class Membrane:
    def __init__(self, n:int=12, hex_radius:float=1.0, label_carbons:bool=False):
        self.n = n
        self.hex_radius = hex_radius
        self.label_carbons = label_carbons

        # Amphiphiles: mono_di in {1,2} gives ~2nm or ~4nm nominal thickness
        self.amph = {}        # (q,r) -> {"carbon":int, "mono_di":int, "pep":bool}

        # Peptides embedded in membrane
        self.peptides = {}    # pid -> {"cent":(q,r), "orient":int, "inside":bool,
                              #        "chir":'L'|'D'|'0', "length":int, "Lp_nm":float, "helical":bool}

        # Bonds between peptides (reversible raft association)
        self.bonds = set()    # set of frozenset({pid1,pid2})

        # Solution peptide pool (counts)
        self.pool = {"L": 0, "D": 0, "0": 0}

        # Cytoplasmic/medium crowding proxy (rejects accumulate here)
        self.crowding_count = 0

        self.pid_counter = 0
        self.event_log = []
        self._init_membrane()

    def _axial_in_rhombus(self, q, r):
        n = self.n
        return (abs(q) <= n) and (abs(r) <= n) and (abs(q+r) <= n)

    def _init_membrane(self):
        for q in range(-self.n, self.n+1):
            for r in range(-self.n, self.n+1):
                if self._axial_in_rhombus(q,r):
                    self.amph[(q,r)] = {
                        "carbon": random.randint(10,20),
                        "mono_di": 1,      # start as 2 nm
                        "pep": False
                    }

    # --- queries ---
    def triad_cells(self, center, orient):
        return rotate_triad(center, orient)

    def peptide_cells(self, pid):
        p = self.peptides[pid]
        return self.triad_cells(p["cent"], p["orient"])

    def cells_free_for_triad(self, cells):
        for c in cells:
            if c not in self.amph: return False
            if self.amph[c]["pep"]: return False
        return True

    def random_free_anchor(self, max_tries=1000):
        keys = list(self.amph.keys())
        for _ in range(max_tries):
            anchor = random.choice(keys)
            orient = random.randrange(6)
            cells = self.triad_cells(anchor, orient)
            if self.cells_free_for_triad(cells):
                return anchor, orient, cells
        return None, None, None

    # --- peptide pool helpers ---
    def pool_add(self, chir, k=1):
        self.pool[chir] = max(0, self.pool[chir] + int(k))

    def pool_take(self, chir)->bool:
        if self.pool[chir] > 0:
            self.pool[chir] -= 1
            return True
        return False

# -----------------------------
# Energetics / gating
# -----------------------------
def site_thickness_nm(amph_site)->float:
    return 2.0 if amph_site["mono_di"] == 1 else 4.0

def peptide_Lp_nm_from_length(length:int, is_helical:bool)->float:
    per_res_nm_helix = 0.28  # helix rise per residue (nm), tune as needed
    per_res_nm_coil  = 0.15  # coil hydrophobic span proxy
    return (per_res_nm_helix if is_helical else per_res_nm_coil) * float(length)

def allowed_insert_hard(mem: Membrane, pdict, cells, eps_nm: float = EPS_MISMATCH_NM) -> bool:
    """Hard window + require helix."""
    if not pdict.get("helical", False):
        return False
    Lp = pdict["Lp_nm"]
    for c in cells:
        d = site_thickness_nm(mem.amph[c])
        if abs(Lp - d) > eps_nm:
            return False
    return True

def insertion_soft_accept(mem: Membrane, pdict, cells, beta: float = BETA) -> bool:
    """Soft acceptance: small over-mismatch allowed with Boltzmann-weighted probability."""
    Lp = pdict["Lp_nm"]; dE = 0.0
    core_eps = EPS_MISMATCH_NM * SOFT_GATE_FRACTION
    for c in cells:
        d = site_thickness_nm(mem.amph[c])
        over = max(0.0, abs(Lp - d) - core_eps)
        dE += over
    p = math.exp(-beta * dE)
    return random.random() < p

# -----------------------------
# Event framework
# -----------------------------
class Event:
    def __init__(self, name, rate_fn, do_fn):
        self.name = name
        self.rate_fn = rate_fn   # (mem, env) -> float
        self.do_fn = do_fn       # (mem, env) -> None

class DiurnalEnv:
    def __init__(self, diurnal:Diurnal, beta:float=1.0):
        self.diurnal = diurnal
        self.beta = beta
        self.k0_poly = 1.0
        self.k0 = 1.0

class Scheduler:
    """Weighted-discrete chooser (approx. Gillespie): normalize rates -> pick one."""
    def __init__(self):
        self.events: List[Event] = []

    def register(self, ev:Event):
        self.events.append(ev)

    def step(self, mem:Membrane, env:DiurnalEnv):
        rates = [max(0.0, float(ev.rate_fn(mem, env))) for ev in self.events]
        total = sum(rates)
        if total <= 0:
            return None
        probs = [r/total for r in rates]
        idx = np.random.choice(len(self.events), p=probs)
        self.events[idx].do_fn(mem, env)
        return self.events[idx].name

# -----------------------------
# Micro-events & rates
# -----------------------------
# A1: lateral swap (self-inverse)
def rate_swap(mem:Membrane, env:DiurnalEnv):
    return env.k0 * len(mem.amph)

def do_swap(mem:Membrane, env:DiurnalEnv):
    cells = list(mem.amph.keys())
    if not cells: return
    a = random.choice(cells)
    q,r = a
    random.shuffle(AX_NEI)
    for dq,dr in AX_NEI:
        b = (q+dq, r+dr)
        if b in mem.amph:
            mem.amph[a]["carbon"], mem.amph[b]["carbon"] = mem.amph[b]["carbon"], mem.amph[a]["carbon"]
            mem.event_log.append({"evt":"A1_swap","a":a,"b":b})
            break

# A3+/A3−: thickness mono<->di (reversible)
def count_thickenable(mem:Membrane):
    return sum(1 for v in mem.amph.values() if v["mono_di"]==1)

def count_thinnable(mem:Membrane):
    return sum(1 for v in mem.amph.values() if v["mono_di"]==2)

def rate_thicken(mem:Membrane, env:DiurnalEnv):
    return env.k0 * count_thickenable(mem)

def rate_thin(mem:Membrane, env:DiurnalEnv):
    return env.k0 * count_thinnable(mem)

def do_thicken(mem:Membrane, env:DiurnalEnv):
    monos = [c for c,v in mem.amph.items() if v["mono_di"]==1]
    if not monos: return
    c = random.choice(monos)
    mem.amph[c]["mono_di"] = 2
    mem.event_log.append({"evt":"A3_plus_thicken","cell":c})

def do_thin(mem:Membrane, env:DiurnalEnv):
    dis = [c for c,v in mem.amph.items() if v["mono_di"]==2]
    if not dis: return
    c = random.choice(dis)
    mem.amph[c]["mono_di"] = 1
    mem.event_log.append({"evt":"A3_minus_thin","cell":c})

# P1+/P1−: polymerize/depolymerize (BIAS via diurnal) + crowding slowdown
def rate_polymerize(mem:Membrane, env:DiurnalEnv):
    base = env.k0_poly * env.diurnal.polymerization_gain()
    C = getattr(mem, "crowding_count", 0)
    f = 1.0 / (1.0 + CROWDING_GAMMA * C)  # PCR-like slowdown with crowding
    return base * f

def rate_depoly(mem:Membrane, env:DiurnalEnv):
    pool_total = sum(mem.pool.values())
    return env.k0_poly * (0.2 + 0.8 * (pool_total > 0))

def do_polymerize(mem:Membrane, env:DiurnalEnv):
    chir = random.choices(["L","D","0"], weights=[0.48,0.48,0.04])[0]
    mem.pool_add(chir, k=1)
    mem.event_log.append({"evt":"P1_plus_polymerize","chir":chir})

def do_depoly(mem:Membrane, env:DiurnalEnv):
    for ch in ["0","L","D"]:
        if mem.pool[ch] > 0:
            mem.pool[ch] -= 1
            mem.event_log.append({"evt":"P1_minus_depolymerize","chir":ch})
            return

# P2+/P2−: insert / desorb (mismatch-gated)
def new_peptide_from_pool(mem: Membrane):
    choices = [ch for ch, cnt in mem.pool.items() if cnt > 0]
    if not choices:
        return None
    chir = random.choice(choices)
    mem.pool_take(chir)

    # Length: biased to <= MAX_PEPT_LEN; hard cap
    raw_len = random.randint(5, MAX_PEPT_LEN + 2)
    length = min(raw_len, MAX_PEPT_LEN)

    # Secondary structure: helical vs coil
    is_helical = (random.random() < 0.7)  # tweak as needed
    Lp = peptide_Lp_nm_from_length(length, is_helical)

    return {"chir":chir, "length":length, "Lp_nm":Lp, "helical":is_helical}

def rate_insert(mem:Membrane, env:DiurnalEnv):
    if sum(mem.pool.values()) == 0: return 0.0
    return env.k0 * float(len(mem.amph))

def rate_desorb(mem:Membrane, env:DiurnalEnv):
    return env.k0 * float(len(mem.peptides)) if mem.peptides else 0.0

def do_insert(mem:Membrane, env:DiurnalEnv):
    if sum(mem.pool.values()) == 0: return
    anchor, orient, cells = mem.random_free_anchor()
    if anchor is None: return
    proto = new_peptide_from_pool(mem)
    if proto is None: return

    ok_hard = allowed_insert_hard(mem, proto, cells)
    ok_soft = insertion_soft_accept(mem, proto, cells)

    if not (ok_hard or ok_soft):
        # peptide does not fit membrane => crowding
        mem.crowding_count += 1
        mem.event_log.append({
            "evt": "P2_reject_to_crowding",
            "chir": proto["chir"], "helical": proto["helical"], "length": proto["length"]
        })
        return

    pid = mem.pid_counter; mem.pid_counter += 1
    mem.peptides[pid] = {
        "cent": anchor, "orient": orient, "inside": bool(random.getrandbits(1)),
        "chir": proto["chir"], "length": proto["length"],
        "Lp_nm": proto["Lp_nm"], "helical": proto["helical"]
    }
    for c in cells: mem.amph[c]["pep"] = True
    mem.event_log.append({
        "evt":"P2_plus_insert","pid":pid,"cells":cells,
        "chir":proto["chir"],"length":proto["length"],"helical":proto["helical"]
    })

def do_desorb(mem: Membrane, env: DiurnalEnv):
    """Remove a random embedded peptide and conserve mass by routing it either
    back to the solution pool (reversible counterpart of insertion) or to a
    'crowding' bucket (experimental hypothesis).

    Logs:
      {"evt":"P2_minus_desorb","pid":..., "cells":[...], "chir":"L|D|0", "dest":"pool|crowding"}
    """
    if not mem.peptides:
        return

    pid = random.choice(list(mem.peptides.keys()))
    cells = mem.peptide_cells(pid)
    chir  = mem.peptides[pid]["chir"]

    # free cells
    for c in cells:
        if c in mem.amph:
            mem.amph[c]["pep"] = False

    # remove incident bonds
    mem.bonds = {b for b in mem.bonds if pid not in b}

    # remove peptide
    del mem.peptides[pid]

    # destination: pool (default) or crowding (toggle)
    if DESORB_TO_CROWDING:
        mem.crowding_count += 1
        dest = "crowding"
    else:
        mem.pool_add(chir, 1)
        dest = "pool"

    mem.event_log.append({
        "evt": "P2_minus_desorb",
        "pid": pid,
        "cells": cells,
        "chir": chir,
        "dest": dest
    })

# P3: diffusion (size-dependent via raft components)
def connected_components(mem:Membrane):
    g = defaultdict(list)
    for b in mem.bonds:
        a,bp = list(b)
        g[a].append(bp); g[bp].append(a)
    seen=set(); comps=[]
    for pid in mem.peptides.keys():
        if pid in seen: continue
        comp=set(); q=deque([pid])
        while q:
            u=q.popleft()
            if u in seen: continue
            seen.add(u); comp.add(u)
            for v in g.get(u, []):
                if v not in seen: q.append(v)
        comps.append(comp)
    return comps

def raft_components_and_sizes(mem: Membrane):
    comps = connected_components(mem)
    sizes = [len(c) for c in comps]
    # include singletons (peptides with no bonds) if not captured
    captured = set().union(*comps) if comps else set()
    singletons = [ {pid} for pid in mem.peptides.keys() if pid not in captured ]
    for s in singletons: comps.append(s)
    sizes = [len(c) for c in comps]  # recompute including singletons
    return comps, sizes

def rate_pept_step(mem:Membrane, env:DiurnalEnv):
    if not mem.peptides:
        return 0.0
    comps, sizes = raft_components_and_sizes(mem)
    total = 0.0
    for s in (sizes or []):
        s = max(1, s)
        D = RAFT_D0 / (s ** RAFT_DIFF_SIZE_EXP)
        total += D
    return env.k0 * total

def do_pept_step(mem:Membrane, env:DiurnalEnv):
    if not mem.peptides: return
    comps, sizes = raft_components_and_sizes(mem)
    if not comps: return
    # pick a raft weighted by its diffusion coefficient
    weights = []
    for s in sizes:
        s = max(1, s)
        weights.append(RAFT_D0 / (s ** RAFT_DIFF_SIZE_EXP))
    total = sum(weights)
    if total <= 0: return
    probs = [w/total for w in weights]
    idx = np.random.choice(len(comps), p=probs)
    comp = list(comps[idx])

    # move one random member (proxy for raft Brownian step)
    pid = random.choice(comp)
    p = mem.peptides[pid]
    q,r = p["cent"]
    dq,dr = random.choice(AX_NEI)
    new_cent = (q+dq, r+dr)
    new_cells = mem.triad_cells(new_cent, p["orient"])
    if any((c not in mem.amph or mem.amph[c]["pep"]) for c in new_cells):
        return
    # Check mismatch for THIS peptide at destination
    pdict = {"Lp_nm": p["Lp_nm"], "helical": p.get("helical", True)}
    if not (allowed_insert_hard(mem, pdict, new_cells) or insertion_soft_accept(mem, pdict, new_cells)):
        return

    # free old cells & occupy new
    for c in mem.triad_cells(p["cent"], p["orient"]):
        if c in mem.amph: mem.amph[c]["pep"] = False
    for c in new_cells: mem.amph[c]["pep"] = True
    p["cent"] = new_cent
    mem.event_log.append({"evt":"P3_step","pid":pid,"from":(q,r),"to":new_cent,"raft_size":len(comp)})

# P4: orientation flip (metadata)
def rate_flip(mem:Membrane, env:DiurnalEnv):
    return env.k0 * (len(mem.peptides) if mem.peptides else 0.0)

def do_flip(mem:Membrane, env:DiurnalEnv):
    if not mem.peptides: return
    pid = random.choice(list(mem.peptides.keys()))
    mem.peptides[pid]["inside"] = not mem.peptides[pid]["inside"]
    mem.event_log.append({"evt":"P4_flip","pid":pid,"inside":mem.peptides[pid]["inside"]})

# R1+/R1−: association/dissociation bonds
def adjacent_peptide_pairs(mem:Membrane):
    cell_to_pid = {}
    for pid, p in mem.peptides.items():
        for c in mem.triad_cells(p["cent"], p["orient"]):
            cell_to_pid[c] = pid
    pairs = set()
    for c, pid in cell_to_pid.items():
        q,r = c
        for dq,dr in AX_NEI:
            nb = (q+dq, r+dr)
            if nb in cell_to_pid:
                pid2 = cell_to_pid[nb]
                if pid2 != pid:
                    pairs.add(tuple(sorted((pid,pid2))))
    return pairs

def rate_assoc(mem:Membrane, env:DiurnalEnv):
    adj = adjacent_peptide_pairs(mem)
    candidate = [frozenset(p) for p in adj if frozenset(p) not in mem.bonds]
    return env.k0 * len(candidate)

def rate_dissoc(mem:Membrane, env:DiurnalEnv):
    if not mem.bonds:
        return 0.0
    if FUSION_IRREVERSIBLE:
        # suppress dissociation of bonds fully within large rafts
        comps = connected_components(mem)
        large_pids = set().union(*(c for c in comps if len(c) >= FUSION_SIZE_THRESH)) if comps else set()
        bonds_kept = [b for b in mem.bonds if not (set(b) <= large_pids)]
        return env.k0 * len(bonds_kept)
    return env.k0 * len(mem.bonds)

def do_assoc(mem:Membrane, env:DiurnalEnv):
    adj = adjacent_peptide_pairs(mem)
    candidate = [frozenset(p) for p in adj if frozenset(p) not in mem.bonds]
    if not candidate: return
    b = random.choice(candidate)
    mem.bonds.add(b)
    mem.event_log.append({"evt":"R1_plus_associate","bond":sorted(list(b))})

def do_dissoc(mem:Membrane, env:DiurnalEnv):
    if not mem.bonds: return
    if FUSION_IRREVERSIBLE:
        comps = connected_components(mem)
        large_pids = set().union(*(c for c in comps if len(c) >= FUSION_SIZE_THRESH)) if comps else set()
        # choose only bonds not fully inside a large raft
        candidates = [b for b in mem.bonds if not (set(b) <= large_pids)]
        if not candidates: return
        b = random.choice(candidates)
    else:
        b = random.choice(list(mem.bonds))
    mem.bonds.remove(b)
    mem.event_log.append({"evt":"R1_minus_dissociate","bond":sorted(list(b))})

# -----------------------------
# Visualization
# -----------------------------
def connected_components_full(mem:Membrane):
    """Connected components over the bond graph including singletons."""
    comps = connected_components(mem)
    captured = set().union(*comps) if comps else set()
    singles = [{pid} for pid in mem.peptides.keys() if pid not in captured]
    comps.extend(singles)
    return comps

def draw_frame(mem:Membrane, path:Path, border:float=1.0, title:Optional[str]=None):
    fig, ax = plt.subplots(figsize=(8,8))
    ax.set_aspect('equal'); ax.axis('off')
    minC, maxC = 8.0, 22.0

    # Amphiphiles (gray by carbon; thicker edge for di-sites)
    for (q,r), info in mem.amph.items():
        x,y = axial_to_xy(q,r, mem.hex_radius)
        t = (info["carbon"] - minC) / (maxC - minC)
        t = min(max(t, 0.0), 1.0)
        base = 0.92 - 0.25*t
        lw = 0.6 if info["mono_di"]==1 else 1.5
        patch = RegularPolygon(
            (x,y), numVertices=6, radius=mem.hex_radius, orientation=math.radians(30),
            facecolor=(base,base,base), edgecolor='k', linewidth=lw
        )
        ax.add_patch(patch)
        if mem.label_carbons:
            ax.text(x, y, str(info["carbon"]), ha='center', va='center', fontsize=6)

    # Peptides colored by component index
    comps = connected_components_full(mem)
    pid_to_comp = {}
    for i, comp in enumerate(comps):
        for pid in comp:
            pid_to_comp[pid] = i
    palette = get_palette(max(1, len(comps)))
    for pid, p in mem.peptides.items():
        cells = mem.triad_cells(p["cent"], p["orient"])
        fc = palette[pid_to_comp[pid] % len(palette)]
        lw = 1.5 if p.get("inside", False) else 0.6
        for (q,r) in cells:
            x,y = axial_to_xy(q,r, mem.hex_radius)
            patch = RegularPolygon(
                (x,y), numVertices=6, radius=mem.hex_radius, orientation=math.radians(30),
                facecolor=fc, edgecolor='black', linewidth=lw
            )
            ax.add_patch(patch)

    # limits
    qmin,qmax = -mem.n, mem.n
    rmin,rmax = -mem.n, mem.n
    corners = [
        axial_to_xy(qmin, rmin, mem.hex_radius),
        axial_to_xy(qmin, rmax, mem.hex_radius),
        axial_to_xy(qmax, rmin, mem.hex_radius),
        axial_to_xy(qmax, rmax, mem.hex_radius),
        axial_to_xy(0, -mem.n, mem.hex_radius),
        axial_to_xy(0, mem.n, mem.hex_radius),
        axial_to_xy(-mem.n, 0, mem.hex_radius),
        axial_to_xy(mem.n, 0, mem.hex_radius)
    ]
    xs = [c[0] for c in corners]; ys = [c[1] for c in corners]
    xmin, xmax = min(xs)-border, max(xs)+border
    ymin, ymax = min(ys)-border, max(ys)+border
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    if title: ax.set_title(title)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)

# -----------------------------
# Harness
# -----------------------------
def run_sim(N0:int, TOTAL_EVENTS:int, OUT:Path, SAVE_STEPS:List[int],
            SEED:int=42, HEX_RADIUS:float=1.0,
            DAY_STEPS:int=10, NIGHT_STEPS:int=10,
            POLY_GAIN_DAY:float=10.0, POLY_GAIN_NIGHT:float=0.2,
            BETA_arg:float=1.0):
    random.seed(SEED); np.random.seed(SEED)

    mem = Membrane(n=N0, hex_radius=HEX_RADIUS)
    env = DiurnalEnv(diurnal=Diurnal(DAY_STEPS, NIGHT_STEPS, POLY_GAIN_DAY, POLY_GAIN_NIGHT),
                     beta=BETA_arg)

    sched = Scheduler()
    # reversible pairs + symmetric moves
    sched.register(Event("A1_swap",            lambda m,e: rate_swap(m,e),            do_swap))
    sched.register(Event("A3_plus_thicken",    lambda m,e: rate_thicken(m,e),         do_thicken))
    sched.register(Event("A3_minus_thin",      lambda m,e: rate_thin(m,e),            do_thin))
    sched.register(Event("P1_plus_polymerize", lambda m,e: rate_polymerize(m,e),      do_polymerize))
    sched.register(Event("P1_minus_depoly",    lambda m,e: rate_depoly(m,e),          do_depoly))
    sched.register(Event("P2_plus_insert",     lambda m,e: rate_insert(m,e),          do_insert))
    sched.register(Event("P2_minus_desorb",    lambda m,e: rate_desorb(m,e),          do_desorb))
    sched.register(Event("P3_step",            lambda m,e: rate_pept_step(m,e),       do_pept_step))
    sched.register(Event("P4_flip",            lambda m,e: rate_flip(m,e),            do_flip))
    sched.register(Event("R1_plus_assoc",      lambda m,e: rate_assoc(m,e),           do_assoc))
    sched.register(Event("R1_minus_dissoc",    lambda m,e: rate_dissoc(m,e),          do_dissoc))

    frames_dir = OUT / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    draw_frame(mem, frames_dir / "frame_0000.png", title="t=0")

    for step in range(1, TOTAL_EVENTS+1):
        _ = sched.step(mem, env)

        # optional crowding decay
        if CROWDING_DECAY > 0 and mem.crowding_count > 0:
            mem.crowding_count = max(0, mem.crowding_count - CROWDING_DECAY)

        if step in SAVE_STEPS:
            draw_frame(mem, frames_dir / f"frame_{step:04d}.png", title=f"t={step}")
        env.diurnal.tick()

    if TOTAL_EVENTS not in SAVE_STEPS:
        draw_frame(mem, frames_dir / f"frame_{TOTAL_EVENTS:04d}.png", title=f"t={TOTAL_EVENTS}")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "event_log.json").write_text(json.dumps(mem.event_log, indent=2))

    # raft (component) size histogram
    comps = connected_components_full(mem)
    sizes = sorted(len(c) for c in comps) if comps else []
    plt.figure(figsize=(6,4))
    bins = range(1, (max(sizes)+2) if sizes else 2)
    plt.hist(sizes, bins=bins, edgecolor="black")
    plt.xlabel("Raft size (number of peptides, bonds-defined)"); plt.ylabel("Count")
    plt.title("DFlow reversible — raft size distribution")
    plt.savefig(OUT / "raft_size_histogram.png", dpi=150, bbox_inches="tight"); plt.close()

    # CSV summary
    with open(OUT / "raft_sizes.csv", "w", newline="") as f:
        f.write("component_index,size\n")
        for i,comp in enumerate(comps):
            f.write(f"{i},{len(comp)}\n")

    (OUT / "peptide_pool.json").write_text(json.dumps(mem.pool, indent=2))
    (OUT / "crowding.json").write_text(json.dumps({"crowding_count": mem.crowding_count}, indent=2))

# -----------------------------
# CLI
# -----------------------------
def parse_args():
    p = argparse.ArgumentParser(description="DFlow reversible (mismatch window + soft acceptance; crowding; raft 2D Stokes-like diffusion)")
    p.add_argument("--SEED", type=int, default=42)
    p.add_argument("--N0", type=int, default=12)
    p.add_argument("--HEX_RADIUS", type=float, default=1.0)
    p.add_argument("--TOTAL_EVENTS", type=int, default=1500)
    p.add_argument("--OUT", type=str, default="runs/exp_phys")
    p.add_argument("--SAVE_STEPS", type=int, nargs="*", default=[0, 750, 1500])

    # Diurnal / bias params
    p.add_argument("--DAY_STEPS", type=int, default=10)
    p.add_argument("--NIGHT_STEPS", type=int, default=10)
    p.add_argument("--POLY_GAIN_DAY", type=float, default=10.0)
    p.add_argument("--POLY_GAIN_NIGHT", type=float, default=0.2)

    # Energetics temp factor (β = 1/kT) placeholder for soft gate
    p.add_argument("--BETA", type=float, default=1.0)
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    out = Path(args.OUT); out.mkdir(parents=True, exist_ok=True)
    run_sim(
        N0=args.N0,
        TOTAL_EVENTS=args.TOTAL_EVENTS,
        OUT=out,
        SAVE_STEPS=list(args.SAVE_STEPS),
        SEED=args.SEED,
        HEX_RADIUS=args.HEX_RADIUS,
        DAY_STEPS=args.DAY_STEPS,
        NIGHT_STEPS=args.NIGHT_STEPS,
        POLY_GAIN_DAY=args.POLY_GAIN_DAY,
        POLY_GAIN_NIGHT=args.POLY_GAIN_NIGHT,
        BETA_arg=args.BETA,
    )
    print(f"Artifacts in: {out.resolve()}")

