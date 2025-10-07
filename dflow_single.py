#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DFlow peptide–membrane simulation (single-file)
- Hex lattice membrane (flat-top)
- Amphiphile swap, peptide insert/diffuse/flip
- Raft adhesion+merge
- Growth when displaced amph ≥ 4n+2
- Frames + event log + raft-size histogram

Usage:
  python dflow_single.py --N0 12 --TOTAL_EVENTS 1000 --SAVE_STEPS 0 500 1000 --OUT runs/exp_0001
"""
import os, json, math, random, argparse
from pathlib import Path
from collections import Counter

import numpy as np
import matplotlib as mpl
mpl.use("Agg")  # headless save-only
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

# -----------------------------
# Geometry helpers (flat-top hex axial coords)
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
    try:
        cmap = mpl.colormaps[name]          # Matplotlib >=3.7
    except AttributeError:
        cmap = plt.get_cmap(name)           # older
    if hasattr(cmap, "colors"):             # listed cmap (tab10/tab20/…)
        base = list(cmap.colors)
        if n <= len(base): return base[:n]
        reps = (n + len(base) - 1)//len(base)
        return (base*reps)[:n]
    if n == 1: return [cmap(0.0)]
    return [cmap(i/(n-1)) for i in range(n)]

# -----------------------------
# State
# -----------------------------
class Membrane:
    def __init__(self, n:int=12, hex_radius:float=1.0, label_carbons:bool=False):
        self.n = n
        self.hex_radius = hex_radius
        self.label_carbons = label_carbons
        self.amph = {}        # (q,r) -> {"carbon":int, "pep":bool}
        self.peptides = {}    # pid -> {"cent":(q,r), "orient":int, "raft":int, "inside":bool}
        self.pid_counter = 0
        self.raft_counter = 0
        self.displaced_amph = 0
        self.event_log = []
        self._init_membrane()

    def _axial_in_rhombus(self, q, r):
        n = self.n
        return (abs(q) <= n) and (abs(r) <= n) and (abs(q+r) <= n)

    def _init_membrane(self):
        for q in range(-self.n, self.n+1):
            for r in range(-self.n, self.n+1):
                if self._axial_in_rhombus(q,r):
                    self.amph[(q,r)] = {"carbon": random.randint(10,20), "pep": False}

    # queries
    def triad_cells(self, center, orient):
        return rotate_triad(center, orient)

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

    # growth
    def _grow_membrane(self):
        new_n = self.n + 2
        for q in range(-new_n, new_n+1):
            for r in range(-new_n, new_n+1):
                if self._axial_in_rhombus(q,r) and (q,r) not in self.amph:
                    self.amph[(q,r)] = {"carbon": random.randint(10,20), "pep": False}
        self.n = new_n

# -----------------------------
# Event registry
# -----------------------------
class EventRegistry:
    def __init__(self):
        self._events = {}   # name -> func(mem)
        self._weights = {}  # name -> float

    def register(self, name, func, weight:float=1.0):
        self._events[name] = func
        self._weights[name] = float(weight)

    def pick(self):
        names = list(self._events.keys())
        w = np.array([self._weights[n] for n in names], dtype=float)
        w = w / (w.sum() if w.sum() else 1.0)
        return np.random.choice(names, p=w)

    def run(self, name, mem:Membrane):
        self._events[name](mem)

# -----------------------------
# Baseline events
# -----------------------------
def event_amph_swap(mem:Membrane):
    cells = list(mem.amph.keys())
    if not cells: return
    a = random.choice(cells)
    q,r = a
    random.shuffle(AX_NEI)
    for dq,dr in AX_NEI:
        b = (q+dq, r+dr)
        if b in mem.amph:
            mem.amph[a]["carbon"], mem.amph[b]["carbon"] = mem.amph[b]["carbon"], mem.amph[a]["carbon"]
            mem.event_log.append({"evt":"amph_swap","a":a,"b":b})
            return

def event_peptide_insert(mem:Membrane):
    if len(mem.amph)==0: return
    anchor, orient, cells = mem.random_free_anchor()
    if anchor is None: return
    displaced = sum(1 for c in cells if c in mem.amph and not mem.amph[c]["pep"])
    mem.displaced_amph += displaced
    pid = mem.pid_counter; mem.pid_counter += 1
    rid = mem.raft_counter; mem.raft_counter += 1
    mem.peptides[pid] = {"cent": anchor, "orient": orient, "raft": rid, "inside": bool(random.getrandbits(1))}
    for c in cells: mem.amph[c]["pep"] = True
    mem.event_log.append({"evt":"pept_insert","pid":pid,"cells":cells,"displaced_inc":displaced})
    needed = 4*mem.n + 2
    if mem.displaced_amph >= needed:
        mem._grow_membrane()
        mem.displaced_amph = 0
        mem.event_log.append({"evt":"grow","new_n":mem.n})

def event_peptide_diffuse(mem:Membrane):
    if not mem.peptides: return
    pid = random.choice(list(mem.peptides.keys()))
    p = mem.peptides[pid]
    q,r = p["cent"]
    dq,dr = random.choice(AX_NEI)
    new_cent = (q+dq, r+dr)
    new_cells = mem.triad_cells(new_cent, p["orient"])
    if not mem.cells_free_for_triad(new_cells): return
    # free old
    for c in mem.triad_cells(p["cent"], p["orient"]):
        if c in mem.amph: mem.amph[c]["pep"] = False
    # occupy new
    for c in new_cells: mem.amph[c]["pep"] = True
    p["cent"] = new_cent
    mem.event_log.append({"evt":"pept_diffuse","pid":pid,"from":(q,r),"to":new_cent})

def event_raft_adhere_merge(mem:Membrane):
    if len(mem.peptides) < 2: return
    # map cell -> pid
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
                    r1 = mem.peptides[pid]["raft"]
                    r2 = mem.peptides[pid2]["raft"]
                    if r1 != r2:
                        pairs.add(tuple(sorted((r1,r2))))
    if not pairs: return
    rA, rB = random.choice(list(pairs))
    target, source = (rA, rB) if rA < rB else (rB, rA)
    for pid, p in mem.peptides.items():
        if p["raft"] == source:
            p["raft"] = target
    mem.event_log.append({"evt":"raft_merge","from":source,"into":target})

def event_pept_flip(mem:Membrane):
    if not mem.peptides: return
    pid = random.choice(list(mem.peptides.keys()))
    p = mem.peptides[pid]
    p["inside"] = not p["inside"]
    mem.event_log.append({"evt":"pept_flip","pid":pid,"inside":p["inside"]})

# --- Optional new/stub event for collaborators (safe placeholder) ---
def event_pept_desorb(mem:Membrane):
    """Detach a random peptide and free its triad (simple version)."""
    if not mem.peptides: return
    pid = random.choice(list(mem.peptides.keys()))
    p = mem.peptides.pop(pid)
    for c in mem.triad_cells(p["cent"], p["orient"]):
        if c in mem.amph: mem.amph[c]["pep"] = False
    mem.event_log.append({"evt":"pept_desorb","pid":pid})

# -----------------------------
# Harness / visualization
# -----------------------------
def get_registry():
    reg = EventRegistry()
    reg.register("amph_swap",          event_amph_swap,          weight=1.0)
    reg.register("pept_insert",        event_peptide_insert,     weight=0.15)
    reg.register("pept_diffuse",       event_peptide_diffuse,    weight=0.8)
    reg.register("raft_adhere_merge",  event_raft_adhere_merge,  weight=0.6)
    reg.register("pept_flip",          event_pept_flip,          weight=0.2)
    # reg.register("pept_desorb",      event_pept_desorb,        weight=0.10)  # enable when ready
    return reg

def draw_frame(mem:Membrane, path:Path, border:float=1.0, title:str|None=None):
    fig, ax = plt.subplots(figsize=(8,8))
    ax.set_aspect('equal'); ax.axis('off')
    minC, maxC = 8.0, 22.0
    # amphiphiles
    for (q,r), info in mem.amph.items():
        x,y = axial_to_xy(q,r, mem.hex_radius)
        t = (info["carbon"] - minC) / (maxC - minC)
        t = min(max(t, 0.0), 1.0)
        base = 0.92 - 0.25*t
        patch = RegularPolygon(
            (x,y), numVertices=6, radius=mem.hex_radius, orientation=math.radians(30),
            facecolor=(base,base,base), edgecolor='k', linewidth=0.5
        )
        ax.add_patch(patch)
        if mem.label_carbons:
            ax.text(x, y, str(info["carbon"]), ha='center', va='center', fontsize=6)

    # peptides
    raft_ids = sorted({p["raft"] for p in mem.peptides.values()}) if mem.peptides else []
    rmap = {rid: get_palette(max(1,len(raft_ids)))[i] for i, rid in enumerate(raft_ids)}
    for pid, p in mem.peptides.items():
        cells = mem.triad_cells(p["cent"], p["orient"])
        fc = rmap.get(p["raft"], (0.8,0.2,0.2))
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

def run_sim(N0:int, TOTAL_EVENTS:int, SAVE_STEPS:list[int], OUT:Path, SEED:int=42, HEX_RADIUS:float=1.0):
    # seeds
    random.seed(SEED); np.random.seed(SEED)
    # state + registry
    mem = Membrane(n=N0, hex_radius=HEX_RADIUS)
    reg = get_registry()

    frames_dir = OUT / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    # t=0 frame
    draw_frame(mem, frames_dir / "frame_0000.png", title="t=0")

    # main loop
    for step in range(1, TOTAL_EVENTS+1):
        evt = reg.pick()
        reg.run(evt, mem)
        if step in SAVE_STEPS:
            draw_frame(mem, frames_dir / f"frame_{step:04d}.png", title=f"t={step}")

    # final snapshot at TOTAL_EVENTS (even if not in SAVE_STEPS)
    if TOTAL_EVENTS not in SAVE_STEPS:
        draw_frame(mem, frames_dir / f"frame_{TOTAL_EVENTS:04d}.png", title=f"t={TOTAL_EVENTS}")

    # logs
    (OUT / "event_log.json").write_text(json.dumps(mem.event_log, indent=2))

    # raft histogram
    rc = Counter(p["raft"] for p in mem.peptides.values())
    sizes = sorted(rc.values())
    plt.figure(figsize=(6,4))
    bins = range(1, (max(sizes)+2) if sizes else 2)
    plt.hist(sizes, bins=bins, edgecolor="black")
    plt.xlabel("Raft size (number of peptides)"); plt.ylabel("Count")
    plt.title("DFlow raft size distribution")
    plt.savefig(OUT / "raft_size_histogram.png", dpi=150, bbox_inches="tight"); plt.close()

    # raft sizes CSV
    with open(OUT / "raft_sizes.csv", "w", newline="") as f:
        f.write("raft_id,size\n")
        for rid, sz in sorted(rc.items()):
            f.write(f"{rid},{sz}\n")

# -----------------------------
# CLI
# -----------------------------
def parse_args():
    p = argparse.ArgumentParser(description="DFlow peptide–membrane (single-file)")
    p.add_argument("--SEED", type=int, default=42)
    p.add_argument("--N0", type=int, default=12, help="initial rhombus half-width n")
    p.add_argument("--HEX_RADIUS", type=float, default=1.0)
    p.add_argument("--TOTAL_EVENTS", type=int, default=1000)
    p.add_argument("--SAVE_STEPS", type=int, nargs="*", default=[500,1000])
    p.add_argument("--OUT", type=str, default="runs/exp_0001")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    out = Path(args.OUT)
    run_sim(
        N0=args.N0,
        TOTAL_EVENTS=args.TOTAL_EVENTS,
        SAVE_STEPS=args.SAVE_STEPS,
        OUT=out,
        SEED=args.SEED,
        HEX_RADIUS=args.HEX_RADIUS,
    )
    print(f"Artifacts in: {out.resolve()}")

