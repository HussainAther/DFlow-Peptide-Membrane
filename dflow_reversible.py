# dflow_reversible.py
# -*- coding: utf-8 -*-
"""
DFlow peptide–membrane simulation (single-file edition)
- Hexagonal lattice (axial coords, flat-top).
- Amphiphiles populate cells with a "carbon length" proxy (10..20 -> 2..4 nm thickness).
- Peptides are triads (three adjacent hexes); optional α-helix segment length 3..12 aa.
- Insertion is gated ONLY by vertical hydrophobic mismatch (no lateral packing).
- Rafts = touching peptides; diffusion scales ~ 1/(eta * raft_size) (2D Stokes-like).
- Desorption routes to pool or crowding per CLI flag; failed insertion 50/50 to pool/crowding.
- Periodic saves; logs; raft-size histogram; mass accounting snapshot.

This is a compact, dependency-light baseline intended to run out-of-the-box.
"""

import os, json, math, random, argparse, csv
from collections import deque, defaultdict
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use("Agg")  # save-only
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

# -----------------------------
# Globals / Config defaults
# -----------------------------
# SEED = None
# random.seed(SEED); np.random.seed(SEED)

N0 = 12                   # rhombus radius for axial coords
HEX_RADIUS = 1.0
OUT = "runs/exp"
TOTAL_EVENTS = 1500
SAVE_STEPS = [0, 750, 1500]

DAY_STEPS = 10
NIGHT_STEPS = 10
POLY_GAIN_DAY = 10.0
POLY_GAIN_NIGHT = 0.2
BETA = 1.0

# New knobs
# DI_CARBOXY_FRAC = 0.35       # dihead fraction in medium
# ETA_VISCOSITY = 1.0          # affects diffusion speed inversely with raft size
# A_INSERT_MISMATCH_TOL = 0.5  # nm, hard gate
# DISCARD_TO_CROWDING_P = 0.5  # failed insertion routing prob to crowding

# DESORB_TO_CROWDING = False   # CLI flag

EVENT_WEIGHTS = {
    "amph_swap": 1.0,            # fluidity
    "amph_addition": 0.5,        # A2: mono pair or di single adjustment
    "pept_insert": 0.35,         # insertion attempt (gated by mismatch)
    "pept_desorb": 0.02,         # remove a peptide (routes per flag)
    "pept_diffuse": 0.8,         # motion (weighted by 1/(eta*raft_size))
    "raft_adhere_merge": 0.6,    # merge touching rafts (ID unification)
    "pept_flip": 0.2,            # flip inside/outside
}

# -----------------------------
# Hex helpers (axial coords)
# -----------------------------
SQRT3 = math.sqrt(3.0)
AX_NEI = [(1,0),(1,-1),(0,-1),(-1,0),(-1,1),(0,1)]

def axial_to_xy(q, r, s=HEX_RADIUS):
    x = s * 1.5 * q
    y = s * SQRT3 * (r + 0.5*q)
    return x, y

# Triad shapes (3 hexes around a vertex)
TRIAD_SHAPES = [
    [(0,0),(1,0),(0,1)],
    [(0,0),(0,1),(-1,1)],
    [(0,0),(-1,1),(-1,0)],
    [(0,0),(-1,0),(0,-1)],
    [(0,0),(0,-1),(1,-1)],
    [(0,0),(1,-1),(1,0)]
]

def rotate_triad(anchor, orient):
    aq, ar = anchor
    offs = TRIAD_SHAPES[orient % 6]
    return [(aq+dq, ar+dr) for (dq,dr) in offs]

# -----------------------------
# Color palette helper
# -----------------------------
def get_palette(n, name="tab20"):
    try:
        cmap = mpl.colormaps[name]
    except AttributeError:
        cmap = plt.get_cmap(name)
    if hasattr(cmap, "colors"):
        base = list(cmap.colors)
        if n <= len(base):
            return base[:n]
        reps = (n + len(base) - 1)//len(base)
        return (base * reps)[:n]
    if n == 1:
        return [cmap(0.0)]
    return [cmap(i/(n-1)) for i in range(n)]

# -----------------------------
# Membrane
# -----------------------------
class Membrane:
    def __init__(self, n=N0):
        self.n = n
        self.amph = {}           # (q,r) -> {"carbon":int, "pep":bool}
        self.peptides = {}       # pid -> {"cent":(q,r), "orient":int, "raft":int, "inside":bool, "length":int, "helical":bool}
        self.pid_counter = 0
        self.raft_counter = 0
        self.displaced_amph = 0
        self.event_log = []
        self.crowding_count = 0
        self.pool = {}           # chirality -> count (placeholder uses "0")
        self._init_membrane()

    # geometry region (rhombus in axial coords)
    def _axial_in_rhombus(self, q, r):
        n = self.n
        return (abs(q) <= n) and (abs(r) <= n) and (abs(q+r) <= n)

    def _init_membrane(self):
        # Start with a filled rhombus of amphiphiles (10..20 carbons)
        for q in range(-self.n, self.n+1):
            for r in range(-self.n, self.n+1):
                if self._axial_in_rhombus(q,r):
                    self.amph[(q,r)] = {"carbon": random.randint(10,20), "pep": False}

    # ---- triad helpers ----
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

    # ---- thickness & hydrophobic length ----
    @staticmethod
    def peptide_hydrophobic_length_nm(length_aa:int)->float:
        # ≈0.15 nm per residue (helix rise)
        return 0.15*max(0, length_aa)

    def local_membrane_thickness_nm(self, cent):
        # map carbon 10..20 -> thickness 2..4 nm
        cells = self.triad_cells(cent, 0)
        vals = []
        for c in cells:
            if c in self.amph:
                C = self.amph[c]["carbon"]
                t = 2.0 + (C-10)*(2.0/10.0)  # 10->2nm, 20->4nm
                vals.append(t)
        return sum(vals)/len(vals) if vals else 3.0

    # -----------------------------
    # Events
    # -----------------------------
    def event_amph_swap(self):
        cells = list(self.amph.keys())
        if not cells: return
        a = random.choice(cells)
        random.shuffle(AX_NEI)
        for dq,dr in AX_NEI:
            b = (a[0]+dq, a[1]+dr)
            if b in self.amph:
                ca, cb = self.amph[a]["carbon"], self.amph[b]["carbon"]
                self.amph[a]["carbon"], self.amph[b]["carbon"] = cb, ca
                self.event_log.append({"evt":"amph_swap","a":a,"b":b})
                return

    # A2: amphiphile addition/adjustment (di single / mono pair)
    def sample_amphiphile_species(self):
        chain = random.choices([10,11,12,13,14,15,16,17,18,19,20],
                       weights=[1,2,3,4,5,5,4,3,2,1,1])[0]
        is_di = (random.random() < DI_CARBOXY_FRAC)
        return ('di' if is_di else 'mono', chain)

    def amph_add_single(self, q, r, chain, tag='di'):
        if (q,r) in self.amph and not self.amph[(q,r)]["pep"]:
            self.amph[(q,r)]["carbon"] = chain
            self.event_log.append({"evt":"A2_add","kind":tag,"cell":(q,r),"chain":chain})
            return True
        return False

    def amph_add_pair(self, chain):
        cells = [c for c,v in self.amph.items() if not v["pep"]]
        random.shuffle(cells)
        for q,r in cells:
            random.shuffle(AX_NEI)
            for dq,dr in AX_NEI:
                nb = (q+dq, r+dr)
                if nb in self.amph and not self.amph[nb]["pep"]:
                    self.amph[(q,r)]["carbon"] = chain
                    self.amph[nb]["carbon"] = chain
                    self.event_log.append({"evt":"A2_add","kind":"mono_pair","cells":[(q,r),nb],"chain":chain})
                    return True
        return False

    def event_amph_addition(self):
        kind, chain = self.sample_amphiphile_species()
        if kind == 'di':
            free = [c for c,v in self.amph.items() if not v["pep"]]
            if not free: return
            q,r = random.choice(free)
            self.amph_add_single(q,r,chain,'di')
        else:
            self.amph_add_pair(chain)

    def event_peptide_insert(self):
        if len(self.amph)==0: return
        anchor, orient, cells = self.random_free_anchor()
        if anchor is None: 
            # cannot place anywhere
            self.event_log.append({"evt":"pept_insert_skip","reason":"no_space"})
            return
        # generate candidate peptide
        length = random.randint(3,12)
        helical = (random.random()<0.9)
        Lp = self.peptide_hydrophobic_length_nm(length)
        tloc = self.local_membrane_thickness_nm(anchor)
        leaflet = 0.5 * tloc
        mismatch = abs(Lp - leaflet)
        # gate by vertical mismatch only
        if helical and (mismatch <= A_INSERT_MISMATCH_TOL) and self.cells_free_for_triad(cells):
            displaced = sum(1 for c in cells if (c in self.amph and not self.amph[c]["pep"]))
            self.displaced_amph += displaced
            pid = self.pid_counter; self.pid_counter += 1
            rid = self.raft_counter; self.raft_counter += 1
            self.peptides[pid] = {"cent":anchor,"orient":orient,"raft":rid,"inside":bool(random.getrandbits(1)),
                                  "length":length,"helical":helical}
            for c in cells:
                self.amph[c]["pep"] = True
            self.event_log.append({"evt":"pept_insert","pid":pid,"cells":cells,"displaced_inc":displaced,
                                   "length":length,"Lp_nm":Lp,"mismatch_nm":mismatch})
            # growth trigger
            needed = 4*self.n + 2
            if self.displaced_amph >= needed:
                self._grow_membrane()
                self.displaced_amph = 0
                self.event_log.append({"evt":"grow","new_n":self.n})
        else:
            # failed insertion: 50% crowding vs pool (independent of DESORB flag)
            to_crowd = (random.random() < DISCARD_TO_CROWDING_P)
            if to_crowd:
                self.crowding_count += 1
            else:
                self.pool["0"] = self.pool.get("0", 0) + 1
            self.event_log.append({"evt":"pept_insert_reject","helical":helical,"length":length,
                                   "Lp_nm":Lp,"mismatch_nm":mismatch,"dest":("crowding" if to_crowd else "pool")})

    def _grow_membrane(self):
        new_n = self.n + 2
        for q in range(-new_n, new_n+1):
            for r in range(-new_n, new_n+1):
                if self._axial_in_rhombus(q,r) and (q,r) not in self.amph:
                    self.amph[(q,r)] = {"carbon": random.randint(10,20), "pep": False}
        self.n = new_n

    def event_peptide_diffuse(self):
        if not self.peptides: return
        # compute raft sizes by current raft ids
        raft_sizes = defaultdict(int)
        for p in self.peptides.values():
            raft_sizes[p["raft"]] += 1
        # weight selection by 1/(eta * size)
        pids = list(self.peptides.keys())
        weights = []
        for pid in pids:
            size = max(1, raft_sizes[self.peptides[pid]["raft"]])
            weights.append(1.0/(ETA_VISCOSITY*size))
        pid = random.choices(pids, weights=weights, k=1)[0]
        p = self.peptides[pid]
        q,r = p["cent"]
        dq,dr = random.choice(AX_NEI)
        new_cent = (q+dq, r+dr)
        new_cells = self.triad_cells(new_cent, p["orient"])
        if not self.cells_free_for_triad(new_cells): return
        # free old
        for c in self.triad_cells(p["cent"], p["orient"]):
            if c in self.amph: self.amph[c]["pep"] = False
        # occupy new
        for c in new_cells:
            self.amph[c]["pep"] = True
        p["cent"] = new_cent
        self.event_log.append({"evt":"pept_diffuse","pid":pid,"from":(q,r),"to":new_cent,
                               "raft_size": raft_sizes[p["raft"]]})

    def event_raft_adhere_merge(self):
        # merge touching rafts (unify raft IDs when triad cells touch)
        if len(self.peptides)<2: return
        # map cell->pid
        cell_to_pid = {}
        for pid, p in self.peptides.items():
            for c in self.triad_cells(p["cent"], p["orient"]):
                cell_to_pid[c] = pid
        pairs = set()
        for c, pid in cell_to_pid.items():
            q,r = c
            for dq,dr in AX_NEI:
                nb = (q+dq, r+dr)
                if nb in cell_to_pid:
                    pid2 = cell_to_pid[nb]
                    if pid2 != pid:
                        r1, r2 = self.peptides[pid]["raft"], self.peptides[pid2]["raft"]
                        if r1 != r2:
                            pairs.add(tuple(sorted((r1,r2))))
        if not pairs: return
        rA, rB = random.choice(list(pairs))
        target, source = (rA, rB) if rA<rB else (rB, rA)
        for pid, p in self.peptides.items():
            if p["raft"] == source:
                p["raft"] = target
        self.event_log.append({"evt":"raft_merge","from":source,"into":target})

    def event_pept_flip(self):
        if not self.peptides: return
        pid = random.choice(list(self.peptides.keys()))
        p = self.peptides[pid]
        p["inside"] = not p["inside"]
        self.event_log.append({"evt":"pept_flip","pid":pid,"inside":p["inside"]})

    def event_pept_desorb(self):
        if not self.peptides: return
        pid = random.choice(list(self.peptides.keys()))
        p = self.peptides.pop(pid, None)
        if p is None: return
        # free cells
        for c in self.triad_cells(p["cent"], p["orient"]):
            if c in self.amph: self.amph[c]["pep"] = False
        # route
        dest = "crowding" if DESORB_TO_CROWDING else "pool"
        if dest == "crowding":
            self.crowding_count += 1
        else:
            self.pool["0"] = self.pool.get("0", 0) + 1
        self.event_log.append({"evt":"P2_minus_desorb","pid":pid,"dest":dest})

    # -----------------------------
    # Analysis helpers
    # -----------------------------
    def raft_components(self):
        """Return list of components (each as set of pids) where peptides touch by cell adjacency."""
        # Build adjacency graph over pids
        pids = list(self.peptides.keys())
        if not pids: return []
        cellsets = {pid: set(self.triad_cells(p["cent"], p["orient"])) for pid,p in self.peptides.items()}
        adj = {pid:set() for pid in pids}
        for i in range(len(pids)):
            a = pids[i]
            for j in range(i+1, len(pids)):
                b = pids[j]
                # check touch: any cell in a adjacent to any cell in b
                touch = False
                for (q,r) in cellsets[a]:
                    for dq,dr in AX_NEI:
                        if (q+dq, r+dr) in cellsets[b]:
                            touch = True; break
                    if touch: break
                if touch:
                    adj[a].add(b); adj[b].add(a)
        # BFS components
        comps = []
        seen = set()
        for pid in pids:
            if pid in seen: continue
            comp = set([pid])
            dq = deque([pid])
            seen.add(pid)
            while dq:
                u = dq.popleft()
                for v in adj[u]:
                    if v not in seen:
                        seen.add(v); comp.add(v); dq.append(v)
            comps.append(comp)
        return comps

# -----------------------------
# Drawing / Saving
# -----------------------------
def raft_color_map(peptides):
    raft_ids = sorted({p["raft"] for p in peptides.values()}) if peptides else []
    palette = get_palette(max(1, len(raft_ids)), "tab20")
    return {rid: palette[i] for i,rid in enumerate(raft_ids)}

def save_frame(mem:Membrane, out_path:Path, idx:int, border:float=1.0, title:str=None):
    fig, ax = plt.subplots(figsize=(10,10))
    ax.set_aspect('equal'); ax.axis('off')
    # Amphiphiles
    minC,maxC = 8.0,22.0
    for (q,r),info in mem.amph.items():
        x,y = axial_to_xy(q,r,HEX_RADIUS)
        t = (info["carbon"] - minC)/(maxC-minC); t = max(0.0, min(1.0, t))
        base = 0.92 - 0.25*t
        patch = RegularPolygon((x,y), numVertices=6, radius=HEX_RADIUS, orientation=math.radians(30),
                               facecolor=(base,base,base), edgecolor='k', linewidth=0.5)
        ax.add_patch(patch)
    # Peptides
    cmap = raft_color_map(mem.peptides)
    for pid,p in mem.peptides.items():
        cells = mem.triad_cells(p["cent"], p["orient"])
        fc = cmap.get(p["raft"], (0.8, 0.2, 0.2))
        lw = 1.5 if p.get("inside",False) else 0.6
        for (q,r) in cells:
            x,y = axial_to_xy(q,r,HEX_RADIUS)
            patch = RegularPolygon((x,y), numVertices=6, radius=HEX_RADIUS, orientation=math.radians(30),
                                   facecolor=fc, edgecolor='black', linewidth=lw)
            ax.add_patch(patch)

    # bounds
    qmin,qmax = -mem.n, mem.n; rmin,rmax = -mem.n, mem.n
    corners = [
        axial_to_xy(qmin,rmin), axial_to_xy(qmin,rmax),
        axial_to_xy(qmax,rmin), axial_to_xy(qmax,rmax),
        axial_to_xy(0,-mem.n), axial_to_xy(0,mem.n),
        axial_to_xy(-mem.n,0), axial_to_xy(mem.n,0)
    ]
    xs = [c[0] for c in corners]; ys = [c[1] for c in corners]
    ax.set_xlim(min(xs)-border, max(xs)+border)
    ax.set_ylim(min(ys)-border, max(ys)+border)
    if title:
        ax.set_title(title)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight"); plt.close(fig)

def write_raft_hist(mem:Membrane, out_dir:Path):
    comps = mem.raft_components()
    sizes = [len(c) for c in comps] if comps else []
    # CSV
    with open(out_dir/"raft_sizes.csv","w",newline="") as f:
        w = csv.writer(f); w.writerow(["size"]); [w.writerow([s]) for s in sizes]
    # Plot
    if sizes:
        plt.figure()
        bins = list(range(1, max(sizes)+2))
        plt.hist(sizes, bins=bins, align='left')
        plt.xlabel("Raft size (# peptides)"); plt.ylabel("Count")
        plt.title("Raft size distribution")
        plt.savefig(out_dir/"raft_size_histogram.png", dpi=150, bbox_inches="tight")
        plt.close()

def mass_snapshot(mem:Membrane):
    pool_total = sum(mem.pool.values()) if mem.pool else 0
    crowd = mem.crowding_count
    embedded = len(mem.peptides)
    return {"pool":pool_total,"crowding":crowd,"embedded":embedded,"sum":pool_total+crowd+embedded}

# -----------------------------
# Case loop
# -----------------------------
def pick_event():
    items = list(EVENT_WEIGHTS.items())
    names = [k for k,_ in items]
    w = np.array([v for _,v in items], dtype=float)
    w = w / w.sum()
    return np.random.choice(names, p=w)

def run_sim(N0:int=N0, HEX_RADIUS_:float=HEX_RADIUS, TOTAL_EVENTS_:int=TOTAL_EVENTS,
            OUT_:str=OUT, SAVE_STEPS_:list=SAVE_STEPS, **kwargs):
    global HEX_RADIUS
    HEX_RADIUS = HEX_RADIUS_
    out_dir = Path(OUT_)
    mem = Membrane(n=N0)
    # log meta
    mem.event_log.append({"meta":"run","seed":SEED,"params":{
        "N0":N0,"TOTAL_EVENTS":TOTAL_EVENTS_,"A_INSERT_MISMATCH_TOL":A_INSERT_MISMATCH_TOL,
        "DI_CARBOXY_FRAC":DI_CARBOXY_FRAC,"ETA_VISCOSITY":ETA_VISCOSITY,
        "DESORB_TO_CROWDING":DESORB_TO_CROWDING
    }})
    # initial save
    if 0 in SAVE_STEPS_:
        save_frame(mem, out_dir/f"frame_{0:04d}.png", 0, title="t=0")
    # loop
    for step in range(1, TOTAL_EVENTS_+1):
        evt = pick_event()
        if evt == "amph_swap":
            mem.event_amph_swap()
        elif evt == "amph_addition":
            mem.event_amph_addition()
        elif evt == "pept_insert":
            mem.event_peptide_insert()
        elif evt == "pept_desorb":
            mem.event_pept_desorb()
        elif evt == "pept_diffuse":
            mem.event_peptide_diffuse()
        elif evt == "raft_adhere_merge":
            mem.event_raft_adhere_merge()
        elif evt == "pept_flip":
            mem.event_pept_flip()
        # saves
        if step in SAVE_STEPS_:
            save_frame(mem, out_dir/f"frame_{step:04d}.png", step, title=f"t={step}")
    # final outputs
    save_frame(mem, out_dir/f"frame_{TOTAL_EVENTS_:04d}.png", TOTAL_EVENTS_, title=f"t={TOTAL_EVENTS_}")
    write_raft_hist(mem, out_dir)
    with open(out_dir/"event_log.json","w") as f:
        json.dump(mem.event_log, f, indent=2)
    with open(out_dir/"mass_check.json","w") as f:
        json.dump(mass_snapshot(mem), f, indent=2)
    print(f"Artifacts in: {out_dir}")

# -----------------------------
# CLI
# -----------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--SEED", type=int, default=42)
    p.add_argument("--N0", type=int, default=N0)
    p.add_argument("--HEX_RADIUS", type=float, default=HEX_RADIUS)
    p.add_argument("--TOTAL_EVENTS", type=int, default=TOTAL_EVENTS)
    p.add_argument("--OUT", type=str, default=OUT)
    p.add_argument("--SAVE_STEPS", type=int, nargs='*', default=SAVE_STEPS)
    p.add_argument("--DAY_STEPS", type=int, default=DAY_STEPS)
    p.add_argument("--NIGHT_STEPS", type=int, default=NIGHT_STEPS)
    p.add_argument("--POLY_GAIN_DAY", type=float, default=POLY_GAIN_DAY)
    p.add_argument("--POLY_GAIN_NIGHT", type=float, default=POLY_GAIN_NIGHT)
    p.add_argument("--BETA", type=float, default=BETA)
    p.add_argument("--DESORB_TO_CROWDING", action="store_true",
                   help="Route desorbed peptides to crowding instead of returning to pool.")
    p.add_argument("--DI_CARBOXY_FRAC", type=float, default=0.35, # <-- CHANGED
                   help="Fraction of dihead amphiphiles in medium [0..1].")
    p.add_argument("--ETA_VISCOSITY", type=float, default=1.0,  # <-- CHANGED
                   help="Effective 2D viscosity factor for raft diffusion.")
    p.add_argument("--A_INSERT_MISMATCH_TOL", type=float, default=1.0, # <-- CHANGED
                   help="Hard vertical mismatch tolerance (nm).")
    p.add_argument("--DISCARD_TO_CROWDING_P", type=float, default=0.5, # <-- CHANGED
                   help="Failed insert: probability to route to crowding (else pool).")
    return p.parse_args()

if __name__ == "__main__":
    global SEED, DESORB_TO_CROWDING, DI_CARBOXY_FRAC, ETA_VISCOSITY, A_INSERT_MISMATCH_TOL, DISCARD_TO_CROWDING_P
    args = parse_args()
    # seeds & knobs
    current_seed = int(args.SEED)
    random.seed(current_seed)
    np.random.seed(current_seed)
    SEED = current_seed # Assign to global SEED *after* seeding is done
    DESORB_TO_CROWDING = bool(args.DESORB_TO_CROWDING)
    DI_CARBOXY_FRAC       = float(args.DI_CARBOXY_FRAC)
    ETA_VISCOSITY         = float(args.ETA_VISCOSITY)
    A_INSERT_MISMATCH_TOL = float(args.A_INSERT_MISMATCH_TOL)
    DISCARD_TO_CROWDING_P = float(args.DISCARD_TO_CROWDING_P)

    out = Path(args.OUT)
    out.mkdir(parents=True, exist_ok=True)
    run_sim(
        N0=args.N0,
        HEX_RADIUS_=args.HEX_RADIUS,
        TOTAL_EVENTS_=args.TOTAL_EVENTS,
        OUT_=args.OUT,
        SAVE_STEPS_=args.SAVE_STEPS,
        DAY_STEPS=args.DAY_STEPS,
        NIGHT_STEPS=args.NIGHT_STEPS,
        POLY_GAIN_DAY=args.POLY_GAIN_DAY,
        POLY_GAIN_NIGHT=args.POLY_GAIN_NIGHT,
        BETA=args.BETA,
    )

