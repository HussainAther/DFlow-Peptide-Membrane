# membrane_case_engine.py
# DFlow Peptide Membrane — 8/27/2025 meeting spec
# - Edge-touching pointy-top hex grid (no gaps/overlaps)
# - Global Membrane state shared by all events (case statement style)
# - PBC wrapping
# - Growth by adding a full border ring (rim tracker)
# - Lipid carbon counts labeled on each lipid hex
# - Peptide insertion as triads (3 touching hexes sharing a vertex)
# - Saves frames only; no GUI rendering

import os, json, random
import numpy as np
import matplotlib
matplotlib.use("Agg")  # save-to-file only
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

# ---------- HEX GEOMETRY (pointy-top, edge-touching) ----------
SQRT3 = np.sqrt(3.0)
HEX_R  = 12.0                   # patch 'radius' (center -> vertex)
ORIENT = np.pi / 6.0            # 30°, pointy-top
# axial neighbor directions (q,r)
DIRS   = [(1,0),(1,-1),(0,-1),(-1,0),(-1,1),(0,1)]

def axial_to_xy(q:int, r:int):
    """
    Exact edge-touching mapping for pointy-top axial coordinates.
      x = s * sqrt(3) * (q + r/2)
      y = s * 3/2 * r
    where s = HEX_R = center->vertex distance.
    """
    x = HEX_R * (SQRT3*q + (SQRT3/2.0)*r)
    y = HEX_R * (1.5*r)
    return x, y

def triad_cells(anchor, orientation=0):
    """
    Three touching cells around a shared vertex:
      {anchor, anchor + DIR[k], anchor + DIR[k-1]}
    """
    aq, ar = anchor
    k  = orientation % 6
    k2 = (orientation + 5) % 6
    d1 = DIRS[k]
    d2 = DIRS[k2]
    return [(aq,ar), (aq+d1[0], ar+d1[1]), (aq+d2[0], ar+d2[1])]

# ---------- MEMBRANE (global state) ----------
class Membrane:
    """
    Rhombus region in axial coords:
      q in [q0..q1], r in [r0..r1]
    Content:
      - lipids: (q,r) -> {"C": carbon_count}
      - peptides: pid -> {"cells":[(q,r)x3], "anchor":int, "rod":int, "side":+1/-1, "ori":0..5}
    """
    def __init__(self, q0, q1, r0, r1, use_pbc=True, seed=42):
        random.seed(seed)
        self.q0, self.q1 = q0, q1
        self.r0, self.r1 = r0, r1
        self.use_pbc = use_pbc

        self.lipids   = {}
        self.peptides = {}
        self.pid_seq  = 1
        self.displaced_lipids = 0  # counter for growth trigger

        # seed amphiphiles with random carbon counts (10..20)
        for q in range(self.q0, self.q1+1):
            for r in range(self.r0, self.r1+1):
                self.lipids[(q,r)] = {"C": random.randint(10,20)}

        # rim tracking (outer ring of the rhombus)
        self._build_rim()
        self.reset_rim_cursor()

    # ---- extents, wrap, occupancy ----
    @property
    def width(self):  return self.q1 - self.q0 + 1
    @property
    def height(self): return self.r1 - self.r0 + 1

    def in_bounds(self, q, r):
        return self.q0 <= q <= self.q1 and self.r0 <= r <= self.r1

    def wrap(self, q, r):
        if not self.use_pbc:
            return q, r
        w = self.width
        h = self.height
        qn = self.q0 + ((q - self.q0) % w)
        rn = self.r0 + ((r - self.r0) % h)
        return qn, rn

    def is_peptide_cell(self, q, r):
        # Linear scan over few peptides is fine at this scale.
        for p in self.peptides.values():
            if (q,r) in p["cells"]:
                return True
        return False

    def cells_free_for_peptide(self, cells):
        """
        A triad may overwrite lipids (they're displaced),
        but cannot overlap any existing peptide cell.
        All cells must be within bounds.
        """
        seen = set()
        for (q,r) in cells:
            if not self.in_bounds(q,r): return False
            if (q,r) in seen:          return False
            if self.is_peptide_cell(q,r):
                return False
            seen.add((q,r))
        return True

    def remove_lipids(self, cells):
        removed = 0
        for c in cells:
            if c in self.lipids:
                del self.lipids[c]
                removed += 1
        self.displaced_lipids += removed
        return removed

    # ---- rim (border) utilities ----
    def _build_rim(self):
        rim = []
        # top edge
        for q in range(self.q0, self.q1+1): rim.append((q, self.r0))
        # right edge
        for r in range(self.r0+1, self.r1+1): rim.append((self.q1, r))
        # bottom edge
        for q in range(self.q1-1, self.q0-1, -1): rim.append((q, self.r1))
        # left edge
        for r in range(self.r1-1, self.r0, -1): rim.append((self.q0, r))
        self.rim = rim

    def reset_rim_cursor(self):
        self._rim_idx = 0

    def next_rim_site(self):
        if not self.rim: return None
        cell = self.rim[self._rim_idx % len(self.rim)]
        self._rim_idx += 1
        return cell

    # ---- growth: add a full border ring and seed with lipids ----
    def grow_one_ring(self):
        self.q0 -= 1; self.q1 += 1
        self.r0 -= 1; self.r1 += 1
        # seed only the new ring cells
        for q in range(self.q0, self.q1+1):
            for r in range(self.r0, self.r1+1):
                if (q == self.q0 or q == self.q1 or r == self.r0 or r == self.r1):
                    if not self.is_peptide_cell(q,r) and (q,r) not in self.lipids:
                        self.lipids[(q,r)] = {"C": random.randint(10,20)}
        self._build_rim()
        self.reset_rim_cursor()

    def has_free_interior_slot(self):
        # True if there's any cell that is neither a peptide cell nor a lipid
        for q in range(self.q0, self.q1+1):
            for r in range(self.r0, self.r1+1):
                if not self.is_peptide_cell(q,r) and (q,r) not in self.lipids:
                    return True
        return False

# ---------- EVENTS (“case statements”) ----------
def event_swap_lipids(m: Membrane, n_swaps=250):
    """
    Simulate liquid membrane by swapping neighboring amphiphiles randomly.
    """
    keys = list(m.lipids.keys())
    if not keys: return "swap:0"
    swaps = 0
    for _ in range(n_swaps):
        q, r = random.choice(keys)
        dq, dr = random.choice(DIRS)
        q2, r2 = m.wrap(q+dq, r+dr)
        if (q2,r2) in m.lipids:
            m.lipids[(q,r)], m.lipids[(q2,r2)] = m.lipids[(q2,r2)], m.lipids[(q,r)]
            swaps += 1
    return f"swap:{swaps}"

def event_insert_peptide(m: Membrane, tries=400):
    """
    Pick a random *cell* in region (lipid or empty, but not a peptide).
    Try an orientation; if triad fits (no peptide overlap), place peptide.
    Displace underlying lipids in those 3 cells.
    """
    for _ in range(tries):
        q = random.randint(m.q0, m.q1)
        r = random.randint(m.r0, m.r1)
        if m.is_peptide_cell(q,r):  # anchor can't be on an existing peptide cell
            continue
        ori    = random.randrange(6)
        cells  = triad_cells((q,r), ori)
        cells  = [m.wrap(*c) for c in cells] if m.use_pbc else cells
        if not m.cells_free_for_peptide(cells):
            continue
        removed = m.remove_lipids(cells)
        pid = f"P{m.pid_seq:04d}"; m.pid_seq += 1
        m.peptides[pid] = {
            "cells":  cells,
            "anchor": random.randint(4,8),
            "rod":    random.randint(3,8),
            "side":   random.choice([-1,1]),  # -1 inside, +1 outside
            "ori":    ori
        }
        return f"peptide:ok:{pid}:rm{removed}"
    return "peptide:fail"

def event_add_lipid(m: Membrane, tries=200):
    """
    Add a lipid. Prefer an interior free slot; else add to rim; if rim saturated, grow.
    """
    # interior attempt
    for _ in range(tries):
        q = random.randint(m.q0, m.q1)
        r = random.randint(m.r0, m.r1)
        if not m.is_peptide_cell(q,r) and (q,r) not in m.lipids:
            m.lipids[(q,r)] = {"C": random.randint(10,20)}
            return "lipid:add:interior"

    # rim attempt
    for _ in range(len(m.rim) or 1):
        c = m.next_rim_site()
        if c is None: break
        q,r = c
        if not m.is_peptide_cell(q,r) and (q,r) not in m.lipids:
            m.lipids[(q,r)] = {"C": random.randint(10,20)}
            return "lipid:add:rim"

    # growth if rim also saturated
    m.grow_one_ring()
    return "lipid:add:grow"

def event_grow_if_saturated(m: Membrane):
    """
    Growth trigger per 4n+2 rule (n ~ current linear size).
    After growth, reset displaced counter.
    """
    n = max(m.width, m.height)
    thresh = 4*n + 2
    if m.displaced_lipids >= thresh:
        m.grow_one_ring()
        m.displaced_lipids = 0
        return f"grow:ring:+1@n{n}"
    return "grow:none"

# ---------- EVENT SCHEDULER ----------
def weighted_choice(weights):
    total = sum(weights)
    if total <= 0:
        return random.randrange(len(weights))
    r = random.random() * total
    acc = 0.0
    for i,w in enumerate(weights):
        acc += w
        if r <= acc:
            return i
    return len(weights)-1

# ---------- RENDERING ----------
def save_frame(m: Membrane, step, outdir="frames_case_engine", annotate=True, dpi=220):
    os.makedirs(outdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 7), dpi=dpi)
    ax.set_aspect("equal")
    ax.axis("off")

    # draw lipids (white fill, carbon count annotation)
    for (q,r), meta in m.lipids.items():
        x,y = axial_to_xy(q,r)
        ax.add_patch(
            RegularPolygon((x,y), numVertices=6, radius=HEX_R,
                           orientation=ORIENT, ec="0.4", fc="white", lw=0.8)
        )
        if annotate:
            ax.text(x, y, str(meta["C"]), ha="center", va="center",
                    fontsize=6, color="0.35")

    # draw peptides (red = outside, blue = inside)
    for pid, p in m.peptides.items():
        fc = "#d62728" if p["side"]>0 else "#1f77b4"
        for (q,r) in p["cells"]:
            x,y = axial_to_xy(q,r)
            ax.add_patch(
                RegularPolygon((x,y), numVertices=6, radius=HEX_R,
                               orientation=ORIENT, ec="black", fc=fc, lw=1.2)
            )
        # label peptide metadata at anchor cell
        ax.text(*axial_to_xy(*p["cells"][0]),
                f"{pid}\n{p['anchor']}/{p['rod']}",
                ha="center", va="center", fontsize=7, color="white", weight="bold")

    # bounds
    xs, ys = zip(*(axial_to_xy(q,r)
                   for q in range(m.q0, m.q1+1)
                   for r in range(m.r0, m.r1+1)))
    pad = HEX_R*2.0
    ax.set_xlim(min(xs)-pad, max(xs)+pad)
    ax.set_ylim(min(ys)-pad, max(ys)+pad)

    fig.savefig(os.path.join(outdir, f"frame_{step:04d}.png"),
                bbox_inches="tight")
    plt.close(fig)

# ---------- MAIN ----------
def main():
    # initial rhombus region and PBC enabled
    m = Membrane(q0=0, q1=16, r0=0, r1=12, use_pbc=True, seed=42)

    # “case list” + probabilities
    CASES = [event_swap_lipids, event_insert_peptide, event_add_lipid, event_grow_if_saturated]
    EVENT_WEIGHTS = [0.80, 0.15, 0.04, 0.01]  # swap, insert peptide, add lipid, growth check

    outdir      = "frames_case_engine"
    N_EVENTS    = 600
    SNAP_EVERY  = 25

    log = []
    save_frame(m, 0, outdir)
    for t in range(1, N_EVENTS+1):
        which = weighted_choice(EVENT_WEIGHTS)
        if CASES[which] is event_swap_lipids:
            msg = CASES[which](m, n_swaps=250)
        elif CASES[which] is event_insert_peptide:
            msg = CASES[which](m, tries=400)
        elif CASES[which] is event_add_lipid:
            msg = CASES[which](m, tries=200)
        else:
            msg = CASES[which](m)

        log.append((t, which, msg))
        if (t % SNAP_EVERY) == 0:
            save_frame(m, t, outdir)

    with open(os.path.join(outdir, "event_log.json"), "w") as f:
        json.dump(log, f, indent=2)

if __name__ == "__main__":
    main()

