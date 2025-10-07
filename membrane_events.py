# membrane_events.py  (extended w/ flip + probabilistic adhesion)
from __future__ import annotations
import json, random, math
from dataclasses import dataclass, field
from typing import Dict, Tuple, Optional, List, Set

# ================= Lattice (axial, pointy-top) ====================
Axial = Tuple[int, int]
DIRS: Tuple[Axial, ...] = ((1,0),(1,-1),(0,-1),(-1,0),(-1,1),(0,1))

def add(a: Axial, b: Axial) -> Axial:
    return (a[0]+b[0], a[1]+b[1])

def neighbors(h: Axial) -> Tuple[Axial, ...]:
    return tuple(add(h, d) for d in DIRS)

def is_neighbor(a: Axial, b: Axial) -> bool:
    dq, dr = b[0]-a[0], b[1]-a[1]
    return (dq, dr) in DIRS

def _rng(seed: Optional[int] = None) -> random.Random:
    return random.Random(seed)

# ===================== Data structures ============================
@dataclass
class Peptide:
    pid: str
    triad: Tuple[Axial, Axial, Axial]      # (center, n_k, n_{k+1})
    anchor_len: int
    tail_len: int
    side: str = "outside"                  # "inside" | "outside"
    chirality: str = "L"                   # "L" | "D"
    orientation_k: int = 0                 # k selects the two neighbors

@dataclass
class State:
    amphiphiles: Set[Axial]
    amph_meta: Dict[Axial, Dict[str, int]]
    peptides: Dict[str, Peptide] = field(default_factory=dict)
    displaced_count: int = 0
    n: int = 0                             # rhombus size (q,r ∈ [0..n-1])
    bonds: Set[Tuple[str,str]] = field(default_factory=set)  # persistent peptide-peptide bonds

    def to_json(self) -> str:
        return json.dumps({
            "amphiphiles": list(self.amphiphiles),
            "amph_meta": {f"{q},{r}": v for (q,r), v in self.amph_meta.items()},
            "peptides": {
                pid: {
                    "triad": list(p.triad),
                    "anchor_len": p.anchor_len,
                    "tail_len": p.tail_len,
                    "side": p.side,
                    "chirality": p.chirality,
                    "orientation_k": p.orientation_k,
                } for pid, p in self.peptides.items()
            },
            "displaced_count": self.displaced_count,
            "n": self.n,
            "bonds": list(sorted(tuple(sorted(x)) for x in self.bonds)),
        })

    @staticmethod
    def from_json(s: str) -> "State":
        o = json.loads(s)
        peptides: Dict[str, Peptide] = {}
        for pid, d in o["peptides"].items():
            triad = tuple(tuple(x) for x in d["triad"])  # type: ignore
            peptides[pid] = Peptide(
                pid=pid,
                triad=triad,
                anchor_len=int(d["anchor_len"]),
                tail_len=int(d["tail_len"]),
                side=str(d["side"]),
                chirality=str(d.get("chirality","L")),
                orientation_k=int(d.get("orientation_k",0)),
            )
        amphiphiles = set(tuple(x) for x in o["amphiphiles"])
        amph_meta: Dict[Axial, Dict[str,int]] = {}
        for k, v in o["amph_meta"].items():
            q, r = map(int, k.split(","))
            amph_meta[(q,r)] = {kk:int(vv) for kk, vv in v.items()}
        bonds = set(tuple(x) for x in o.get("bonds", []))
        return State(
            amphiphiles=amphiphiles,
            amph_meta=amph_meta,
            peptides=peptides,
            displaced_count=int(o.get("displaced_count",0)),
            n=int(o.get("n",0)),
            bonds=bonds,
        )

# ===================== Initialization helpers =====================
def _rand_C(rr: random.Random) -> int:
    return rr.randint(8, 20)  # fatty-acid carbon count

def init_state(n: int, seed: Optional[int] = None) -> State:
    rr = _rng(seed)
    amph: Set[Axial] = set()
    meta: Dict[Axial, Dict[str,int]] = {}
    for q in range(n):
        for r in range(n):
            h = (q,r)
            amph.add(h)
            meta[h] = {"C": _rand_C(rr)}
    return State(amphiphiles=amph, amph_meta=meta, n=n)

# ===================== Triad geometry =============================
def triad_from_center_and_k(center: Axial, k: int) -> Tuple[Axial, Axial, Axial]:
    n1 = add(center, DIRS[k % 6])
    n2 = add(center, DIRS[(k+1) % 6])
    return (center, n1, n2)

def random_triad_and_k(center: Axial, rr: random.Random) -> Tuple[Tuple[Axial,Axial,Axial], int]:
    k = rr.randrange(6)
    return triad_from_center_and_k(center, k), k

def triad_valid_for_insert(st: State, tri: Tuple[Axial,Axial,Axial]) -> bool:
    occupied = set().union(*(set(p.triad) for p in st.peptides.values()))
    if any(h in occupied for h in tri):
        return False
    return all(h in st.amphiphiles for h in tri)

# ============================ Events ==============================
def insert_peptide(st: State,
                   pid: str,
                   anchor_len: int,
                   tail_len: int,
                   side: str = "outside",
                   chirality: str = "L",
                   seed: Optional[int] = None,
                   max_tries: int = 5000) -> bool:
    rr = _rng(seed)
    candidates = list(st.amphiphiles)
    rr.shuffle(candidates)

    for _ in range(max_tries):
        if not candidates:
            candidates = list(st.amphiphiles)
            rr.shuffle(candidates)
            if not candidates:
                return False
        c0 = candidates.pop()
        tri, k = random_triad_and_k(c0, rr)
        if triad_valid_for_insert(st, tri):
            for h in tri:
                st.amphiphiles.discard(h)
                st.amph_meta.pop(h, None)
            st.displaced_count += 3
            st.peptides[pid] = Peptide(
                pid=pid, triad=tri, anchor_len=anchor_len, tail_len=tail_len,
                side=side, chirality=chirality, orientation_k=k
            )
            return True
    return False

def rotate_peptide(st: State, pid: str, direction: int = +1, seed: Optional[int] = None) -> bool:
    if pid not in st.peptides:
        return False
    rr = _rng(seed)
    p = st.peptides[pid]
    center, *_ = p.triad
    old_cells = set(p.triad)
    new_k = (p.orientation_k + (1 if direction >= 0 else -1)) % 6
    new_triad = triad_from_center_and_k(center, new_k)
    new_cells = set(new_triad)

    occupied_other = set().union(*(set(pp.triad) for kpp, pp in st.peptides.items() if kpp != pid))
    for h in new_cells:
        if h in occupied_other:
            return False
        if (h not in old_cells) and (h not in st.amphiphiles):
            return False

    leaving = old_cells - new_cells
    joining = new_cells - old_cells
    for h in leaving:
        st.amphiphiles.add(h)
        st.amph_meta[h] = {"C": _rand_C(rr)}
    for h in joining:
        if h in st.amphiphiles:
            st.amphiphiles.discard(h)
            st.amph_meta.pop(h, None)
            st.displaced_count += 1
    p.triad = new_triad
    p.orientation_k = new_k
    return True

def flip_peptide_side(st: State, pid: str, seed: Optional[int] = None) -> bool:
    """
    Toggle inside/outside and try a 180° triad pivot (k -> k+3 mod 6).
    Only applies the pivot if legal; side flips regardless (so we can log biochemical state).
    """
    if pid not in st.peptides:
        return False
    p = st.peptides[pid]
    # toggle side label first
    p.side = "inside" if p.side == "outside" else "outside"

    center, *_ = p.triad
    old_cells = set(p.triad)
    new_k = (p.orientation_k + 3) % 6   # 180° around center
    new_triad = triad_from_center_and_k(center, new_k)
    new_cells = set(new_triad)

    occupied_other = set().union(*(set(pp.triad) for kpp, pp in st.peptides.items() if kpp != pid))
    # check feasibility
    ok = True
    for h in new_cells:
        if h in occupied_other:
            ok = False; break
        if (h not in old_cells) and (h not in st.amphiphiles):
            ok = False; break
    if not ok:
        return True  # side flipped, geometry unchanged

    # apply pivot
    rr = _rng(seed)
    leaving = old_cells - new_cells
    joining = new_cells - old_cells
    for h in leaving:
        st.amphiphiles.add(h)
        st.amph_meta[h] = {"C": _rand_C(rr)}
    for h in joining:
        if h in st.amphiphiles:
            st.amphiphiles.discard(h)
            st.amph_meta.pop(h, None)
            st.displaced_count += 1
    p.triad = new_triad
    p.orientation_k = new_k
    return True

def swap_amphiphiles(st: State, attempts: int = 1, seed: Optional[int] = None) -> int:
    rr = _rng(seed)
    amph_list = list(st.amphiphiles)
    if not amph_list:
        return 0
    swaps = 0
    for _ in range(attempts):
        a = rr.choice(amph_list)
        nbrs = [n for n in neighbors(a) if n in st.amphiphiles]
        if not nbrs:
            continue
        b = rr.choice(nbrs)
        st.amph_meta[a], st.amph_meta[b] = st.amph_meta.get(b, {}), st.amph_meta.get(a, {})
        swaps += 1
    return swaps

def grow_if_threshold(st: State) -> bool:
    threshold = 4*st.n + 2
    if st.displaced_count < threshold:
        return False
    old_n = st.n
    new_n = old_n + 2
    st.n = new_n

    rr = _rng(None)
    peptide_cells = set().union(*(set(p.triad) for p in st.peptides.values()))
    for q in range(new_n):
        for r in range(new_n):
            h = (q,r)
            if h not in st.amphiphiles and h not in peptide_cells:
                st.amphiphiles.add(h)
                st.amph_meta[h] = {"C": _rand_C(rr)}
    st.displaced_count = 0
    return True

# ====================== Rafts & adhesion ==========================
def peptides_touch(p1: Peptide, p2: Peptide) -> bool:
    for a in p1.triad:
        for b in p2.triad:
            if is_neighbor(a, b):
                return True
    return False

def _canon_edge(a: str, b: str) -> Tuple[str,str]:
    return (a,b) if a < b else (b,a)

def probabilistic_adhesion(st: State,
                           p_stick_base: float = 0.2,
                           same_chirality_bonus: float = 0.3,
                           seed: Optional[int] = None) -> int:
    """
    For each touching peptide pair, form a persistent bond with probability
    p = p_stick_base + (same_chirality_bonus if same chirality else 0).
    Returns number of new bonds formed.
    """
    rr = _rng(seed)
    pids = list(st.peptides.keys())
    new_bonds = 0
    for i in range(len(pids)):
        for j in range(i+1, len(pids)):
            a, b = pids[i], pids[j]
            pa, pb = st.peptides[a], st.peptides[b]
            if not peptides_touch(pa, pb):
                continue
            e = _canon_edge(a,b)
            if e in st.bonds:
                continue
            p = p_stick_base + (same_chirality_bonus if pa.chirality == pb.chirality else 0.0)
            p = max(0.0, min(1.0, p))
            if rr.random() < p:
                st.bonds.add(e)
                new_bonds += 1
    return new_bonds

def label_rafts_from_bonds(st: State) -> Dict[str,int]:
    """
    Rafts are connected components of the bond graph (persistent adhesion).
    If no bonds exist, isolated peptides are each their own raft.
    """
    pids = list(st.peptides.keys())
    adj: Dict[str, List[str]] = {pid: [] for pid in pids}
    for a,b in st.bonds:
        adj[a].append(b)
        adj[b].append(a)

    raft_map: Dict[str,int] = {}
    raft_id = 0
    for pid in pids:
        if pid in raft_map:
            continue
        raft_id += 1
        stack = [pid]
        raft_map[pid] = raft_id
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v not in raft_map:
                    raft_map[v] = raft_id
                    stack.append(v)
    return raft_map

# ================= JSON wrappers (for Mathematica) ================
def init_state_json(n: int, seed: Optional[int] = None) -> str:
    return init_state(n, seed).to_json()

def insert_peptide_json(state_json: str, pid: str, anchor_len: int, tail_len: int,
                        side: str = "outside", chirality: str = "L", seed: Optional[int] = None) -> str:
    st = State.from_json(state_json)
    ok = insert_peptide(st, pid, anchor_len, tail_len, side, chirality, seed)
    return json.dumps({"ok": ok, "state": json.loads(st.to_json())})

def rotate_peptide_json(state_json: str, pid: str, direction: int = +1, seed: Optional[int] = None) -> str:
    st = State.from_json(state_json)
    ok = rotate_peptide(st, pid, direction, seed)
    return json.dumps({"ok": ok, "state": json.loads(st.to_json())})

def flip_peptide_side_json(state_json: str, pid: str, seed: Optional[int] = None) -> str:
    st = State.from_json(state_json)
    ok = flip_peptide_side(st, pid, seed)
    return json.dumps({"ok": ok, "state": json.loads(st.to_json())})

def swap_amphiphiles_json(state_json: str, attempts: int = 1, seed: Optional[int] = None) -> str:
    st = State.from_json(state_json)
    swaps = swap_amphiphiles(st, attempts, seed)
    return json.dumps({"swaps": swaps, "state": json.loads(st.to_json())})

def grow_if_threshold_json(state_json: str) -> str:
    st = State.from_json(state_json)
    grew = grow_if_threshold(st)
    return json.dumps({"grew": grew, "state": json.loads(st.to_json())})

def probabilistic_adhesion_json(state_json: str,
                                p_stick_base: float = 0.2,
                                same_chirality_bonus: float = 0.3,
                                seed: Optional[int] = None) -> str:
    st = State.from_json(state_json)
    nb = probabilistic_adhesion(st, p_stick_base, same_chirality_bonus, seed)
    rafts = label_rafts_from_bonds(st)
    return json.dumps({"new_bonds": nb, "rafts": rafts, "state": json.loads(st.to_json())})

# ==================== Headless figure saver =======================
def save_frame_png(st: State, path: str, dpi: int = 180) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import RegularPolygon

    size = 0.5
    w = math.sqrt(3) * size
    h = 1.5 * size

    def axial_to_xy(q: int, r: int) -> Tuple[float, float]:
        x = w * (q + r/2.0)
        y = h * r
        return x, y

    fig, ax = plt.subplots(figsize=(9, 10))
    ax.set_aspect('equal')
    ax.axis("off")

    # Amphiphiles
    for a in st.amphiphiles:
        x, y = axial_to_xy(*a)
        patch = RegularPolygon((x, y), numVertices=6, radius=size, orientation=math.radians(30),
                               facecolor="#f2f2f2", edgecolor="#444", linewidth=0.75)
        ax.add_patch(patch)
        C = st.amph_meta.get(a, {}).get("C", None)
        if C is not None:
            ax.text(x, y, f"{C}", ha="center", va="center", fontsize=6, color="#333")

    # Peptides
    colors = ("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
              "#e377c2","#7f7f7f","#bcbd22","#17becf")
    centers: Dict[str, Tuple[float,float]] = {}
    for i, p in enumerate(st.peptides.values()):
        col = colors[i % len(colors)]
        for hcell in p.triad:
            x, y = axial_to_xy(*hcell)
            patch = RegularPolygon((x, y), numVertices=6, radius=size, orientation=math.radians(30),
                                   facecolor=col, edgecolor="#222", linewidth=1.1)
            ax.add_patch(patch)
        cx = sum(axial_to_xy(*h)[0] for h in p.triad)/3.0
        cy = sum(axial_to_xy(*h)[1] for h in p.triad)/3.0
        centers[p.pid] = (cx, cy)
        ax.text(cx, cy, f"{p.pid}\n{p.anchor_len}/{p.tail_len}\n{p.side}\n{p.chirality}",
                ha="center", va="center", fontsize=7, color="white")

    # Draw bonds (rafts)
    for a,b in st.bonds:
        if a in centers and b in centers:
            (x1,y1), (x2,y2) = centers[a], centers[b]
            ax.plot([x1,x2],[y1,y2], linewidth=2.0, alpha=0.7, color="#111")

    # Tight bounds
    all_cells = list(st.amphiphiles) + [h for p in st.peptides.values() for h in p.triad]
    if all_cells:
        xs, ys = zip(*(axial_to_xy(*h) for h in all_cells))
        pad = 2*size
        ax.set_xlim(min(xs)-pad, max(xs)+pad)
        ax.set_ylim(min(ys)-pad, max(ys)+pad)

    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

# ======================= CLI smoke-test ===========================
if __name__ == "__main__":
    st = init_state(28, seed=1)
    # Insert peptides
    for k in range(1, 9):
        insert_peptide(st, pid=f"P{k}",
                       anchor_len=random.randint(2,3),
                       tail_len=random.randint(8,14),
                       side=("outside" if k%2 else "inside"),
                       chirality=("L" if k%2 else "D"))
    # Rotate and flip some
    rotate_peptide(st, "P1", +1)
    rotate_peptide(st, "P2", -1)
    flip_peptide_side(st, "P2")
    # Amphiphile mixing and growth
    swap_amphiphiles(st, attempts=2500)
    grow_if_threshold(st)
    # Probabilistic adhesion
    probabilistic_adhesion(st, p_stick_base=0.25, same_chirality_bonus=0.35)
    # Save figure
    save_frame_png(st, "frame_test.png", dpi=180)
    print("Saved frame_test.png")

