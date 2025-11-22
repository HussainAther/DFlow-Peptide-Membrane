# tests/test_r2.py
import types
import numpy as np
import random
import pytest

import dflow_reversible as sim   

# -------------------------
# Helpers
# -------------------------

def mk_cfg(**over):
    args = types.SimpleNamespace(
        SEED=123,
        N0=4,
        HEX_RADIUS=1.0,
        TOTAL_EVENTS=20,
        OUT="runs/_tests",
        SAVE_STEPS=[0],
        DAY_STEPS=10,
        NIGHT_STEPS=10,
        POLY_GAIN_DAY=0.0,
        POLY_GAIN_NIGHT=0.0,
        BETA=0.0,
        MAX_PEPT_LEN=8,
        A_INSERT_MISMATCH_TOL=10.0,
        CROWDING_GAMMA=0.0,
        CROWDING_DECAY=0.0,
        DESORB_TO_CROWDING=False,
        DISCARD_TO_CROWDING_P=0.0,
        RAFT_D0=1.0,
        RAFT_DIFF_SIZE_EXP=0.0,
        FUSION_IRREVERSIBLE=False,
        FUSION_SIZE_THRESH=999,
        FISSION_MIN_SIZE=2,
        FISSION_HEURISTIC="bridge",
        FISSION_K0=1.0,
        FISSION_SIZE_SCALING="none",
        FISSION_ENABLED=True,
        FISSION_ALLOW_PURE_CHIRAL=False,
        FUSION_K0=1.0,
    )
    for k, v in over.items():
        setattr(args, k, v)
    return sim.cfg_from_args(args)

def env(beta=0.0):
    return sim.DiurnalEnv(sim.Diurnal(10, 10, 0.0, 0.0), beta=beta, step_idx=0)

def add_tri(mem, cent, orient, raft, chir='L'):
    """Insert a peptide triad at exact coordinates without using pool/insert."""
    pid = mem.pid_counter
    mem.pid_counter += 1
    mem.peptides[pid] = {
        "cent": cent,
        "orient": orient,
        "raft": raft,
        "inside": False,
        "nres": 6,
        "chir": chir,
    }
    for c in sim.rotate_triad(cent, orient):
        mem.amph[c]["pep"] = True
    return pid

def comp_sizes(mem):
    comps = sim.connected_components_full(mem)
    return sorted(len(c) for c in comps)

# -------------------------
# Tests
# -------------------------

def test_fission_bridge_split_preserves_mass():
    """
    Verify:
      - peptide count is conserved
      - raft partition increases (>= 2 distinct raft ids)
      - event is logged
    Use a guaranteed-bridge 2-triad setup with *no overlapping cells* but
    with an *edge-touch* so adjacency has an edge (bridge).
    """
    cfg = mk_cfg(
        FISSION_HEURISTIC="bridge",
        FISSION_ALLOW_PURE_CHIRAL=True,  
    )
    random.seed(cfg.SEED); np.random.seed(cfg.SEED)
    mem = sim.Membrane(cfg.N0, cfg.HEX_RADIUS, cfg)
    r = 0

    add_tri(mem, (0, 0), 0, r)

    add_tri(mem, (2, 1), 3, r)

    e = env()
    before = len(mem.peptides)
    sim.do_fission(mem, e)
    after = len(mem.peptides)

    assert after == before

    raft_ids = {p["raft"] for p in mem.peptides.values()}
    assert len(raft_ids) >= 2

    assert any(ev.get("event") == "R2_plus_fission" for ev in mem.event_log)

def test_no_fission_pure_chiral_when_disabled():
    cfg = mk_cfg(FISSION_ALLOW_PURE_CHIRAL=False)
    mem = sim.Membrane(cfg.N0, cfg.HEX_RADIUS, cfg)
    for x in range(3):
        add_tri(mem, (x, 0), 3, 0, chir='L')
    rate = sim.rate_fission(mem, env())
    assert rate == 0.0

def test_fusion_merges_components_and_bonds():
    cfg = mk_cfg()
    mem = sim.Membrane(cfg.N0, cfg.HEX_RADIUS, cfg)
    r0, r1 = 0, 1
    add_tri(mem, (0, 0), 0, r0)
    add_tri(mem, (1, 0), 3, r0)
    add_tri(mem, (1, -1), 0, r1)
    add_tri(mem, (2, -1), 3, r1)
    e = env()
    sim.do_fusion(mem, e)
    sizes = comp_sizes(mem)
    assert sizes[-1] >= 4
    for i, j in mem.bonds:
        assert mem.peptides[i]["raft"] == mem.peptides[j]["raft"]

def test_detailed_balance_sanity():
    """
    With BETA=0 and symmetric K0 and size scaling=none,
    on a tiny system we should observe both R2+ and R2− under long sampling.
    Start from the same non-overlapping, edge-touching 2-triad raft to
    ensure the first split is available; once split, fusion is available.
    """
    cfg = mk_cfg(
        FISSION_SIZE_SCALING="none",
        FISSION_K0=1.0,
        FUSION_K0=1.0,
        FISSION_HEURISTIC="bridge",
        FISSION_ENABLED=True,
        FISSION_ALLOW_PURE_CHIRAL=True,
        BETA=0.0,
    )
    mem = sim.Membrane(cfg.N0, cfg.HEX_RADIUS, cfg)

    add_tri(mem, (0, 0), 0, 0)
    add_tri(mem, (2, 1), 3, 0)

    e = env()
    sch = sim.Scheduler()
    sch.register(sim.Event("R2_plus_fission", sim.rate_fission, sim.do_fission))
    sch.register(sim.Event("R2_minus_fusion", sim.rate_fusion, sim.do_fusion))

    for t in range(1, 1200):
        e.step_idx = t
        sch.step(mem, e)
        e.diurnal.tick()

    fissions = sum(1 for ev in mem.event_log if ev.get("event") == "R2_plus_fission")
    fusions  = sum(1 for ev in mem.event_log if ev.get("event") == "R2_minus_fusion")
    assert fissions > 0 and fusions > 0

    ratio = max(fissions, fusions) / max(1, min(fissions, fusions))
    assert ratio < 4.0
