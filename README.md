# DFlow — Peptide–Membrane Simulation (Reversible, Diurnal-Driven)

**DFlow** models a 2D amphiphile membrane on a hexagonal lattice with embedded peptide “rafts.”  
All **membrane** processes respect **microscopic reversibility** (forward ⇄ reverse). The sole non-equilibrium driver is **diurnal-biased peptide polymerization** in solution (wet/dry or hot/cold cycling). Peptide **insertion** into the membrane is **hydrophobic-mismatch gated**.

This repo contains a single-file runner (`dflow_reversible.py`) plus modular sources under `src/` (WIP).

---

## 🧠 Scientific Principles (this build)

- **Microscopic reversibility:** Every allowed membrane event has a reverse (or is self-inverse).  
- **Single external bias:** Only **polymerization in solution** is biased by day/night; depolymerization remains possible (non-zero).  
- **Hydrophobic mismatch gating (window + soft acceptance):**  
  - Hard window: insert if \( |L_p - d_{\text{site}}| \le \varepsilon \).  
  - Soft gate: slight over-mismatch may insert with Boltzmann weight \( \exp(-\beta\,\Delta) \).  
- **Helix requirement & length cap:** Peptides are generated with a **helical/coil** flag. Only **helical** peptides can insert. Length is **capped at 12 aa** (Deamer note; primary reference pending).  
- **Crowding feedback:** Peptides that don’t fit the membrane go to a **crowding pool** which **slows polymerization** (PCR-like) via \( f=1/(1+\gamma\,C) \).  
- **Rafts & 2D Stokes-like diffusion:** Peptides can form **reversible bonds** when touching; raft diffusion attempt rate scales as \( D \sim D_0/\text{size}^\alpha \).  
- **Optional fusion exception (off by default):** You can toggle **quasi-irreversible fusion** for large rafts (no dissociation above a size threshold).

---

## 📦 Layout

```

.
├─ dflow_reversible.py        # main physically-grounded runner (use this)
├─ dflow_single.py            # minimal demo (optional)
├─ src/
│  ├─ membrane_model.py       # reusable helpers (WIP)
│  └─ peptide_generator.py    # reusable helpers (WIP)
├─ examples/
│  └─ legacy/
│     └─ membrane_events_case_loop.py
├─ configs/                   # optional small configs
├─ requirements.txt
└─ README.md

````

> Outputs (frames, logs) go under `runs/...` and should be git-ignored.

---

## 🛠️ Install

```bash
python -V              # Python 3.8+ (3.10+ recommended)
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
````

`requirements.txt`

```
numpy>=1.20
matplotlib>=3.4
```

---

## 🚀 Quick Start

Run the reversible, diurnal-driven simulation:

```bash
python dflow_reversible.py \
  --N0 12 \
  --TOTAL_EVENTS 1500 \
  --OUT runs/exp_phys \
  --SAVE_STEPS 0 750 1500
```

Artifacts:

* `runs/exp_phys/frames/frame_0000.png` (+ frames at save steps)
* `runs/exp_phys/event_log.json`
* `runs/exp_phys/raft_size_histogram.png`
* `runs/exp_phys/raft_sizes.csv`
* `runs/exp_phys/peptide_pool.json`
* `runs/exp_phys/crowding.json`

Minimal demo:

```bash
python dflow_single.py --N0 12 --TOTAL_EVENTS 1000 --OUT runs/exp_0001 --SAVE_STEPS 0 500 1000
```

---

## ⚙️ Key CLI Options

`dflow_reversible.py`

* `--N0` *(int)*: half-width of initial rhombus (lattice size).
* `--TOTAL_EVENTS` *(int)*: number of stochastic events.
* `--OUT` *(path)*: output directory.
* `--SAVE_STEPS` *(ints)*: steps at which frames are saved (e.g., `0 750 1500`).
* **Diurnal driver:**

  * `--DAY_STEPS`, `--NIGHT_STEPS` *(int)*: steps per phase.
  * `--POLY_GAIN_DAY`, `--POLY_GAIN_NIGHT` *(float)*: polymerization bias multipliers.
* `--BETA` *(float)*: inverse temperature used in the **soft mismatch gate**.

Example (strong day bias, shorter days):

```bash
python dflow_reversible.py --N0 14 --TOTAL_EVENTS 2000 \
  --DAY_STEPS 6 --NIGHT_STEPS 6 --POLY_GAIN_DAY 12.0 --POLY_GAIN_NIGHT 0.15 \
  --OUT runs/day_short --SAVE_STEPS 0 1000 2000
```

---

## 🧩 Events (forward ⇄ reverse)

**Amphiphiles / membrane**

* `A1_swap` — lateral carbon-length swap (self-inverse).
* `A3_plus_thicken` / `A3_minus_thin` — site thickness 2 nm ⇄ 4 nm (mono⇄di).

**Peptides (solution & membrane)**

* `P1_plus_polymerize` / `P1_minus_depolymerize` — **biased** polymerization (diurnal), reverse allowed. Crowding slows the forward rate.
* `P2_plus_insert` / `P2_minus_desorb` — **mismatch-gated** insertion/ejection (hard window + soft acceptance; helix required).
* `P3_step` — diffusion; **raft-scaled attempt rate** (2D Stokes-like).
* `P4_flip` — inside↔outside (metadata only).

**Raft association**

* `R1_plus_assoc` / `R1_minus_dissoc` — reversible bonds between edge-touching peptides.
* *(Optional)* **Fusion exception**: when enabled, bonds in large rafts don’t dissociate.

---

## 🔧 Parameters (in code)

* **Mismatch window:** `EPS_MISMATCH_NM = 0.5` (hard), with soft gate fraction `SOFT_GATE_FRACTION = 0.6` and `BETA = 1.0`.
* **Peptide length cap:** `MAX_PEPT_LEN = 12`; helical probability ~0.7 (tunable).
* **Crowding feedback:** `CROWDING_GAMMA = 0.02` (slowdown), `CROWDING_DECAY = 0.0` (no decay by default).
* **Raft diffusion scaling:** `RAFT_D0 = 1.0`, `RAFT_DIFF_SIZE_EXP = 1.0`.
* **Fusion toggle:** `FUSION_IRREVERSIBLE = False`, `FUSION_SIZE_THRESH = 6`.

> All of the above are defined near the top of `dflow_reversible.py` for easy tuning.

---

## 🧱 Data Model (key fields)

* `amph[(q,r)] = { "carbon": int, "mono_di": 1|2, "pep": bool }`
* `peptides[pid] = { "cent": (q,r), "orient": int, "inside": bool, "chir": "L"|"D"|"0", "length": int, "Lp_nm": float, "helical": bool }`
* `bonds = { frozenset({pid1,pid2}), ... }`
* `pool = {"L": int, "D": int, "0": int}`
* `crowding_count` *(int)*

---

## 📊 Outputs

* **Frames**: membrane snapshots, peptides colored by **bond-connected components** (rafts).
* **Histogram**: `raft_size_histogram.png` (component sizes).
* **Logs**: `event_log.json` (all events), `peptide_pool.json`, `crowding.json`, `raft_sizes.csv`.

---

## 🤝 Contributing

1. Create a feature branch: `feat/event-<name>`.
2. Add your event as a **reversible pair** (or self-inverse) and register it in the scheduler.
3. Preserve insertion **mismatch gating** and mass accounting (pool ↔ membrane ↔ crowding).
4. Don’t commit outputs/binaries; rely on `.gitignore`.

**Good starter ideas:**

* Temperature-dependent rates (affect both forward/reverse symmetrically).
* Curvature/protein coupling proxies (keep reverse path available).
* Alternative crowding models or compartmentalization.

---

## 🧪 Quick Checks (manual)

* Insertions are mostly **helical** and within the **mismatch window**.
* Rejected insertions increase `crowding_count`; polymerization rate drops as crowding grows.
* Raft histogram forms; larger rafts **diffuse less** (fewer `P3_step` on big clusters).
* Optional: enable fusion toggle and verify large rafts stop dissociating.

---

## 🐞 Troubleshooting

* **TypeError with `str|None`** → upgrade to Python ≥ 3.10 or use this repo (already 3.8-compatible).
* **No frames** → ensure `--SAVE_STEPS` includes `0` or desired steps; check `runs/.../frames`.
* **Backend warnings** → `matplotlib` uses `Agg` for headless saves; this is expected.

---

## 📚 Notes & Acknowledgments

* Peptide length cap (~12 aa) follows a note in David Deamer’s *First Life* (Kindle p.84). We’re seeking a primary source; if you have one, please share.
* Project contributors include Richard Gordon, Syed Hussain Ather, Savino Longo, and collaborators.

*If you use this code, please cite the DFlow project and related manuscripts in preparation.*

