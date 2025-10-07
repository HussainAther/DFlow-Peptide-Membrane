# DFlow — Peptide–Membrane Simulation (Reversible, Diurnal-Driven)

**DFlow** models a 2D amphiphile membrane on a hex lattice with embedded peptide “rafts.”  
All membrane processes obey **microscopic reversibility** (forward ⇄ reverse), with a **single non-equilibrium driver**: diurnal-biased **peptide polymerization** (wet/dry or hot/cold cycling as the Earth rotates).  
Peptide **insertion/ejection** is **hydrophobic-mismatch gated** against local membrane thickness.

---

## 🧠 Scientific principles (what makes this “physical”)

- **Microscopic reversibility:** For every allowed membrane event, the reverse event exists and is accessible.  
- **Single bias:** Only **polymerization in solution** is externally biased by the day–night cycle; depolymerization remains possible (non-zero).  
- **Hydrophobic mismatch gating:** A peptide’s effective hydrophobic length \(L_p\) must match site thickness (2 nm for mono-carboxylic, 4 nm for di-carboxylic) within a tolerance \(\varepsilon\) to insert.  
- **Rafts as bonds:** Peptides form **reversible association bonds** (edges touching). Rafts grow or split naturally as bonds associate/dissociate.  
- **Mass accounting:** A **solution pool** stores peptides. Polymerization adds to the pool; insertion consumes; desorption returns.

---

## 📦 Repo layout

```

.
├─ dflow_reversible.py        # main physically-grounded runner (use this)
├─ dflow_single.py            # minimal baseline runner (simple demo)
├─ src/
│  ├─ membrane_model.py       # reusable helpers (optional, WIP)
│  └─ peptide_generator.py    # reusable helpers (optional, WIP)
├─ examples/
│  └─ legacy/
│     └─ membrane_events_case_loop.py  # earlier prototype (optional)
├─ configs/                   # (optional) small text/yaml configs
├─ requirements.txt
└─ README.md

````

> Outputs (frames, logs) are written under `runs/…` and should be git-ignored.

---

## 🛠️ Install

```bash
python -V              # Python 3.8+ (3.10+ recommended)
python -m venv .venv && source .venv/bin/activate  # or your preferred env
pip install -r requirements.txt
````

`requirements.txt`

```
numpy>=1.20
matplotlib>=3.4
```

---

## 🚀 Quick start

Run the reversible, physically grounded simulation:

```bash
python dflow_reversible.py \
  --N0 12 \
  --TOTAL_EVENTS 1500 \
  --OUT runs/exp_phys \
  --SAVE_STEPS 0 750 1500
```

Artifacts (created automatically):

* `runs/exp_phys/frames/frame_0000.png` (and other frames)
* `runs/exp_phys/event_log.json`
* `runs/exp_phys/raft_size_histogram.png`
* `runs/exp_phys/raft_sizes.csv`
* `runs/exp_phys/peptide_pool.json`

Minimal demo (simpler visuals, fewer physics):

```bash
python dflow_single.py --N0 12 --TOTAL_EVENTS 1000 --OUT runs/exp_0001 --SAVE_STEPS 0 500 1000
```

---

## ⚙️ CLI options (key ones)

`dflow_reversible.py`

* `--N0` *(int)*: half-width of rhombus region (initial lattice size).
* `--TOTAL_EVENTS` *(int)*: number of stochastic events to apply.
* `--OUT` *(path)*: output directory for frames + logs.
* `--SAVE_STEPS` *(ints)*: steps at which to save frames (e.g., `0 750 1500`).
* **Diurnal driver:**

  * `--DAY_STEPS`, `--NIGHT_STEPS` *(int)*: steps per phase.
  * `--POLY_GAIN_DAY`, `--POLY_GAIN_NIGHT` *(float)*: bias multipliers for polymerization.
* `--BETA` *(float)*: inverse temperature (β = 1/kT) placeholder for rate helpers.

Example (strong day bias, short days):

```bash
python dflow_reversible.py --N0 14 --TOTAL_EVENTS 2000 \
  --DAY_STEPS 6 --NIGHT_STEPS 6 --POLY_GAIN_DAY 12.0 --POLY_GAIN_NIGHT 0.15 \
  --OUT runs/day_short --SAVE_STEPS 0 1000 2000
```

---

## 🧩 Events (forward ⇄ reverse)

**Amphiphiles & thickness**

* `A1_swap` — lateral carbon-length swap (self-inverse).
* `A3_plus_thicken` / `A3_minus_thin` — mono⇄di thickness (2 nm ⇄ 4 nm).

**Peptides (solution & membrane)**

* `P1_plus_polymerize` / `P1_minus_depolymerize` — **biased** polymerization; reverse stays > 0.
* `P2_plus_insert` / `P2_minus_desorb` — **mismatch-gated** insertion/ejection.
* `P3_step` — triad diffusion (self-inverse if symmetric).
* `P4_flip` — inside↔outside orientation flip (metadata).

**Raft association**

* `R1_plus_assoc` / `R1_minus_dissoc` — reversible bonds between edge-touching peptides.

---

## 🧱 Data model (key fields)

* `amph[(q,r)] = { "carbon": int, "mono_di": 1|2, "pep": bool }`
* `peptides[pid] = { "cent": (q,r), "orient": int, "inside": bool, "chir": "L"|"D"|"0", "length": int, "Lp_nm": float }`
* `bonds = { frozenset({pid1,pid2}), ... }`
* `pool = {"L": int, "D": int, "0": int}`  (solution peptides)

---

## 📊 Outputs

* **Frames**: membrane snapshots with peptides colored by **bond-connected components** (rafts).
* **Histogram**: `raft_size_histogram.png` (component sizes).
* **Logs**: `event_log.json` (sequence of events), `peptide_pool.json`, `raft_sizes.csv`.

---

## 🤝 Contributing

1. Open a branch: `feat/event-pept-desorb` (or similar).
2. Add your event as a **reversible pair** (or self-inverse) and wire into the scheduler.
3. Keep insertion rules **mismatch-gated**.
4. Don’t commit outputs/binaries; rely on `.gitignore`.

Recommended new events to try:

* `pept_desorb` (done) — refine kinetics/energetics.
* `raft_fission` — when weakly bonded clusters split (via dissociation rules).
* `diurnal_modulators` — let day/night modulate additional **rates**, but keep forward/reverse symmetry intact.

---

## 🧪 Tests (optional starter)

* Run a short sim and assert that:

  * `P2_plus_insert` increases covered sites; `P2_minus_desorb` restores them.
  * Bonds appear/disappear over time (`R1_plus_assoc`/`R1_minus_dissoc`).
  * Pool counts move in expected directions with polymerization bias.

(We’ll add a tiny `pytest` suite later.)

---

## 🐞 Troubleshooting

* **TypeError with `str|None`**: you’re on Python < 3.10. Use this repo (3.8-compatible) or upgrade Python.
* **No frames appear**: check `--SAVE_STEPS` includes `0` or your chosen steps; confirm `runs/.../frames/` exists.
* **Matplotlib backend**: we force `Agg` for headless save; this is normal on servers.

---

## 📚 Citation / Attribution

If you use this in a paper or report, please cite the DFlow project and the associated manuscripts in preparation (Longo et al.; Gordon & Ather, 2025).
Contributors: R. Gordon, S. H. Ather, S. Longo, et al.

