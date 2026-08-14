# Rebuild / Vibe-Coding Guide — PDAC coexistence-control framework

This document lets you re-create this project on another machine, either by **cloning the
repo** (fast path) or by **rebuilding it from scratch with an AI coding assistant**
(vibe-coding path). It is written so you can paste sections directly to a coding agent
(Claude Code, Cursor, etc.) as the driving spec.

> **What the program is.** A spatial agent-based model (ABM) of the PDAC
> tumor–myCAF–immune ecosystem plus a set of analyses (phase map, Sobol, Pareto, Monte
> Carlo, observation model, spatial statistics, ablations) that produce the manuscript
> figures. In silico, hypothesis-generating; not a clinical predictor.
>
> **Authoritative model spec:** `docs/manuscript/S3_ODD_protocol.md` (ODD protocol) —
> the single source of truth for entities, state variables, the per-step schedule, and
> every submodel equation. When rebuilding `abm.py`, follow the ODD, not this summary.

---

## 0. Two paths

**Path A — clone and run (minutes).**
```bash
git clone https://github.com/kusi81/pdac-coexistence-control.git
cd pdac-coexistence-control
python -m venv .venv && . .venv/Scripts/activate   # Windows Git Bash; use bin/activate on Linux/mac
pip install -r requirements.txt
pip install anndata scanpy squidpy SALib pyyaml     # optional deps for spatial + Sobol
python pipeline/phase_map_dense.py                   # smoke test one analysis
```

**Path B — rebuild from scratch with an AI agent.** Follow §3–§7 below in order. The
agent implements one module at a time and you verify against the numbers in §8.

---

## 1. Environment

- **Python 3.13** (3.11+ works; 3.13 is what results were produced on).
- Core deps (`requirements.txt`): `numpy scipy pandas matplotlib streamlit python-docx`.
- Analysis deps (install separately): `SALib` (Sobol), `anndata scanpy squidpy`
  (real spatial data), `pyyaml` (validate CITATION.cff).
- **OS note (important):** parallel analyses use `concurrent.futures.ProcessPoolExecutor`.
  On Windows this uses *spawn*, so every parallel script must be run as a file
  (`python pipeline/x.py`), never piped via stdin, and must guard work under
  `if __name__ == "__main__":`. Cap workers at ~6 (`WORKERS = 6`).
- **Determinism:** all randomness flows through a single `numpy.random.default_rng(seed)`.
  Same seed + same library versions ⇒ identical results. Do **not** use `Math.random`,
  global `np.random`, or wall-clock seeding.
- **Fonts/encoding:** `pipeline/fonts.py` sets a font that renders a Korean-capable face
  and `matplotlib.rcParams["axes.unicode_minus"] = False`. Use ASCII `-` in labels, not
  U+2212. Reconfigure stdout to UTF-8 at the top of every script
  (`sys.stdout.reconfigure(encoding="utf-8")`).

---

## 2. Module map (what to build, in dependency order)

| Module | Role | Depends on |
|---|---|---|
| `pipeline/fonts.py` | matplotlib font + unicode-minus setup | matplotlib |
| `pipeline/synthetic.py` | `make_tissue(mode, seed)` → contained / diffuse tissue (identical counts) | numpy |
| `pipeline/abm.py` | **core ABM**: `DEFAULT_PARAMS`, `SUBSTANCES`, `TOXICITY`, `TumorABM`, `simulate()`, `control_metrics()` | numpy, scipy.spatial.cKDTree |
| `pipeline/spatial_core.py` | `barrier_score()`, `rim_enrichment()` (matched-null spatial metrics) | numpy, scipy |
| `pipeline/analysis.py` | `run_rim_panel()` (per-sample rim z by cell type) | spatial_core |
| — analyses below all import abm/synthetic — | | |
| `phase_map_dense.py` | Fig 4: 6×6 confinement×exclusion grid, P(keep) uncertainty | abm, synthetic |
| `sobol_analysis.py` | Fig S7: global Sobol (SALib), output = final burden | abm, synthetic, SALib |
| `pareto_seeds.py` | Fig S10: 30-seed progression-constrained Pareto | abm, synthetic |
| `mc_uncertainty.py` | Fig S14: Monte Carlo over compound assumptions | abm, synthetic |
| `obs_model_analysis.py` | Fig S12: CA19-9 observation model | abm, synthetic |
| `gem_schedule.py` | Fig S13: single-agent schedule comparison | abm, synthetic |
| `predict_geometry.py` | Fig S15: geometry vs abundance; invasion front | abm, synthetic |
| `predict_sequence.py` | Fig S16: treatment order | abm, synthetic |
| `predict_biomarker.py` | Fig S17: spatial biomarker predicts immune regime | abm, synthetic, spatial_core |
| `predict_wellmixed.py` | Fig S18: mean-field ablation (ODE, no ABM) | numpy |
| `phase_map_protumor.py` | Fig S11: pro-tumor axis robustness | abm, synthetic |
| `scotia_rim_stats.py` | Fig 3: patient-level CRT stats (needs real data) | analysis, anndata |
| `patient_grounded.py` | Fig S6: run ABM on real CosMx tissue | abm, anndata |
| `build_docx.py` | assemble manuscript.docx from manuscript.md + assets | python-docx |

---

## 3. Build stage 1 — the core ABM (`abm.py`)

**Do this first and get it exactly right; everything else depends on it.** Implement
strictly from `docs/manuscript/S3_ODD_protocol.md`. Key facts the agent must honor:

- **Off-lattice**, 2D field `field_um = 1500`, step `dt_days = 0.5`.
- Cell types: `Tumor` (with a boolean resistance flag), `myCAF`, `iCAF`, `CD8_T`,
  `Macrophage`. Non-tumor cells are never resistant.
- **Global state:** `t, n0, cum_tox, drug_on, phase_label, biomarker, obs_ratio,
  last_obs_t, last_switch_t, low_reads, crit_time`.
- **Per-step schedule (fixed order, method `step()`):** (1) update latent biomarker;
  (2) apply regimen → set `drug_on`, effective params `self.eff`, accrue toxicity;
  (3) tumor proliferation with containment/pressure/drug-block/survival/pro-tumor;
  (4) myCAF activation; (5) myCAF turnover; (6) CD8 migration (barrier-attenuated) +
  killing; (7) tumor apoptosis (drug-induced on sensitive only); (8) CD8 turnover;
  (9) apply removals/additions in a batch; (10) clip to field; (11) record.
- **Update semantics (critical):** division uses the tumor k-d tree built at the *start*
  of the step (daughters this step don't affect same-step divisions = synchronous);
  removals/additions applied together at step end; CD8 move in place then kill.
- **Submodel equations:** copy them verbatim from ODD §3.3 (proliferation probability,
  `eff_cap`, `pen`/`surv`/confinement, activation, turnover, `mob = exp(-α·corr)`,
  killing with `resistant_immune_evasion`, apoptosis, `Poisson(cd8_recruit·dt)`
  recruitment, the CA19-9 observation model, the adaptive band).
- **`DEFAULT_PARAMS`** and their baseline values: see `docs/manuscript/S1_parameters.md`
  (24 params). Notable: `k_prolif=0.11, k_kill=1.0, cd8_recruit=30, cd8_barrier_alpha=0.9,
  caf_confine=0.8, caf_confine_ref=3, caf_pressure=1.2, caf_drug_block=0.6,
  caf_survival=0.0, caf_protumor=0.0, resistance_cost=0.24, resistant_immune_evasion=0.45,
  init_resistant_frac=0.01, mutation_rate=0.001`. Observation-model params:
  `obs_model=False, obs_interval=28, obs_noise_cv=0.25, obs_lag_days=14, min_on_days=14,
  min_off_days=21, obs_confirm=2, obs_safety_mult=1.4, nonsecretor=False`. Runtime cap:
  `max_tumor=0` (off by default).
- **`SUBSTANCES`**: dict of `{name: {label, evidence, effects:{param: multiplier},
  rationale}}`. Effects are multipliers on perturbable params (`k_prolif, k_caf_activate,
  k_kill, cd8_recruit, cd8_speed_um, k_tumor_apoptosis`). Includes the resolved APIs
  (`sac, eupatilin, rg3_20s`) and generics (`generic_cytotoxic, generic_immunostim,
  generic_antifibrotic`). Compose via dose-interpolation `1+(m-1)·dose`, multiply across
  compounds, optional `synergy` lowers `k_prolif`/`k_caf_activate`.
- **`TOXICITY`**: `{name: weight}`; `regimen_toxicity(subs) = Σ weight·dose`.
- **`control_metrics(history, n0, crit_mult=1.5)`** returns `ttp_days,
  progression_censored, peak_frac, final_frac, final_resistant_frac, cum_toxicity,
  auc_burden`. **Do not** add a composite `control_score` — it was removed as flawed;
  rank by `(progression_censored, cum_toxicity)`.
- **`simulate(coords, labels, days, params, regimen_subs=None, schedule=None,
  synergy=0.0, adaptive=False, adapt_on=1.2, adapt_off=0.6, snapshots=(...))`** →
  `(history, snapshots)`.

Sanity check after building: a single 150-day adaptive run of `generic_cytotoxic` at the
escape regime (below) should control the tumor (`final_frac` near 0) at exposure far below
continuous.

---

## 4. Build stage 2 — synthetic tissue (`synthetic.py`)

`make_tissue(mode, seed)` → `(coords[N,2] µm, labels[N], islet_centers)`.
- **`contained`**: tumor islets; myCAF in rings around islets; iCAF distal; CD8 only
  outside the ring (excluded); macrophages in stroma.
- **`diffuse`**: **identical cell counts** but myCAF/iCAF/CD8 scattered (no ring).
- The contained-vs-diffuse pair is the backbone of the geometry ablations (Figs S15/S18),
  so the counts **must** match between modes — verify with a `Counter`.

---

## 5. Build stage 3 — analyses (exact settings)

Each script: subsamples tissue for speed (`CAP`), runs a grid/ensemble in parallel
(`ProcessPoolExecutor`, chunksize ~4), writes a CSV to `data/` and a PNG to `assets/`,
and supports a `replot` mode that redraws from the CSV without re-simulating. Common
"escape regime" baseline used across analyses:
`k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10, init_resistant_frac=0.03,
mutation_rate=0.003, resistant_immune_evasion=0.35`.

| Script | Key settings |
|---|---|
| `phase_map_dense.py` | grid `caf_pressure∈{0,0.6,1.2,1.8,2.4,3.0}` × `cd8_barrier_alpha∈{0.2,0.6,1.0,1.4,1.8,2.4}`; sweep `k_caf_activate∈{0,0.06,0.12,0.20,0.30}`; **5 seeds/cell**; DAYS=90; CAP=2000; drug=`generic_cytotoxic@0.30`. Report **P(keep)** = fraction of seeds where argmin over myCAF is >0 and benefit>0.02. |
| `sobol_analysis.py` | SALib Saltelli, **D=8**, **N=32**, `calc_second_order=False` ⇒ 320 evals; output = **final_frac** (adaptive arm); DAYS=130, CAP=1800, 1 seed, `max_tumor=4500`. Params/bounds: k_prolif[.075,.30], k_kill[.25,1], cd8_recruit[5,20], k_caf_activate[.05,.20], resistance_cost[.10,.40], resistant_immune_evasion[.20,.70], cd8_barrier_alpha[.5,3], caf_confine[.2,1.5]. |
| `pareto_seeds.py` | **30 seeds** (tissue+sim varied together), **synergy=0**, DAYS=150, CAP=2500, adapt 1.1/0.7, `max_tumor=6000`. Rank by `(not censored, exposure)`; report median + 95% interval + per-seed rank + Pareto frontier. |
| `mc_uncertainty.py` | **100 draws**; per draw sample per-compound effect×LogNormal(CV 0.35), exposure weight×LogNormal(CV 0.40), **bioavailability~Beta**: sac(9,1), eupatilin(5,3), rg3_20s(2,6), curcumin(1.5,10), gemcitabine(50,1); synergy~U[0,0.3]. Apply as `m'=1+(m-1)·(effCV·bioav)`; temporarily mutate module `SUBSTANCES`/`TOXICITY` per draw, **restore in finally**. DAYS=120, CAP=1500, `max_tumor=2800`. |
| `obs_model_analysis.py` | arms: continuous, ideal adaptive, observed CA19-9, observed non-secretor; +interval sweep {14,28,56}d; 20 seeds; obs params as §6/ODD. |
| `gem_schedule.py` | one agent (gemcitabine), fixed efficacy; schedules: continuous, intermittent q28 (3wk on/1wk off via `schedule=[{days:21,subs:GEM},{days:7,subs:[]}]`), ideal adaptive, observed CA19-9; 20 seeds. |
| `predict_geometry.py` | contained vs diffuse; modalities {cytotoxic, anti-fibrotic, immune, combos}; +invasion front = 90th-pct tumor radius, untreated vs anti-fibrotic; 5 seeds. |
| `predict_biomarker.py` | panel: contained thinned {0,35,60,80%} + diffuse; measure **peritumoral myCAF density** (myCAF within 30µm of tumor / #tumor) + abundance; run fixed immune therapy; correlate. |
| `predict_wellmixed.py` | pure ODE mean-field (no ABM): immune kill `∝ exp(-α·M)`, no confinement; show contained==diffuse and optimum always M=0. |
| `scotia_rim_stats.py` | per-sample rim z (`run_rim_panel`, shell=30µm, n_perms=300); patient-level **Mann-Whitney U + Cliff's delta + bootstrap 95% CI + Benjamini-Hochberg** over ~10 cell types. Needs the CosMx h5ad (§7). |

---

## 6. Observation model (spec, since it's subtle)

Latent biomarker is a first-order lag of the true burden ratio:
`biomarker += (n_tumor/n0 - biomarker) * (dt / max(obs_lag_days, dt))`.
Every `obs_interval` days the controller reads `obs = biomarker * η`, η log-normal
(median 1, CV `obs_noise_cv`). Decision: force ON if `obs ≥ obs_safety_mult`; else
de-escalate to OFF only after `obs_confirm` consecutive reads `≤ adapt_off` **and**
`min_on_days` elapsed; re-escalate ON when `obs ≥ adapt_on` **and** `min_off_days`
elapsed. If `nonsecretor`, the marker is uninformative ⇒ default to continuous dosing.

---

## 7. Data acquisition (for Figs 2, 3, S6)

- **Xenium PDAC:** GEO **GSE274673** (480-gene panel). Used for the metric positive-control
  that *fails* on targeted panels (Fig 2c).
- **CosMx SMI PDAC (SCOTIA):** Mendeley Data **doi:10.17632/kx6b69n3cb.1** (Shiau et al.,
  Nat Genet 2024), 1009-gene panel, author-provided annotations. Place the `.h5ad` at
  `data/scotia/raw_meta_data_final.h5ad`. Pixel→µm scale **0.12028 µm/px**. Cell-type
  labels from `annotation_majortypes` / `annotation_subtypes`; treatment from
  `treatment_status` (Untreated / CRT; the single **CRTL** sample is excluded from the
  untreated-vs-CRT comparison).
- The synthetic-tissue analyses (most figures) need **no** downloads.

---

## 8. Verification — expected results (sanity targets)

A faithful rebuild should reproduce these to within seed noise:

| Analysis | Expected |
|---|---|
| Sobol (final burden) | `k_prolif` S_T≈0.73 and `k_kill` S_T≈0.61 dominate; barrier/resistance terms S_T≲0.12 |
| Dense phase map | keep-stroma robust (P≥0.8) at low immune-exclusion; ~20/36 cells P(keep)≥0.6 |
| 30-seed Pareto | combinations 30/30 progression-free; single garlic 28/30; wild ginseng 1/30; mugwort/curcumin 0/30 |
| Monte Carlo | natural-API regimens control in **0–2%** of draws; gemcitabine **96%** |
| Observation model | exposure ideal ≈21 → observed ≈67 vs continuous ≈120; 56-day interval ≈90; non-secretor → continuous |
| Gem schedule | exposure 128 → 98 → 71 → 18 (continuous → intermittent → observed → ideal), all control |
| Geometry (S15) | immune therapy ~12× weaker in contained (0.73×) vs diffuse (0.06×); invasion Δ +78µm vs +51µm |
| Biomarker (S17) | peritumoral density predicts immune outcome |r|≈0.92; abundance only |r|≈0.56 |
| Well-mixed (S18) | contained == diffuse (identical); optimum always M=0 |
| Patient-level rim (Fig 3) | **no** untreated-vs-CRT difference significant after BH (myCAF p≈0.53) |

If a number is far off, check: seed plumbing (single `default_rng`), the per-step update
order, and whether `max_tumor` / `CAP` differ from the table in §5.

---

## 9. Reproduce all figures + manuscript

```bash
# synthetic-only figures (no data download needed)
python pipeline/phase_map_dense.py
python pipeline/sobol_analysis.py
python pipeline/pareto_seeds.py
python pipeline/mc_uncertainty.py
python pipeline/obs_model_analysis.py
python pipeline/gem_schedule.py
python pipeline/predict_geometry.py
python pipeline/predict_sequence.py
python pipeline/predict_biomarker.py
python pipeline/predict_wellmixed.py
python pipeline/phase_map_protumor.py
# data-dependent (need the CosMx h5ad from §7)
python pipeline/scotia_rim_stats.py
python pipeline/patient_grounded.py
# assemble the Word manuscript from manuscript.md + assets/
python pipeline/build_docx.py
```
Each script also accepts `replot` to redraw from its saved CSV without re-simulating
(e.g., `python pipeline/pareto_seeds.py replot`).

Interactive dashboard: `streamlit run app.py` (or `./run.ps1` on Windows).

---

## 10. Vibe-coding prompt (paste to an AI agent)

> You are rebuilding a spatial agent-based model of PDAC and its analysis pipeline. Work
> in stages and stop for my verification after each. **Stage 1:** implement `abm.py`
> exactly per the ODD protocol I will paste (entities, state, the fixed 11-step schedule
> with start-of-step-tree synchronous division, and every submodel equation); expose
> `DEFAULT_PARAMS, SUBSTANCES, TOXICITY, TumorABM, simulate(), control_metrics()`; do NOT
> add a composite control_score. **Stage 2:** implement `synthetic.make_tissue` with
> contained/diffuse modes at identical cell counts. **Stage 3:** implement the analyses
> in §5 with the exact settings shown, each parallelized with ProcessPoolExecutor
> (Windows spawn-safe, `if __name__=="__main__"`, ~6 workers), each writing a CSV to
> data/ and PNG to assets/ and supporting a `replot` mode. All randomness must flow
> through a single `numpy.random.default_rng(seed)`. After each stage, run the smoke test
> and compare to the verification targets in §8; if a number is off, diagnose the update
> order and seed plumbing before proceeding.

Paste, in order, into that agent: this file, then `docs/manuscript/S3_ODD_protocol.md`,
then `docs/manuscript/S1_parameters.md`, then `docs/manuscript/S2_drug_products.md`.

---

*In silico, hypothesis-generating research. Outputs are testable hypotheses, not evidence
of clinical efficacy or safety. Code archived at Zenodo doi:10.5281/zenodo.21521806.*
