# Project handoff — continue here (for a new Claude Code / dev session)

This file lets a fresh session (e.g., Claude Code on another machine) pick up the project
without prior context. It is technical and public-safe. (Personal strategy/funding notes
are kept **outside** this repo.)

## What this project is
A spatial agent-based model (ABM) of the PDAC tumor–myCAF–immune ecosystem plus an
analysis pipeline that produces a manuscript. Central thesis: **condition-dependent
stromal control** — the optimal target is the stromal *state* (regime-dependent), not
depletion vs preservation. Food-medicine-homology compounds are a hypothesis-generating
case study. In silico, hypothesis-generating; not a clinical predictor.

- **Author:** Seung-Il Kim · ORCID 0009-0007-5965-9212 · Independent Researcher
- **GitHub:** https://github.com/kusi81/pdac-coexistence-control (this is the backbone;
  clone it to continue — it carries full history and the `origin` remote)
- **Zenodo archive:** doi:10.5281/zenodo.21521806
- **Manuscript (canonical):** `docs/manuscript/manuscript.md` (build Word with
  `python pipeline/build_docx.py`)

## Where everything is
| Thing | Path |
|---|---|
| Rebuild/verify guide | `docs/REBUILD_GUIDE.md` |
| ABM spec (authoritative) | `docs/manuscript/S3_ODD_protocol.md` |
| Parameter table | `docs/manuscript/S1_parameters.md` |
| Drug-product table | `docs/manuscript/S2_drug_products.md` |
| Core model | `pipeline/abm.py`, `pipeline/synthetic.py`, `pipeline/spatial_core.py`, `pipeline/analysis.py` |
| Figure scripts | `pipeline/*.py` (one per figure; each supports `replot`) |
| Figures | `assets/*.png` · Results CSVs | `data/*.csv` |
| Cover letter | `docs/manuscript/cover_letter.md` |

## Current state (as of this handoff)
The manuscript has been through a full external-review response. Done and committed:
- Figure 4 recomputed on a **6×6 grid, 5 seeds**, with probabilistic phase uncertainty
  P(keep-stroma) (`phase_map_dense.py`).
- **control_score deleted** everywhere; replaced by a progression-constrained Pareto
  evaluation (`pareto_seeds.py`, 30 seeds).
- §3.6 leads with the **30-seed** result; single-seed ranking demoted to Fig S8
  (exploratory/superseded).
- **Monte Carlo** over compound effect / exposure weight / bioavailability / synergy
  (`mc_uncertainty.py`) — natural-compound control collapses to 0–2% vs gemcitabine 96%.
- **Sobol** global SA (`sobol_analysis.py`) — reports **final tumor burden** only
  (k_prolif ≈0.73, k_kill ≈0.61 dominant); barrier params govern *which strategy*, not burden level.
- **CA19-9 observation model** (`obs_model_analysis.py`), **single-agent schedule**
  comparison (`gem_schedule.py`), **spatial predictions** (`predict_geometry.py`,
  `predict_sequence.py`, `predict_biomarker.py`), **well-mixed ablation**
  (`predict_wellmixed.py`), **pro-tumor robustness** (`phase_map_protumor.py`).
- **Patient-level CRT stats** (`scotia_rim_stats.py`): Mann–Whitney U + Cliff's delta +
  bootstrap CI + Benjamini–Hochberg → **no untreated-vs-CRT difference is significant**
  (honest; the paper rests on the pattern common to both groups).
- ODD protocol, Table S2, Zenodo metadata (`.zenodo.json`, `CITATION.cff`), and an
  **AI-assisted-technology disclosure** (Declarations) all added.
- Editorial pass complete (cover-page, CRT causal→cross-sectional, novelty softening,
  seed table §2.8, Fig captions, DOI recorded).

## Local-only files NOT in GitHub (gitignored — carry separately or re-download)
- **`data/scotia/raw_meta_data_final.h5ad`** (~287 MB, CosMx SMI; needed for Fig 3 & S6).
  Source: Mendeley Data **doi:10.17632/kx6b69n3cb.1** (SCOTIA / Shiau et al., Nat Genet 2024).
- `data/xenium/`, `data/zhou/`, `data/pdb/`, `data/gse272362/` (~1.7 GB, earlier
  exploration; mostly superseded by the synthetic-tissue analyses). Re-downloadable per
  `docs/REBUILD_GUIDE.md` §7 and `data/README.md`.
- `.venv/` (recreate: `python -m venv .venv` then `pip install -r requirements.txt` plus
  `anndata scanpy squidpy SALib`), generated `*.docx` (regenerate with `build_docx.py`).

## How to continue on the new machine
1. `git clone https://github.com/kusi81/pdac-coexistence-control.git`
2. Recreate the venv and install deps (above).
3. Drop `raw_meta_data_final.h5ad` into `data/scotia/` (from the transfer bundle or Mendeley).
4. Smoke test: `python pipeline/phase_map_dense.py replot` (redraws Fig 4 from the saved CSV).
5. Verify a full analysis against the targets in `docs/REBUILD_GUIDE.md` §8.

## Open / possible next steps
- Finalize target journal → citation style + figure numbering to house style; export
  figures as vector (PDF/SVG) for submission; trim abstract to word limit.
- Post the bioRxiv preprint (code DOI already minted; recommended order: Zenodo → bioRxiv → journal).
- Longer-term: 3D serial-section barrier pipeline; experimental validation ladder
  (organoid/CAF co-culture) — needs a wet-lab collaborator.
