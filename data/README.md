# Data — download guide

Raw and large datasets are **not committed** to the repository (see `.gitignore`).
Download them here to reproduce the analyses. Small derived result tables
(`sens_*.csv`, `scotia/scotia_*.csv`) are tracked for figure reproducibility.

## 1. CosMx SMI PDAC — SCOTIA (primary spatial dataset)
- **Source:** Shiau et al., *Nat Genet* 2024 (PMID 39227743); Mendeley Data
  **doi:10.17632/kx6b69n3cb.1**, folder `SMI/`.
- **File needed:** `raw_meta_data_final.h5ad` (300 MB; 717,493 cells × 1,009 genes,
  author cell-type/CAF-subtype annotations, treatment status).
- **Place at:** `data/scotia/raw_meta_data_final.h5ad`
- **SHA-256:** `7c65c0719d4242dcb4cc6e6d2d5798c4ad7bed64008ac66e791395251f3a4157`
- Used by: `pipeline/scotia_posctrl.py`, `pipeline/scotia_rim.py`, `pipeline/inspect_scotia.py`.

## 2. Xenium PDAC (targeted-panel comparison)
- **Source:** GEO **GSE274673** (10x Xenium, 480-gene panel). Download the per-sample
  bundles (`output-XETG...`) and extract under `data/xenium/`.
- Samples used: 31076, 39928, 28429 (naive); 38245, 35406 (CRT).
- Used by: `pipeline/data_loader.py`, `pipeline/diag_mycaf.py`,
  `pipeline/refine_annotate.py`, `pipeline/compare_rim_6samples.py`.

## 3. PDB structures (molecular viewer)
- Fetch the following PDB IDs from RCSB (https://www.rcsb.org) into `data/pdb/`:
  `1P62` (dCK+gemcitabine), `1M17` (EGFR+erlotinib), `1BG1` (STAT3),
  `5IKR` (COX-2), `5A1F` (KDM5B), `2ITY` (EGFR+gefitinib), `1VKX` (NF-κB),
  `5C3T` (PD-L1).
- Used by: `pipeline/mol_page.py`, `pipeline/render_compare_structures.py`.

## 4. Supplementary / control datasets (optional)
- **MIBI-TOF** colorectal (negative-architecture control): via `squidpy.datasets.mibitof()`.
- **Zhou Visium** PDAC (spot-level, limitation example): figshare (recovered via
  `pipeline/salvage_zip.py`). Not required for main results.

## Tracked derived tables
- `data/sens_rcost.csv`, `data/sens_tornado.csv` — sensitivity analysis (Fig. S3).
- `data/scotia/scotia_posctrl.csv` — myCAF positive control per sample (Fig. 2b).
- `data/scotia/scotia_rim.csv` — peritumoral rim composition (Fig. 3).
