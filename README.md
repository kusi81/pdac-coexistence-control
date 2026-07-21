# Coexistence-control of pancreatic cancer with food-medicine-homology compounds

A spatially grounded **agent-based framework** that models the pancreatic ductal
adenocarcinoma (PDAC) tumor–myCAF–immune ecosystem and uses **food-medicine-homology
compounds** (foods that are also traditional East-Asian medicines — garlic, ginseng,
Platycodon, mugwort, curcumin, …) to explore **tumor coexistence/control rather than
eradication**, treating the myofibroblastic-CAF (myCAF) barrier as a *controllable
resource*.

> **Status:** in-silico, hypothesis-generating research. Manuscript draft in
> [`docs/manuscript/manuscript.md`](docs/manuscript/manuscript.md).
> **Author:** Seung-Il Kim (김승일) · ORCID [0009-0007-5965-9212](https://orcid.org/0009-0007-5965-9212) · Independent Researcher.

## What this is

Three research lines are individually mature but were not connected: (i) *static*
network pharmacology of food-medicine-homology compounds, (ii) wet-lab studies of
herbal modulation of PDAC CAFs, and (iii) spatial/agent-based models of the tumor
microenvironment. This project bridges them: it grounds a dynamic spatial ABM in real
single-cell spatial transcriptomics (Xenium, CosMx), encodes compounds as mechanistic
parameter perturbations across a molecule→cell→tissue hierarchy, and optimizes
combinations, sequences, and doses against a **coexistence-control** objective with
explicit drug-resistance dynamics.

Outputs are **ranked, testable hypotheses** — low-toxicity regimens and schedules —
intended to focus, not replace, experimental validation.

## Key results (see the manuscript)

- The containment metric reproduces the myCAF-proximal positive control in **15/16**
  author-annotated CosMx tumors; the same metric fails on targeted-panel Xenium module
  scores, localizing the problem to *annotation, not metric* (Fig. 2).
- myCAFs form a tumor-adjacent barrier that **tightens under chemoradiation** while
  cytotoxic CD8⁺ T cells remain excluded (Fig. 3).
- Adaptive scheduling achieves coexistence at **~1/5 the toxicity** of continuous
  dosing (Fig. 4); food-medicine-homology regimens control the tumor at **~1/25 the
  toxicity** of gemcitabine, with anti-CAF backbones minimizing resistance (Fig. 5–6).
- The adaptive advantage is robust to ±50% parameter variation (Fig. S3).

## Repository structure

```
pipeline/        analysis & simulation code (ABM, spatial metrics, figures)
docs/manuscript/ manuscript.md + section drafts + Supplementary Table S1
docs/literature_search/  systematic PubMed survey (queries, results, novelty)
assets/          figures (PNG)
data/            small derived CSVs (large raw data: see data/README.md)
app.py, run.ps1  interactive Streamlit dashboard
```

## Install

```bash
python -m venv .venv
# Windows: .\.venv\Scripts\python.exe -m pip install -r requirements.txt
# Unix:    ./.venv/bin/pip install -r requirements.txt
```
Python 3.13; core: numpy, scipy, pandas, matplotlib. Optional (spatial data):
anndata, scanpy, squidpy, pyarrow.

## Reproduce the figures

Download data first (see [`data/README.md`](data/README.md)), then:

| Figure | Script |
|---|---|
| Fig 1 — novelty landscape | `python pipeline/fig1.py` |
| Fig 2 — metric validation | `python pipeline/fig2.py` |
| Fig 3 — rim composition (SCOTIA) | `python pipeline/scotia_rim.py` |
| Fig 4 — control strategies | `python pipeline/compare_control.py` |
| Fig 5 — natural regimens | `python pipeline/optimize_natural.py` |
| Fig 6 — dose × drug-holiday | `python pipeline/optimize_dose_band.py` |
| Fig 7 — molecular structures | `python pipeline/render_compare_structures.py` |
| Fig S3 — sensitivity analysis | `python pipeline/sens_analysis.py` |
| Positive control (15/16) | `python pipeline/scotia_posctrl.py` |
| Systematic literature survey | `python pipeline/lit_search.py` |

## Interactive dashboard

```powershell
./run.ps1   # or: .\.venv\Scripts\python.exe -m streamlit run app.py
```
Four modes: spatial analysis, perturbation simulation, control optimization, and a
molecular binding viewer. Opens at <http://localhost:8501>.

## Data

Large datasets are not committed; download and place them per
[`data/README.md`](data/README.md): SCOTIA CosMx (Mendeley doi:10.17632/kx6b69n3cb.1),
Xenium PDAC (GEO GSE274673), and PDB structures (RCSB).

## Citation

If you use this code or framework, please cite the associated manuscript (in
preparation). BibTeX will be added on preprint/publication.

## License

Code: MIT (see [`LICENSE`](LICENSE)). Manuscript text and figures: research content by
the author — please cite the paper if reused.

---

### 한글 요약
췌장암(PDAC)에서 **myCAF 장벽을 통제 자원으로 활용**해 종양을 박멸이 아닌 **공존·통제**로
다루는 공간 agent-based 프레임워크입니다. 약식동원(식품이자 약재인) 화합물을 분자→세포→조직
다중스케일 섭동으로 인코딩하고, 적응요법·내성 동역학과 함께 조합·순차·용량을 최적화합니다.
결과는 **검증 대상 가설**(저독성 배합·스케줄)이며, 실험 검증을 대체하지 않습니다.
논문 초고: `docs/manuscript/manuscript.md`.
