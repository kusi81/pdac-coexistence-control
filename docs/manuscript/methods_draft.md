# Methods (draft v1)

> **작성 메모:** 실제 구현(pipeline/*.py)에 정확히 대응. 파라미터·문헌접지는 abm.py
> DEFAULT_PARAMS 주석 기반. 정직: 현상학적 ABM(데이터 fit 아님·문헌접지), 2D, 독성 임의단위.
> 전체 파라미터·인용은 Supp Table S1. [n]=Introduction 참조번호.

---

## 2. Methods

### 2.1 Systematic novelty survey
We queried PubMed via the NCBI E-utilities API with ten Boolean queries spanning
the intersecting dimensions of the approach (adaptive therapy, agent-based/
computational modeling, network pharmacology, food-medicine-homology compounds,
myCAF/CAF biology, PDAC, spatial transcriptomics, resistance/coexistence), plus
two "integrated" queries requiring three axes simultaneously. For each query we
retrieved hit counts and up to 60 records (relevance-sorted), fetching abstracts
for the union (132 unique articles). Queries, counts, and retrieved records are in
Supplementary Data. This is a targeted novelty search, not a PRISMA systematic
review, and is restricted to PubMed.

### 2.2 Spatial transcriptomic data and cell annotation
We used two single-cell spatial PDAC datasets. (i) **Xenium** (10x Genomics;
GSE274673): five resected human PDAC samples (three treatment-naive, two
chemoradiation-treated [CRT]) with a 480-gene panel; cell coordinates (µm) from
`cells.parquet` and counts from `cell_feature_matrix.h5`. Cells were typed by a
two-stage module-score scheme (coarse type by argmax over marker modules, then CAF
subtype among CAF-assigned cells) using literature marker panels. (ii) **CosMx SMI**
(NanoString; Shiau et al. [20]; Mendeley Data kx6b69n3cb.1): 16 PDAC tumors
(9 untreated, 6 CRT, 1 CRTL) with a 1,009-gene panel and **author-provided**
major-type and CAF-subtype annotations (717,493 cells); global pixel coordinates
were converted to microns (0.12028 µm/px). MIBI-TOF (colorectal, squidpy) served
as a negative-architecture control. The two PDAC datasets play complementary
roles: the author-annotated CosMx cohort provides the annotation ground truth for
metric validation (§3.2), while the Xenium cohort illustrates the limits of
marker-based typing.

### 2.3 Barrier and containment metrics
All spatial metrics operate on (coordinates, cell-type labels) and use
permutation nulls that preserve tissue geometry. (i) **Barrier score**: the
fraction of straight-line tumor→immune paths that pass within a corridor
(default 30 µm) of a barrier-cell (myCAF), compared against a matched-null that
resamples barrier positions at fixed density; reported as a z-score, so the metric
tests *interposition geometry* rather than abundance. (ii) **Rim enrichment**: for
a peritumoral shell (30 µm), the enrichment of each cell type relative to a
label-permutation null, as a z-score (>2 tumor-adjacent, <−2 excluded).
(iii) **Proximity test**: whether cell type A is closer than cell type B to the
nearest tumor cell (median nearest-tumor distance; label-permutation null over the
A/B pool), used as the myCAF-proximal/iCAF-distal positive control [3,20].
Permutations: 200–500. For CosMx, distances are reported in microns after pixel
conversion.

### 2.4 Agent-based control model
The core model is a phenomenological, off-lattice agent-based model (ABM) of the
PDAC tumor–myCAF–immune ecosystem on a 2D field (default 1,500 µm), advanced in
0.5-day steps. Agents are tumor cells (sensitive or resistant), myCAFs, and CD8⁺ T
cells (with iCAF/macrophage as context). Update rules per step:

- **Tumor proliferation.** Each tumor cell divides with probability `k_prolif`
  (baseline 0.11/day; doubling ~6 days), gated by local contact inhibition
  (≤`tumor_density_cap` neighbors within the kill radius); baseline apoptosis
  `k_tumor_apoptosis`.
- **myCAF barrier.** Peri-tumoral stroma activates to myCAF with probability
  `k_caf_activate` within `caf_ring_um` (150 µm) of tumor, up to
  `caf_cap_per_tumor`, and reverts (turnover) at `k_caf_death`; anti-fibrotic
  perturbation lowers activation so turnover reduces barrier mass. An optional
  `caf_protumor` term boosts local proliferation (0 for PDAC—pure barrier; large
  for the HCC context, reflecting cirrhosis-as-pro-tumor soil).
- **CD8 migration/killing.** CD8 cells move toward the nearest tumor at
  `cd8_speed_um`/day, attenuated by myCAF density in the corridor via a gate
  parameterized by `cd8_barrier_alpha` (the key immune-exclusion lever); they kill
  tumor within `kill_radius_um` at `k_kill`/day. Recruitment `cd8_recruit`/day at
  the tissue margin; turnover `k_cd8_death`.
- **Resistance dynamics.** A small pre-existing resistant fraction
  (`init_resistant_frac` 0.01) plus phenotypic switching (`mutation_rate`) model
  the predominantly non-genetic gemcitabine resistance of PDAC. Resistant cells
  (a) are unaffected by drug, (b) evade CD8 killing with probability
  `resistant_immune_evasion` (0.45; antigen/MHC-I loss, JEM 2017), and (c) carry a
  fitness cost `resistance_cost` (0.24; NSCLC-derived ~76% of sensitive growth,
  Science Advances), which enables competitive suppression of resistance during
  drug holidays.

Parameters are literature-grounded rather than fit to data; full values and
citations are in Supplementary Table S1. Two organ contexts (PDAC, HCC) share the
stellate/TGF-β axis and differ only in `caf_protumor`.

### 2.5 Encoding food-medicine-homology compounds
Seventeen agents—food-medicine-homology compounds (e.g., garlic, ginseng,
Platycodon, mugwort, curcumin, ginsenoside-Rg3) alongside conventional references
(gemcitabine, erlotinib) and a hypothesis-stage repurposing candidate
(entecavir)—are encoded as multiplicative perturbations of a defined set of
perturbable parameters (proliferation, CAF activation, immune killing/recruitment,
immune motility, apoptosis), scaled by dose. Each compound carries an evidence
badge and a toxicity weight used to accumulate the treatment burden. Regimens are
composed from compound–dose lists (`compose_regimen`), optionally with a
literature-motivated synergy term, and can be arranged into timed sequential
cycles (`build_cycle_schedule`: myCAF-weakening → anti-proliferative → stroma-reset
phases) via an automatic anti-fibrotic/immune vs cytotoxic/targeted classification.

### 2.6 Multiscale molecular grounding
Each compound is linked to its principal molecular target(s) with an explicit
evidence tier—experimental co-crystal, computational docking, or pathway-level
mechanistic inference (no solved structure). Structures are rendered from
verified PDB entries (e.g., gemcitabine–dCK 1P62, erlotinib–EGFR 1M17
[experimental]; curcumin–STAT3 1BG1 [docking]); the viewer distinguishes when the
displayed ligand is the compound itself (experimental) versus a reference ligand
or predicted residues (docking/mechanistic). Targets map to the ABM parameter
class they perturb, providing the molecule→cell→tissue link.

### 2.7 Control objective and metrics
The objective is durable control at minimal burden, not eradication. From each
simulated trajectory we compute: time-to-progression (first time tumor reaches
1.5× baseline; if never, the horizon is used and the run is *censored*), peak and
final burden, final resistant fraction, cumulative toxicity, burden AUC, and a
composite `control_score` = time-to-progression / (cumulative toxicity + 1).
Because `control_score` is artifactually inflated for zero-toxicity (untreated)
runs, strategies are ranked by progression control first (censored vs progressed)
and then by toxicity. Adaptive scheduling toggles treatment on above an upper band
(`adapt_on`) and off below a lower band (`adapt_off`) around the tumor-burden
reference.

### 2.8 In-silico experimental design
Simulations were run on synthetic tissues generated with matched cell counts in
two architectures—`contained` (myCAF rings around tumor islets, CD8 excluded) and
`diffuse` (identical counts, no architecture)—so that abundance is held constant
and only geometry differs. Experiments comprised: (E1) untreated vs continuous
maximum-dose vs adaptive scheduling; (E2) single agents vs combinations vs timed
sequential cycles; (E3) a seed-averaged dose × drug-holiday grid search; and a
sensitivity analysis combining a `resistance_cost` sweep with one-at-a-time ±50%
perturbation of nine key parameters. Stochastic results are averaged over three
seeds (42, 7, 123).

### 2.9 Implementation and availability
The framework is implemented in Python (3.13) using scanpy/anndata, squidpy,
pyarrow, NumPy/SciPy, and Matplotlib, with an interactive Streamlit dashboard
exposing the spatial, simulation, optimization, and molecular-viewer modes.
Spatial data: Xenium GSE274673; CosMx from Mendeley Data
(doi:10.17632/kx6b69n3cb.1) [20]. All analysis code is available at
[repository URL].

---

## 편집 체크리스트
- [ ] Supplementary Table S1 작성(전 파라미터·기본값·인용·perturbable 여부)
- [ ] 독성 "arbitrary units" 정의 명확화(물질 toxicity 가중 × dose-day 누적)
- [ ] barrier_score 수식(경로 샘플링·matched-null) 부록에 정식 기술
- [ ] repository URL·라이선스 확정
- [ ] CRTL 정의 각주(§2.2), MIBI-TOF/Zhou Visium 보조데이터 위치
- [ ] 인용번호 [3][20] 등 Introduction과 일치 확인
