# Condition-dependent stromal control of pancreatic cancer in a spatially grounded agent-based model

*Subtitle: food-medicine-homology compounds as a hypothesis-generating application*

*In-silico, hypothesis-generating study.*

**Author:** Seung-Il Kim (김승일 / 金昇日)¹
**Corresponding author:** Seung-Il Kim · kusi81kim@gmail.com · **ORCID:** [0009-0007-5965-9212](https://orcid.org/0009-0007-5965-9212)
¹ Independent Researcher (unaffiliated).

---

## Abstract

Pancreatic ductal adenocarcinoma (PDAC) combines dismal survival with a dense,
myofibroblastic cancer-associated fibroblast (myCAF) stroma that both supports and
physically constrains the tumor. Because attempts to ablate this stroma have
paradoxically accelerated disease, we ask whether the myCAF barrier can instead be
treated as a *controllable resource* for restraining, rather than eradicating, the
tumor. We present an in-silico framework that integrates three previously
disconnected research lines—static network pharmacology of "food-medicine
homology" compounds, wet-laboratory studies of herbal modulation of PDAC CAFs, and
agent-based modeling of the tumor microenvironment. A targeted PubMed survey
found that, while each line is individually well populated, no directly matching study
integrates them. Our framework builds a spatial agent-based model of the
PDAC tumor–myCAF–immune ecosystem, grounded in single-cell spatial transcriptomics
(Xenium and CosMx), and encodes food-medicine-homology compounds (e.g., garlic,
ginseng, Platycodon, mugwort) as mechanistic parameter perturbations across a
molecular-to-cellular-to-tissue hierarchy. Critically, the model represents the myCAF
barrier physically—local stroma confines tumor expansion but also excludes cytotoxic
T cells and impairs drug delivery—so preserving stroma is not universally beneficial.
Under a coexistence-control objective with explicit resistance dynamics, a phase
analysis shows that treating myCAF as a resource to modulate (rather than deplete) is
favored only where physical confinement outweighs immunosuppression; where
immunosuppression dominates, stromal reduction is preferred, making the optimal
stromal state—not depletion or preservation per se—the central control target and the
principal contribution of this work. A global (Sobol) sensitivity analysis identifies
tumor proliferation and the immune-exclusion barrier as the dominant controls. As a
hypothesis-generating application of the framework, we encode food-medicine-homology
compounds (a case study, not an efficacy claim) and find that adaptively scheduled,
low modeled-exposure regimens achieved simulated control under the specified
assumptions; because inputs such as compound exposure weights are assigned rather than
measured, these are strategy comparisons rather than evidence of clinical efficacy or
safety. We stress that all outputs are testable hypotheses—candidate stromal targets,
combinations, and schedules—intended to focus, not replace, experimental validation.
The framework offers a reproducible route from real spatial data to condition-dependent
stromal-control principles for PDAC and a prioritization engine for the experiments
that must follow.

**Keywords:** pancreatic ductal adenocarcinoma; cancer-associated fibroblasts;
myCAF; adaptive therapy; tumor coexistence; agent-based model; spatial
transcriptomics; food-medicine homology; network pharmacology; drug resistance

## Author Summary

Pancreatic cancer is one of the deadliest cancers, and patients often lose healthy
time very quickly. Most treatments aim to destroy the tumor outright, but this
frequently fails and carries heavy toxicity. The tumor is wrapped in a dense
"scar-like" tissue built by cells called myofibroblastic cancer-associated
fibroblasts (myCAFs). Surprisingly, experiments that tried to remove this tissue
made the cancer worse—suggesting the barrier also holds the tumor back. We asked a
different question: instead of destroying the tumor, can we *manage* it—slowing its
growth and using the myCAF barrier to keep it contained, so patients feel better
and live longer? To explore this safely and cheaply, we built a computer simulation
of pancreatic tumor tissue based on real spatial measurements of where each cell
sits. The main finding is a principle: whether it helps to keep or to thin the
barrier depends on the individual tumor—on whether the barrier is mostly walling the
tumor in, or mostly keeping immune cells and drugs out. As a worked example of how the
model can be used, we added compounds that are eaten as food and used as traditional
medicine in Northeast Asia and let the model rank promising, low-exposure strategies.
These rankings are illustrative hypotheses generated under assumed inputs, not claims
that these foods treat cancer, and are meant to guide future laboratory experiments,
not to replace them.

---

## 1. Introduction

Pancreatic ductal adenocarcinoma (PDAC) remains among the most lethal human
malignancies, with a five-year survival of roughly 13% and a median survival
measured in months for advanced disease [1]. A defining feature of PDAC is its
extensive desmoplastic stroma, in which cancer-associated fibroblasts (CAFs) and a
dense extracellular matrix can constitute the majority of the tumor mass [2]. For
patients, the clinical reality is not only high mortality but a steep, often rapid
decline in remaining quality-adjusted life: the pace at which PDAC erodes patient
time and comfort is itself a therapeutic target. Interventions that meaningfully
*slow* tumor expansion—even without eradication—could translate into reduced
symptom burden and extended functional survival.

CAFs in PDAC are not a monolithic population. Single-cell and spatial studies have
resolved at least three functional states—myofibroblastic CAFs (myCAF),
inflammatory CAFs (iCAF), and antigen-presenting CAFs (apCAF)—that occupy distinct
spatial niches and are, importantly, plastic and interconvertible [3,4,5]. myCAFs,
enriched immediately adjacent to tumor epithelium, deposit the fibrotic matrix that
both supports and physically constrains the tumor. This duality has a critical
therapeutic corollary: attempts to *ablate* the stroma in preclinical PDAC
paradoxically accelerated disease and worsened survival [6,7], reframing the
field's question from "should CAFs be deleted?" to "should they be reeducated?"
[8]. We take this further and ask whether the myCAF-derived barrier can act as a
*controllable resource*—a containment structure whose value is conditional: its
physical confinement of the tumor is a benefit, but the same density excludes
cytotoxic T cells and impairs drug delivery, which are costs. Whether preserving or
reducing stroma aids control therefore depends on which force dominates locally—a
trade-off we make explicit in the model and map, rather than assume.

This perspective aligns naturally with a shift in therapeutic philosophy from
maximum-tolerated-dose (MTD) eradication toward evolutionary *control*. Adaptive
therapy, which maintains a population of treatment-sensitive cells to competitively
suppress resistant subclones through intermittent or modulated dosing, has shown
that stabilizing tumor burden can outperform attempts to minimize it, provided
resistant cells carry a fitness cost [9,10,11,12]. Related control-theoretic
analyses now frame tumor management explicitly in terms of *controllability* under
co-evolving resistance [13]. In PDAC—where MTD regimens impose severe toxicity for
limited benefit—an approach oriented toward durable coexistence, low toxicity, and
preserved quality of life is especially compelling, yet the spatial dimension of
CAF-mediated containment has not been integrated into such control frameworks.

A complementary opportunity lies in the compounds available for long-term,
low-toxicity control. "Food-medicine homology" (medicine-food homology)
substances—plant and animal materials used in Northeast Asian traditions as both food
and remedy (e.g., garlic, ginseng, Platycodon, mugwort, hawthorn)—are reputed to
offer multi-target activity with a favorable safety and accessibility profile, making
them candidate backbones for chronic control regimens rather than acute cytotoxic
bursts [14,15]; we treat that reputation as a hypothesis to be tested, not an
established fact. However, the computational study of these compounds in oncology has
been dominated by *static* network pharmacology: compound–target–pathway mapping
and molecular docking that predict mechanism but do not simulate dynamics, space,
or scheduling [14,15,16]. Our own systematic PubMed survey (Methods) confirmed this
imbalance quantitatively—food-medicine-homology network-pharmacology cancer studies
are abundant, whereas their integration with dynamic, spatial, or adaptive-control
modeling is essentially absent.

Indeed, the relevant prior art separates into three camps that, to our knowledge,
have not been bridged. First, static network-pharmacology studies of
food-medicine-homology compounds against cancer establish molecular targets but
remain time- and space-agnostic [16]. Second, wet-laboratory studies show that
specific herbal compounds and formulas modulate PDAC CAFs—suppressing CAF
activation [17] and reversing CAF-induced gemcitabine resistance [18]—providing
biological grounding but no predictive or optimization framework. Third,
agent-based and multiscale models of the tumor microenvironment, increasingly
informed by spatial transcriptomics, capture tumor–stroma–immune dynamics but do
not incorporate natural-compound perturbations or a coexistence objective [19,20].
No existing work, to our knowledge, connects the compound knowledge of the first
camp to the dynamic spatial models of the third in service of the CAF-containment
biology of the second.

Here we present an in-silico framework that occupies precisely this gap. We build a
spatial agent-based model of the PDAC tumor–myCAF–immune ecosystem, grounded in
real single-cell spatial transcriptomics (Xenium and CosMx) [20], in which
food-medicine-homology compounds are encoded as mechanistic parameter perturbations
acting across a molecular-to-cellular-to-tissue hierarchy. The model represents the
myCAF barrier physically—local stromal density confines tumor expansion, excludes
cytotoxic T cells, and impairs drug penetration—so that stromal modulation carries an
explicit benefit–cost trade-off. Because it represents a single contained focus, the
framework's scope is localized and locally advanced disease, where local control and
relief of mass-effect complications are the goals; it does not model disseminated or
peritoneal metastatic disease (Discussion §4.4). Rather than optimizing for eradication,
we define a *coexistence-control* objective with explicit resistance dynamics and a
resistance
fitness cost, and we first map, over the confinement-versus-immunosuppression
trade-off, the conditions under which preserving/modulating the stroma outperforms
depleting it. This condition-dependent stromal-control principle—that the optimal
stromal state is regime-dependent rather than uniformly "preserve" or "deplete"—is the
paper's central result. As a downstream application that demonstrates how the framework
turns into concrete, testable predictions, we then use it to prioritize compound
combinations, sequences, and doses under adaptive scheduling; the food-medicine-homology
compounds serve as a case study for this prioritization, not as a claim of therapeutic
efficacy or safety. We emphasize throughout that the framework is a
hypothesis-*generating* and hypothesis-*prioritizing* engine: its outputs are testable
predictions—candidate stromal targets, compound combinations, and schedules—intended
to focus, not replace, subsequent experimental validation.

---

## 2. Methods

### 2.1 Systematic novelty survey
We queried PubMed via the NCBI E-utilities API with ten Boolean queries spanning
the intersecting dimensions of the approach (adaptive therapy, agent-based/
computational modeling, network pharmacology, food-medicine-homology compounds,
myCAF/CAF biology, PDAC, spatial transcriptomics, resistance/coexistence), plus two
"integrated" queries requiring three axes simultaneously. For each query we
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
(NanoString; Shiau et al. [20]; Mendeley Data doi:10.17632/kx6b69n3cb.1): 16 PDAC
tumors (9 untreated, 6 CRT, and 1 specimen labeled CRTL—chemoradiation plus an additional
line of therapy in the source metadata—which we excluded from the untreated-versus-CRT
rim comparison for group clarity) with a 1,009-gene panel and **author-provided**
major-type and CAF-subtype annotations (717,493 cells); global pixel coordinates
were converted to microns (0.12028 µm/px). MIBI-TOF (colorectal, squidpy) served as
a negative-architecture control. The two PDAC datasets play complementary roles: the
author-annotated CosMx cohort provides the author-provided reference annotation used
as a positive control for the metric (§3.2), while the Xenium cohort illustrates the
limits of marker-based typing.

### 2.3 Barrier and containment metrics
All spatial metrics operate on (coordinates, cell-type labels) and use permutation
nulls that preserve tissue geometry. (i) **Barrier score** (a *myCAF spatial
interposition index*, not a measured impediment to cell movement): the fraction of
straight-line tumor→immune paths that pass within a corridor (default 30 µm) of a
barrier-cell (myCAF), compared against a matched-null that resamples barrier
positions at fixed density; reported as a z-score, so the metric tests
*interposition geometry* rather than abundance. The straight-line construction is a
2D approximation; true CD8 trafficking additionally depends on vascular entry points,
chemokine gradients, and ECM fiber orientation (Discussion §4.5). (ii) **Rim enrichment**: for a
peritumoral shell (30 µm), the enrichment of each cell type relative to a
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
cells (with iCAF/macrophage as context). A complete ODD-protocol specification (entities
and state variables, the fixed per-step process schedule and update semantics, design
concepts, initialization, and every submodel equation) is provided in the Supplement
(`S3_ODD_protocol.md`); the essential update rules per step are:

- **Tumor proliferation.** Each tumor cell divides with probability `k_prolif`
  (baseline 0.11/day; doubling ~6 days), gated by local contact inhibition
  (≤`tumor_density_cap` neighbors within the kill radius); baseline apoptosis
  `k_tumor_apoptosis`.
- **myCAF barrier (activation/turnover).** Peri-tumoral stroma activates to myCAF with
  probability `k_caf_activate` within `caf_ring_um` (150 µm) of tumor, up to
  `caf_cap_per_tumor`, and reverts (turnover) at `k_caf_death`; anti-fibrotic
  perturbation lowers activation so turnover reduces barrier mass. An optional
  `caf_protumor` term boosts local proliferation (baseline 0; enabled in the pro-tumor
  robustness analysis of §3.4/Fig. S11, and available for co-opted-stroma settings such as
  a cirrhotic liver).
- **myCAF physical containment.** Local myCAF density ρ (within the kill radius) exerts
  three effects on the tumor. (i) *Confinement*: a dividing cell's daughter is blocked
  from being placed into a stroma-dense location with probability
  `caf_confine`·min(ρ/`caf_confine_ref`, 1), so the tumor cannot expand through the
  stromal wall. (ii) *Mechanical pressure*: the local carrying capacity is scaled by
  exp(−`caf_pressure`·ρ/`caf_confine_ref`), so the tumor cannot pack where stroma is
  dense. (iii) *Impaired drug delivery*: the drug's anti-proliferative effect on
  sensitive cells is attenuated by exp(−`caf_drug_block`·ρ/`caf_confine_ref`). Setting
  these to zero recovers a model in which myCAF affects only immunity.
- **CAF pro-tumor support (paracrine growth and drug tolerance).** To avoid biasing the
  model toward containment, CAFs can also *support* the tumor via two axes used in a
  robustness analysis (§3.4, Fig. S11): a paracrine proliferation/carrying-capacity boost
  (`caf_protumor`, the "pro-tumor soil" term) and a signaling-based drug-tolerance term
  (`caf_survival`) by which local myCAF further attenuates the drug's effect on sensitive
  cells—multiplicatively with the physical `caf_drug_block` term but representing distinct
  biology (IL-6/JAK–STAT-type therapy-induced survival signaling rather than physical
  exclusion). Both default to zero in the baseline and are switched on to test robustness.
- **CD8 migration/killing.** CD8 cells move toward the nearest tumor at `cd8_speed_um`
  /day, attenuated by myCAF density in the corridor via a gate parameterized by
  `cd8_barrier_alpha` (the immune-exclusion lever); they kill tumor within
  `kill_radius_um` at `k_kill`/day. Recruitment `cd8_recruit`/day at the tissue
  margin; turnover `k_cd8_death`.
- **Resistance dynamics.** A small pre-existing resistant fraction (`init_resistant_frac`
  0.01) plus phenotypic switching (`mutation_rate`) model the predominantly
  non-genetic gemcitabine resistance of PDAC. Resistant cells (a) are unaffected by
  drug, (b) evade CD8 killing with probability `resistant_immune_evasion` (0.45;
  antigen/MHC-I loss, JEM 2017), and (c) carry a fitness cost `resistance_cost` (0.24;
  NSCLC-derived ~76% of sensitive growth, Science Advances), which enables competitive
  suppression of resistance during drug holidays.

The myCAF barrier thus confers one benefit—physical confinement (effects i–ii above)—
against costs: immune exclusion (via `cd8_barrier_alpha`), impaired drug delivery
(effect iii), and, when the pro-tumor axes are enabled, paracrine growth support and
signaling-based drug tolerance. Their balance determines whether preserving or reducing
stroma aids control, which we map in §3.4 (and re-map with the pro-tumor axes on in
Fig. S11). Parameters are literature-grounded rather than fit to
data; full values and citations are in Supplementary Table S1. (The stellate/TGF-β axis
generalizes to other desmoplastic organs; results here are for PDAC only.)

### 2.5 Encoding food-medicine-homology compounds
Seventeen agents—food-medicine-homology compounds (e.g., garlic, ginseng, Platycodon,
mugwort, curcumin, ginsenoside-Rg3) alongside conventional references (gemcitabine,
erlotinib) and a hypothesis-stage repurposing candidate (entecavir)—are encoded as
multiplicative perturbations of a defined set of perturbable parameters
(proliferation, CAF activation, immune killing/recruitment, immune motility,
apoptosis), scaled by dose. Each compound carries an evidence badge and a toxicity
weight used to accumulate the treatment burden. Regimens are composed from
compound–dose lists, optionally with a literature-motivated synergy term, and can be
arranged into timed sequential cycles (myCAF-weakening → anti-proliferative →
stroma-reset phases) via an automatic anti-fibrotic/immune vs cytotoxic/targeted
classification.

### 2.6 Multiscale molecular grounding
Each compound is linked to its principal molecular target(s) with an explicit
evidence tier—experimental co-crystal, computational docking, or pathway-level
mechanistic inference (no solved structure). Structures are rendered from verified
PDB entries (e.g., gemcitabine–dCK 1P62, erlotinib–EGFR 1M17 [experimental];
curcumin–STAT3 1BG1 [docking]); the viewer distinguishes when the displayed ligand is
the compound itself (experimental) versus a reference ligand or predicted residues
(docking/mechanistic). Targets map to the ABM parameter class they perturb, providing
the molecule→cell→tissue link.

### 2.7 Control objective and metrics
The objective is durable control at minimal burden, not eradication. From each
simulated trajectory we compute a panel of outcomes that we treat as a
multi-objective (Pareto) problem rather than collapsing to one number:
time-to-progression (first time tumor reaches 1.5× baseline; if never, the horizon is
used and the run is *censored*), peak and final burden, final resistant fraction,
cumulative modeled exposure, and burden AUC. We deliberately do **not** collapse these
into a single composite score: a score that rewards low exposure (e.g., time-to-progression
divided by exposure) would rank the zero-exposure untreated arm highest despite
uncontrolled growth. Instead we use a **progression-constrained multi-objective
evaluation**—rank first by progression control (censored vs progressed), then, among
progression-free strategies, by cumulative exposure, i.e., by position on the
progression-versus-exposure Pareto frontier (§3.6, Fig. S10)—with resistant fraction and
burden-AUC as secondary read-outs. Adaptive scheduling
toggles treatment on above an upper band (`adapt_on`) and off below a lower band
(`adapt_off`) around the tumor-burden reference. Because this feedback uses the true
simulated burden with zero measurement delay or noise, it is an idealized controller;
we therefore also implement an explicit CA19-9-based observation model (interval
sampling, measurement noise, biomarker-to-burden lag, non-secretors, minimum treatment
durations, and a confirm-before-stop safety rule) and test whether the adaptive advantage
survives it (§3.5, Fig. S12; Discussion §4.5).

### 2.8 In-silico experimental design
Simulations were run on synthetic tissues generated with matched cell counts in two
architectures—`contained` (myCAF rings around tumor islets, CD8 excluded) and
`diffuse` (identical counts, no architecture)—so that abundance is held constant and
only geometry differs. Experiments comprised: (E1) untreated vs continuous
maximum-dose vs adaptive scheduling—a within-model comparison of exposure schedules,
where "continuous maximum-dose" is a modeling reference arm and not a representation
of clinical standard-of-care; (E2) single agents vs combinations vs timed
sequential cycles; (E3) a seed-averaged dose × drug-holiday grid search; and a
sensitivity analysis combining a `resistance_cost` sweep with one-at-a-time ±50%
perturbation of nine key parameters, complemented by a variance-based global (Sobol)
analysis over eight parameters (Saltelli sampling, first- and total-order indices;
Fig. S7). Exploratory results use three seeds (42, 7, 123); the confirmatory analyses use
larger replication as summarized below.

**Replication per experiment.** Seed/sample counts vary by analysis and are reported with
each result:

| Analysis | Figure | Replication |
|---|---|---|
| Exploratory sims, phase map | Figs 4, 5, S3, S8, S9 | 3 seeds |
| Fair 2×2 (agent × schedule) | Fig. S5 | 5 seeds |
| GV1001-style barrier gating | Fig. S4 | (regime sweep) |
| Global Sobol sensitivity | Fig. S7 | Saltelli N=32, 1 seed/run |
| Multi-seed Pareto ranking | Fig. S10 | 30 seeds |
| Pro-tumor phase-map robustness | Fig. S11 | 3 seeds/cell, 2 conditions |
| CA19-9 observation model | Fig. S12 | 20 seeds |
| Single-agent schedule | Fig. S13 | 20 seeds |
| Monte Carlo epistemic uncertainty | Fig. S14 | 100 draws |
| Non-trivial spatial predictions | Figs S15–S17 | 5 seeds |

### 2.9 Implementation and availability
The framework is implemented in Python (3.13) using scanpy/anndata, squidpy, pyarrow,
NumPy/SciPy, and Matplotlib, with an interactive Streamlit dashboard exposing the
spatial, simulation, optimization, and molecular-viewer modes. Spatial data: Xenium
GSE274673; CosMx from Mendeley Data (doi:10.17632/kx6b69n3cb.1) [20]. All analysis
code is available at https://github.com/kusi81/pdac-coexistence-control; a version-tagged
release archived with a DOI (e.g., Zenodo) will accompany publication so that every figure
is reproducible from the exact commit, including per-figure generation scripts and the
random seeds listed in §2.8.

---

## 3. Results

### 3.1 The components of our approach are individually well-established, but no directly matching integration was found in a targeted survey
To position the framework, we ran a targeted PubMed survey (twelve Boolean queries—ten
spanning the intersecting dimensions of our approach plus two "integrated" queries;
Methods §2.1) and quantified
how much prior literature occupies each axis versus their integration (Fig. 1a). The
individual axes are well populated: coexistence/control with natural compounds (56
records), network-pharmacology of food-medicine-homology compounds against cancer
(33), and spatial/agent-based modeling of PDAC (31) each returned substantial
literature. In contrast, queries requiring the *simultaneous* combination of dynamic/
computational modeling, natural compounds, and CAF biology returned only 2–3 records,
and on inspection all were off-target—reviews or wet-laboratory studies in other tumor
types—none of which couples a spatial/dynamic model to natural-compound perturbation
of PDAC CAF biology.

The relevant prior art therefore resolves into three camps that our survey shows are
individually mature but mutually disconnected (Fig. 1b): (i) static network
pharmacology of food-medicine-homology compounds, which maps compound–target
relationships but is time- and space-agnostic [14,16]; (ii) wet-laboratory studies of
herbal modulation of PDAC CAFs, which provide biology but no predictive or
optimization framework [17,18]; and (iii) agent-based and spatial models of the PDAC
microenvironment, which capture tumor–stroma–immune dynamics but incorporate neither
natural compounds nor a coexistence-control objective [19,20]. The triple
intersection—a dynamic, spatially grounded, adaptive-control framework that uses
food-medicine-homology compounds to modulate the myCAF barrier—returned no directly
matching study in our targeted survey, indicating the gap this work addresses. We note the limits of this analysis: it
queries PubMed only and is a targeted novelty search rather than a PRISMA systematic
review; nonetheless, the near-absence of integrated hits amid abundant single-axis
literature supports that the contribution lies in the integration.

### 3.2 The containment metric is supported by an author-annotated positive control, and localizes annotation—not metric—as the limiting factor
On synthetic tissues, the barrier score cleanly separated an architecture in which
stroma is interposed between tumor and immune cells (contained; z ≈ +22.3) from one in
which the same cell counts are randomly intermixed (diffuse; z ≈ +0.9), confirming
that the matched-null construction responds to spatial *geometry* rather than
abundance (Fig. 2a). We next tested a known biological positive control—myCAFs lie
closer to malignant cells than iCAFs [3,20]—on the author-annotated CosMx dataset
(717,493 cells across 16 tumors). myCAF was significantly closer to malignant cells
than iCAF in 15 of 16 tumors (the single exception showed a 0.2 µm difference,
p = 0.279, a tie), and the gap was larger in CRT-treated than untreated specimens
(−2.6 µm untreated vs −18.8 µm CRT; cross-sectional, not paired) (Fig. 2b). Given faithful
annotations, our metric recovers the expected architecture.

Critically, the *same* metric failed this positive control on our targeted-panel
Xenium data typed by marker-based module scores: myCAF scored as marginally *farther*
from tumor than iCAF (positive control failed in 4–5 of 5 tumors). We traced this to
the annotation: the 480-gene panel contained only 4 of 8 canonical myCAF markers, and
the two present pan-myofibroblast markers (ACTA2, TAGLN) are shared with perivascular/
smooth-muscle cells; a refined annotation adding gold-standard myCAF markers (POSTN,
LRRC15, CTHRC1) and a dedicated perivascular class still failed (Fig. 2c). Together
these establish that our metrics are sound—reproducing the literature positive control
on full-panel author-annotated data—and that robust CAF-subtype spatial analysis
requires full-transcriptome or author-level annotations rather than targeted-panel
module scoring, consistent with our hypothesis-generating framing.

### 3.3 A tumor-adjacent myCAF barrier with excluded cytotoxic T cells is present in both untreated and CRT specimens
We characterized the peritumoral rim (30 µm shell) composition in the CosMx cohort with
the **patient as the unit of analysis** (9 untreated, 6 CRT), computing each patient's
per-cell-type rim-enrichment z, then comparing groups with a nonparametric test
(Mann–Whitney U), a Cliff's-delta effect size, patient-bootstrap 95% confidence intervals,
and Benjamini–Hochberg correction across the ten cell types (Fig. 3). This is a
cross-sectional comparison of different patients, not a paired longitudinal measurement.

The **robust, group-independent** finding is architectural: myCAFs are strongly
rim-enriched (untreated mean z = +9.1 [95% CI −0, 19]; CRT +15.7 [9, 24]) while cytotoxic
CD8⁺ T cells (−5.6 [−8, −3] vs −4.0 [−6, −2]) and other lymphocytes are excluded—in *both*
groups. In contrast, **no untreated-versus-CRT difference was statistically significant**:
although several point estimates trend as expected (denser myCAF rim, more excluded iCAF,
Cliff's delta −0.63, and rim-ward pericyte shift under CRT), the largest raw effects reach
only p ≈ 0.05 and **none survives multiple-testing correction** (all BH-adjusted p > 0.19;
myCAF p = 0.53). With n = 9 versus 6 and wide patient-level scatter, the cohort is
underpowered to establish a treatment effect on rim composition. We therefore report the
untreated-versus-CRT contrasts as **suggestive, non-significant associations**, not
treatment effects, and rest the spatial claim on the pattern common to both groups: a
persistent myCAF-rim barrier coincident with CD8⁺ T-cell exclusion.

This pattern refines the intuitive expectation that therapy relieves immune exclusion: in
neither group are cytotoxic T cells at the tumor rim, and the myCAF barrier is present in
both. Rather than contradicting the containment hypothesis, this supports it—the myCAF
barrier and immune exclusion co-occur regardless of treatment status, maintaining spatial
separation between malignant cells and cytotoxic effectors, precisely the configuration
our control framework seeks
to exploit (§3.5–§3.7). These are associations in a modest cohort (n = 15); our use of
the authors' annotations makes the myCAF/iCAF positive control confirmatory by
construction, but the immune and vascular compositional findings are independent of it.

### 3.4 The myCAF barrier controls the tumor only under specified conditions
Having implemented myCAF as a physical barrier—confinement, mechanical pressure, and
impaired drug delivery—alongside its immune-exclusion effect (Methods §2.4), we asked
when preserving or modulating the stroma, rather than depleting it, best controls the
tumor. Under a fixed sub-maximal treatment, we swept the myCAF barrier level (via
`k_caf_activate`, the target of anti-fibrotic modulation) and located the stromal level
that minimized tumor burden, across a grid spanning confinement strength (`caf_pressure`)
and immune-exclusion strength (`cd8_barrier_alpha`) (Fig. 4).

The optimal stromal level was regime-dependent. On a denser 6×6 grid with five seeds per
cell, we report the phase boundary probabilistically—the fraction of seeds for which
keeping stroma is optimal, P(keep-stroma)—rather than as a hard line. Where confinement was
strong and immune exclusion weak, the optimum was interior (a non-zero stromal level) and
this was robust across seeds (P(keep) ≥ 0.8), with preservation reducing tumor burden by up
to ~0.5-fold relative to full stromal depletion—the myCAF-as-controllable-resource regime.
Where immune exclusion dominated, the optimum collapsed to zero (P(keep) ≤ 0.2; full
stromal reduction was best), because the immune- and drug-shielding costs of stroma
outweighed its confinement benefit. A transition band (P(keep) ≈ 0.4–0.6) separates the
two, quantifying where the strategy is genuinely uncertain. Across the grid the keep-stroma
regime held in 20 of 36 sampled conditions at P(keep) ≥ 0.6 (Fig. 4a), with the benefit of
keeping stroma concentrated in the low-immune-exclusion region (Fig. 4b).

This reframes the central question from "deplete or preserve the stroma?" to "what is the
optimal stromal state, and when?" It also formalizes the CAF-biology debate: the same
peritumoral myCAF barrier can be net-protective or net-harmful depending on the local
balance, so both the failure of stromal ablation [6,7] and reports that CAF reduction
improves immunotherapy emerge as different regimes of one trade-off. We stress that the
phase boundary depends on parameters that are mechanistically motivated but not
data-fitted; where real PDAC lies on it is an empirical question this framework poses
rather than answers (Discussion §4.5).

**Robustness to restored CAF pro-tumor biology.** Because a baseline that omits the
tumor-supporting roles of CAFs could bias this map toward "keep stroma," we recomputed it
with two such axes switched on—a paracrine proliferation boost and a signaling-based
drug-tolerance term by which local myCAF attenuates the drug's effect on sensitive cells
(IL-6/JAK–STAT-type therapy-induced resistance; Methods §2.4, Fig. S11). Restoring this
biology shrank the resource regime modestly (7/16 → 6/16 sampled conditions) and reduced
the average benefit of keeping stroma, but did not abolish it: preservation still won in
the strong-confinement, weak-immune-exclusion corner. The condition-dependent principle
is therefore qualitatively robust to the previously omitted pro-tumor biology, while the
favorable region is somewhat narrower (Discussion §4.5).

**Application—immune priming is gated by the stromal barrier.** The phase map predicts
that immune-directed agents should help only where the stroma does not exclude T cells.
We illustrate this with an immune-priming agent (an illustrative immune-priming agent
under a strongly immune-excluding stromal regime; a clinical telomerase-vaccine instance
and the trial evidence motivating it are given in the Supplement, Fig. S4). Encoded as an
immune-priming perturbation in a strongly immune-excluding regime, the agent alone barely
reduced tumor burden (1.89× vs 2.28× untreated)—the dense barrier excluded even the
boosted T cells—whereas co-administration with a barrier-opening anti-fibrotic agent
controlled the tumor (0.76×); the anti-fibrotic alone was insufficient (1.66×)
(Fig. S4). Notably, sustained co-administration outperformed a pulsed
"open-then-rest" schedule (1.49×), because letting the barrier re-form between pulses
re-excluded the T cells—so, unlike drug holidays that exploit competition among resistant
clones, immune access must be held open. This is a concrete instance of the framework's
message: an immune-priming agent's value is contingent on the stromal state, and pairing
it with barrier-opening is predicted to be necessary specifically in the
immune-exclusion-dominant regime.

### 3.5 Adaptive scheduling achieves simulated coexistence at a fraction of the modeled exposure of continuous dosing
On a controllable-but-not-eradicable tumor, we compared no treatment, continuous
full-intensity therapy (a within-model reference arm, not the clinical standard of care),
and adaptive on/off dosing (Fig. 5). Continuous dosing drove the
sensitive population to extinction (final burden 0.00× baseline) at the highest
cumulative exposure (120 units), leaving a purely resistant residue. Adaptive dosing
held the tumor low (0.16× baseline) at roughly one-fifth the exposure (22 units), with a
modest resistant fraction (0.14); both treated arms remained below the progression
threshold over 150 days, whereas the untreated tumor progressed at day 127. Rather than
collapse these into a single composite score (a naive score that rewards low exposure
would rank the zero-exposure untreated arm highest, which is meaningless), we evaluate
outcomes as a **progression-constrained objective**: among strategies that keep the tumor
progression-free, prefer the one with the lowest cumulative exposure, with resistant
fraction and burden-AUC as secondary read-outs. Under this constraint adaptive scheduling
dominates continuous dosing (both progression-free, but adaptive at ~1/5 the exposure);
the untreated arm is excluded because it fails the progression constraint. This
progression-constrained, Pareto-based evaluation (developed fully in §3.6, Fig. S10)
replaces any single composite metric throughout.

**Separating schedule from agent (fair 2×2).** To avoid conflating adaptive scheduling
with the choice of agent, we ran a factorial comparison—{gemcitabine, the lead natural
combination} × {continuous, adaptive}, with compound synergy disabled (Fig. S5). All four
arms controlled the tumor (final burden 0.00–0.17×); they differed chiefly in cumulative
exposure. Adaptive scheduling was the robust, agent-independent lever: it cut exposure
~6–7-fold for both agents (gemcitabine 128→20; natural 52→7) at a minor increase in
residual tumor. The natural combination's lower exposure than gemcitabine at matched
schedule (52 vs 128 continuous; 7 vs 20 adaptive) follows largely from the assumed
compound toxicity weights rather than a model-discovered efficacy difference, since tumor
control was comparable across agents. The large gap between natural-adaptive and
continuous-gemcitabine therefore reflects the schedule effect plus the exposure weighting,
not agent efficacy per se—an important caveat for interpreting §3.6.

**Schedule effect within a single agent.** To isolate the schedule effect from any
compound assumption, we varied only the *schedule* of one encoded agent—gemcitabine at a
fixed efficacy coefficient—across a continuous full-intensity arm (a model reference, not
the clinical standard), an intermittent clinical-approximation arm (three weeks on, one
week off, q28), and burden-triggered adaptive dosing (Fig. S13). All arms held the tumor
below the progression threshold in every seed (final burden ~0.00×); they differed only
in cumulative modeled exposure, which fell monotonically with less-continuous scheduling
(continuous 128 → intermittent q28 98 → adaptive 18, 20 seeds). The message is
schedule-driven and agent-agnostic: even for the same encoded agent at the same efficacy,
adaptive scheduling can maintain control at a substantially lower modeled exposure than
continuous full-intensity dosing. This is a within-model comparison of dosing strategies,
not a claim that any schedule is clinically superior.

**Does the adaptive advantage survive realistic monitoring?** The adaptive arms above
assume an idealized controller that observes true tumor burden instantaneously and without
error. In the clinic, feedback would come from a serum biomarker such as CA19-9, sampled
at intervals, with measurement noise, a biomarker-to-burden lag, a non-informative
non-secretor subgroup, minimum treatment durations, and a safety rule against premature
discontinuation. We re-ran the adaptive controller under such an observation model
(Fig. S12): a latent CA19-9 signal lagging true burden, read every 28 days with 25%
log-normal noise, de-escalating only after two consecutive confirmed low reads and a
minimum on-duration, and forced back on above a safety ceiling. The adaptive advantage
survived but attenuated: control was still maintained in every seed, and exposure remained
well below continuous (median 67 vs 120), but roughly half the idealized saving was lost
to observation delay and noise (ideal 21 → observed 67). More frequent sampling recovered
more of the benefit (14-day interval → 45; 56-day → 90), and Lewis-antigen-negative
non-secretors—for whom CA19-9 is uninformative—defaulted to continuous therapy and lost
the adaptive benefit entirely. Thus the schedule advantage is real but monitoring-limited,
and its clinical realization would depend on a usable, sufficiently frequent readout.

### 3.6 Food-medicine-homology regimens control tumor burden at low modeled exposure, but only combinations do so robustly
As an illustrative case study of how the framework prioritizes candidates—not an efficacy
claim—we evaluated single agents and combinations under adaptive scheduling against a
continuous-gemcitabine reference. We present the **multi-seed (n = 30) analysis as the
primary result**; a single-seed exploratory ranking is provided only in the Supplement
(Fig. S8) and is superseded here, because single-seed orderings proved unreliable (below).

**Multi-seed, progression-constrained Pareto view (n = 30).** We re-ran the regimens over
30 seeds (varying both tissue architecture and simulation stochasticity, with combination
synergy switched off) and summarized the progression-versus-exposure trade-off as a Pareto
frontier with rank-stability statistics (Fig. S10). Only the
adaptive *combinations* controlled the tumor in every seed—curcumin + garlic + Rg3
(progression-free in 30/30; exposure 6, 95% interval [5–7]), garlic + mugwort (30/30;
7 [7–8]), and ginseng + garlic + mugwort (30/30; 8 [7–9])—together with continuous
gemcitabine (30/30) but at ~20× their exposure (128). Single garlic gave the lowest
exposure of all (3 [3–4]) and the best median rank, but was less reliable (progression-
free in 28/30). Critically, several single agents that scored well on one seed were *not*
robust: wild ginseng controlled the tumor in only 1/30 seeds, and mugwort and curcumin in
none, so the single-seed exploratory ordering (Fig. S8) does not survive replication. The
anti-fibrotic pair danshen + astragaloside likewise failed to control across seeds
(median final burden 1.23×), consistent with anti-fibrotic action helping only in the
regimes the phase map identifies (§3.4) rather than universally. The robust, defensible
reading is therefore narrow: adaptively scheduled low-exposure *combinations* can control
this simulated tumor at a small fraction of the modeled exposure of continuous
gemcitabine, whereas rankings among individual low-exposure agents are seed-sensitive and
should not be over-interpreted. We carry the most robust combination (curcumin + garlic +
Rg3, progression-free in 30/30 at the lowest combination exposure) forward as the
illustrative lead for dose optimization (§3.7). This multi-seed, progression-constrained
Pareto analysis is the primary compound result; the single-seed ordering of Fig. S8 is
exploratory only.

**Epistemic uncertainty dominates: control collapses when compound assumptions are
sampled.** The 30-seed analysis varies architecture and stochasticity but holds the
*compound assumptions* fixed—yet those coefficients are assigned, not measured. We
therefore ran a Monte Carlo over the compound assumptions themselves, on the
API-resolved regimens (§3.8): for each of 100 draws we sampled every compound's effect
coefficient (log-normal), its exposure weight (log-normal), its oral bioavailability
(from API-specific priors grounded in Table S2—high for SAC, very low for curcumin—used
to scale the effect), and the combination synergy (Uniform[0, 0.3]), then re-ranked all
regimens under those shared assumptions (Fig. S14). The result sharply qualifies the
nominal-parameter picture. Whereas at nominal coefficients the combinations and
20(S)-Rg3 control the tumor, under sampled uncertainty the natural-API regimens control
in only 0–2% of draws, while gemcitabine—whose pharmacology is defined and whose
bioavailability prior is near-certain—controls in 96%. The collapse is driven chiefly by
bioavailability: an orally administered natural product whose systemic/target exposure is
low and uncertain cannot reliably deliver its encoded effect, and no number of
architecture seeds repairs that. The honest conclusion is thus stronger than "rankings
are seed-sensitive": the *control benefit* of the natural-compound regimens is not robust
to epistemic uncertainty in their pharmacology, and their apparent advantage over a
defined cytotoxic is contingent on assumptions the model cannot currently justify. This
is the sharpest statement of why these outputs are hypotheses to be tested—by measuring
the effect and exposure parameters—rather than rankings to be trusted.

### 3.7 Dose and drug-holiday optimization further lowers modeled exposure (illustrative)
For the lead combination (curcumin + garlic + ginsenoside-Rg3), a seed-averaged grid
search over dose intensity and off-threshold under adaptive scheduling (200 days;
Fig. S9) found an optimum at a 40% simulated dose fraction with an off-threshold of 0.4,
achieving control at cumulative exposure 17 (resistant 0.04, final 0.8× baseline). This
40% is a model-parameter fraction under the tested parameterization, not a clinical dose. This markedly undercut
full-dose adaptive scheduling (exposure 27 at comparable control), indicating that a
sub-maximal, adaptively timed dose preserves control while further trimming exposure—
reinforcing that, under a control objective, less drug delivered adaptively can be
better.

### 3.8 Multiscale molecular grounding links compounds to model parameters under a transparent evidence hierarchy
Each compound is annotated with its principal target(s) and an explicit evidence tier—
experimental co-crystal, docking prediction, or pathway-level mechanistic inference
(Fig. 6). Only the two conventional agents are supported by experimental co-crystals in
which the compound itself is resolved—gemcitabine–deoxycytidine kinase (PDB 1P62) and
erlotinib–EGFR kinase domain (PDB 1M17). A minority rest on docking (e.g.,
curcumin–STAT3 1BG1; the entecavir–KDM5B repurposing hypothesis 5A1F), where the
displayed ligand is a reference rather than our compound. Most food-medicine-homology
anti-fibrotic compounds act through the TGF-β/Smad and NF-κB axes for which no
compound-bound structure exists and are represented at the mechanistic tier only. Each
target maps to the ABM parameter class it perturbs—anti-proliferative to proliferation
rate, anti-CAF/anti-fibrotic (predominantly TGF-β/Smad) to CAF activation, and
immunomodulatory to CD8 recruitment/killing—providing the molecule-to-cell-to-tissue
link. This hierarchy is a transparency device, not a claim of binding validation: the
mechanistic-tier compounds that dominate our low modeled-exposure regimens are grounded
in pathway pharmacology, and their molecular engagement remains a hypothesis for
experimental confirmation.

Encoding a food such as garlic as one perturbation also conflates dozens of
constituents of variable content and stability, so a distinct food-medicine entity is
not yet a defined molecular species. Where a standardized, bioavailable active
ingredient exists, the framework can be re-grounded at that resolution. We illustrate
this with S-allylcysteine (SAC), the most abundant and chemically stable organosulfur
constituent of aged garlic extract—odorless and orally bioavailable, unlike the labile
allicin that gives raw garlic its acute effects [25]. SAC's characterized activity is
anti-fibrotic: suppression of TGF-β/SMAD3 and STAT3 signaling with reduced α-SMA and
collagen deposition in hepatic [23] and pulmonary [24] myofibroblasts. Accordingly, and
unlike the coarse "garlic" entry (encoded on the anti-proliferative axis), SAC is
resolved onto the CAF-activation (barrier-modulating) parameter—placing it in the
anti-fibrotic arm that the phase map (§3.4) identifies as beneficial precisely in the
immune-exclusion-dominant regime, and as counter-productive where the barrier is
containment-dominant. Two caveats keep this honest: SAC's anti-fibrotic evidence is from
liver and lung, so its direct action on PDAC myCAF is an untested hypothesis (a natural
target for the organoid/CAF co-culture assay of §4.6); and its net benefit is
regime-conditional rather than universal. SAC is offered as a concrete template for the
active-ingredient resolution the framework supports, not as a validated PDAC therapeutic.

**From food entities to defined drug products.** SAC is not the only such resolution.
Each food-level entity in our regimens can, where a standardized active ingredient exists,
be re-grounded onto a defined molecular species with a documented development history
(Table S2): mugwort (*Artemisia argyi*) → eupatilin, a flavonoid that is the standardized
active ingredient of a marketed gastroprotective drug and suppresses stellate-cell/
myofibroblast activation via β-catenin/PAI-1 [26], mapped—like SAC—onto the anti-fibrotic
CAF axis; and ginseng → 20(S)-ginsenoside Rg3, a single stereoisomer (distinct in activity
and oral pharmacokinetics from its 20(R) epimer [27]) that is standardized as a clinical
capsule product [28]. This matters because a food and a drug product differ on axes our
parameter perturbation does not by itself capture—active ingredient, extraction method,
standardized content, purity, stereochemistry, formulation, oral bioavailability, tissue
distribution, and batch-to-batch consistency. Table S2 makes these explicit for each
resolved API and, crucially, marks where they remain undefined: curcumin, for instance, is
chemically well-defined but has notoriously poor and formulation-dependent oral
bioavailability [29], so its encoded effect presumes an exposure that a conventional
preparation may not achieve. The framework therefore resolves *which* molecule and *which*
mechanistic axis, but the pharmacokinetic and chemistry-manufacturing-controls layer that
turns a molecule into a drug product is out of its current scope and is the explicit
subject of Table S2 and §4.6.

### 3.9 The control framework runs directly on real patient tissue as its initial condition
The mechanistic sweeps above use a single synthetic architecture so that one variable can
be changed at a time; to show that the pipeline is not confined to idealized geometry, we
seeded the ABM directly from real CosMx tissue. From the SCOTIA cohort we cropped a
1500-µm tumor-centered window from a treatment-naive patient (U7-a) and a
chemoradiation-treated patient (T4-a) and used each cell's measured position and type as
the model's initial condition (Fig. S6). The two tissues are architecturally distinct:
the naive window is tumor- and myCAF-dominated (1,418 malignant, 775 myCAF, 184 iCAF),
whereas the treated window is depleted of tumor and shifted toward an
iCAF/macrophage-rich stroma (236 malignant, 936 iCAF, 748 macrophage, 580 myCAF)—the
CAF-subtype and myeloid remodeling expected after therapy. Running the same untreated
versus adaptive-natural-combo contrast on each real tissue reproduced the qualitative
control result—the adaptive regimen held tumor burden below the naive tissue's untreated
trajectory—while the quantitative dynamics differed between the two architectures,
confirming that patient-specific structure, not a hand-built lattice, drives the outcome.
This grounds the framework's initial conditions in measured tissue; the parameter sweeps
themselves remain synthetic, and patient-specific *calibration* of dynamics is future work
(§4.6).

### 3.10 The model makes non-trivial spatial predictions, not just restatements of its rules
A fair objection is that once myCAF is coded to confine daughter cells and block drug and
immune access, "confinement matters" is built in. We therefore asked whether the model
yields predictions that its local rules do not make obvious—predictions that depend on
*emergent* spatial structure rather than on any single encoded rule.

First, spatial *geometry* changes therapeutic outcome at *fixed* cell counts. Comparing
the contained architecture (myCAF in peritumoral rings) with a diffuse one (identical
myCAF, iCAF, CD8, and macrophage counts, scattered) under matched sub-maximal therapy,
every modality was less effective against the confining geometry, and the effect was
modality-specific: immune-directed therapy was ~12-fold weaker in the contained tissue
(residual 0.73× vs 0.06× in diffuse) because the ring barrier physically excludes CD8,
whereas cytotoxic therapy was less affected (Fig. S15a). Abundance alone is therefore not
predictive—the same myCAF, arranged as a barrier versus scattered, changes the value of
immunotherapy by an order of magnitude, so a measurable descriptor of *arrangement*
(not amount) is what would stratify immune-therapy benefit. Second, and relatedly, the
model predicts *when* reducing CAF widens the invasion front: anti-fibrotic treatment
increased tumor spatial spread more in the confining geometry than in the diffuse one
(90th-percentile radius +78 µm vs +51 µm; Fig. S15b), i.e., releasing a barrier that was
actually containing the tumor—an emergent, spatially explicit trade-off, not a coded rule.
Third, treatment *order* matters in a barrier-gated regime: with a strongly
immune-excluding stroma, immune-alone and anti-fibrotic-alone both failed, combining them
in any order roughly halved residual burden, and opening the barrier before adding
immunotherapy was marginally better than the reverse (Fig. S16)—a modest but
directionally consistent sequence effect (using both mechanisms matters more than their
order). Fourth, these geometry effects are captured by a *measurable* spatial biomarker,
not by abundance: across a panel of architectures (contained tissue thinned to varying
degrees, plus diffuse), the peritumoral myCAF density (myCAF within 30 µm of tumor, per
tumor cell) predicted immune-therapy efficacy strongly (|r| = 0.92), whereas total myCAF
abundance predicted it only weakly (|r| = 0.56) and failed on the decisive contrast—the
contained and diffuse tissues have *identical* myCAF abundance yet opposite immune
outcomes, differing 7-fold in peritumoral density (Fig. S17). The actionable, testable
prediction is thus a patient-stratification hypothesis: a peritumoral (rim) myCAF-density
readout, not a bulk stromal-content measure, should identify the tumors in which
immune-directed therapy is barrier-limited and would benefit from stromal opening.
Together these are predictions about emergent geometry, timing, and measurable structure
that motivate spatial biomarkers and scheduling hypotheses, moving the framework beyond a
restatement of its inputs.

Finally, to show that spatial structure is *necessary*—not incidental—we ran an ablation
series that isolates each ingredient and asks the decisive question: does spatial geometry
change the *selected strategy* (keep vs reduce stroma), or only the tumor number? (i) A
**well-mixed (mean-field)** model with the same mechanisms but with myCAF acting only
through its aggregate fraction (immune attenuation ∝ exp(−α·M)), no geometry or local
confinement, is blind to the effect that drives our results: because contained and diffuse
tissues share the same abundance it predicts *identical* outcomes, whereas the ABM predicts
an ~7-fold difference (Fig. S18a); and lacking any confinement benefit, its optimum is
*always* full stromal reduction (M = 0) at every immune-exclusion level, so it never
produces the keep-stroma regime (Fig. S18b). The mean-field is thus the "no-confinement"
ablation, and it inverts the strategy. (ii) **Randomized geometry at matched cell counts**
(contained vs diffuse) separates abundance from arrangement: same myCAF amount, ~12-fold
different immune-therapy outcome (§3.10 above, Fig. S15). (iii) The **immune-exclusion
contribution** is isolated by the phase map's vertical axis: increasing `cd8_barrier_alpha`
alone flips the optimum from keep to reduce (Fig. 4a), so removing that cost would leave
only the confinement benefit and favor preservation. Together these ablations show that
both spatial ingredients are load-bearing and, crucially, that geometry changes the
*selected stromal strategy*—not merely the numerical burden—which a non-spatial model
cannot reproduce.

---

## 4. Discussion

### 4.1 A framework that bridges three disconnected research lines
We present an in-silico framework that connects three individually mature but mutually
disconnected bodies of work—static network pharmacology of food-medicine-homology
compounds, wet-laboratory studies of herbal modulation of PDAC CAFs, and spatial/
agent-based modeling of the tumor microenvironment. Our systematic survey quantified
this gap: single-axis literature is abundant, yet queries requiring all three axes
returned only off-target hits (§3.1, Fig. 1). The contribution is therefore
integrative—a spatially grounded, adaptive-control model in which food-medicine-
homology compounds are encoded as mechanistic perturbations and optimized against a
coexistence objective. Crucially, the model implements the myCAF barrier physically,
with an explicit trade-off between its confinement of the tumor and its
immunosuppressive and drug-shielding costs; rather than assuming the stroma is a
resource, it maps the conditions under which preserving/modulating versus reducing it
best controls the tumor (§3.4, Fig. 4).

### 4.2 Biological plausibility
Several results align with independent biology. On author-annotated CosMx data, our
metric reproduced the myCAF-proximal/iCAF-distal architecture in 15 of 16 tumors and,
at the patient level, showed a rim-enriched myCAF barrier coincident with cytotoxic
T-cell exclusion in *both* untreated and CRT specimens (§3.3); the untreated-versus-CRT
differences trend as expected but are not statistically significant in this small cohort
(all BH-adjusted p > 0.19), so we lean only on the pattern common to both groups. We are
careful in interpreting even that: a peritumoral myCAF rim with persistent T-cell
exclusion is equally consistent with an *immunosuppressive* barrier as with a
tumor-*containing* one—the observation is a necessary, not sufficient, condition for
containment. This ambiguity is precisely
what our phase analysis formalizes: whether such a barrier aids or harms control
depends on the balance of physical confinement against immunosuppression, and both
the paradoxical worsening of PDAC after stromal ablation [6,7] and reports that CAF
depletion can improve immunotherapy response fall out as different regimes of the same
trade-off. In simulation, adaptive scheduling reduced tumor burden at lower cumulative
exposure than continuous dosing under the tested conditions (§3.5), and our sensitivity
analysis showed that the resistance fitness cost governs competitive suppression of
resistant clones (§S3), mirroring the adaptive-therapy rationale [9–12].

### 4.3 Positioning relative to prior art
The framework advances each of the three neighboring camps on the axis it lacks:

| Prior camp | What it provides | What it lacks | Our addition |
|---|---|---|---|
| Static network pharmacology of food-medicine-homology compounds [14,16] | compound–target–pathway maps | time, space, scheduling | dynamic, spatial, adaptive-control simulation |
| Wet-lab herbal–CAF studies in PDAC [17,18] | real CAF-modulating biology | prediction/optimization | in-silico ranking of combinations, sequences, doses |
| Spatial/ABM models of PDAC TME [19,20] | tumor–stroma–immune dynamics | natural compounds; coexistence objective | food-medicine-homology perturbations + control objective |

Crucially, we do not claim novelty in *identifying* molecular targets of these
compounds—that space is saturated (§3.1)—but in embedding that knowledge in a dynamic,
spatial, control-theoretic framework.

### 4.4 Translational implications
Two design principles emerge, both stated as model-generated hypotheses. First, the
stromal target is not depletion but the optimal stromal state: in the regime where
confinement outweighs immunosuppression, partial modulation—rather than maximal
anti-fibrotic pressure—best controls the tumor, whereas in the immunosuppression-
dominant regime stromal reduction is favored (§3.4). Which regime a given tumor occupies
is an empirical question our framework poses but cannot yet answer. Second, within a
control objective, pairing a low-exposure anti-CAF/immunomodulatory backbone with a
resistance-sparing anti-proliferative partner controlled tumor burden at a fraction of
the *modeled treatment-exposure* of continuous gemcitabine under assumed compound
weights (§3.6–§3.7)—a statement about exposure burden in the model, not a clinical
toxicity prediction. For a disease in which conventional therapy imposes severe toxicity
for limited benefit, a strategy oriented toward slowing proliferation and preserving
quality of life—rather than eradication—may be a more realistic near-term goal, but only
if efficacy and the assumed low-toxicity profile are confirmed experimentally and
clinically.

**Scope—disease setting.** Our model represents a single tumor focus contained by a
myCAF barrier, a geometry that pertains to localized and locally advanced PDAC
(stages I–III). In this setting, slowing and containing the primary tumor plausibly
delays local mass-effect complications—biliary and pancreatic-duct obstruction,
duodenal/gastric-outlet compression, and vascular encasement—benefits aligned with a
quality-of-life and functional-survival objective. The rationale weakens for
metastatic disease (stage IV): disseminated cells have, by definition, escaped the
local barrier, and metastatic-site stroma (e.g., hepatic) differs biologically. In
particular, malignant ascites in peritoneal carcinomatosis arises from peritoneal
vascular permeability, subperitoneal lymphatic obstruction, and portal
hypertension—mechanisms a primary-tumor stromal-control framework does not address and
should not be expected to prevent. The predictions here therefore apply to local
control and complication palliation, not to systemic or peritoneal metastatic disease.

### 4.5 Limitations
Our study has important limitations that bound its claims. First and foremost, it is
entirely in silico: outputs are ranked, testable hypotheses, not evidence of efficacy,
and require experimental validation. Second, the ABM is phenomenological with
literature-inferred parameters rather than data-fitted; although the adaptive advantage
was robust to ±50% parameter variation (§S3), outcomes were most sensitive to tumor
proliferation and immune-barrier parameters, which therefore require empirical
calibration, and control was not guaranteed under the most aggressive proliferation
regimes. A variance-based global (Sobol) analysis over eight parameters concurs and
sharpens this: tumor proliferation dominates the control outcome (total-order
$S_T\approx0.93$), the immune-exclusion barrier is the next most influential and acts
largely through interactions ($S_T\approx0.43\gg S_1\approx0.12$), and—consistent with
the confinement-versus-immunosuppression framing—the physical-confinement strength and
the resistance parameters have negligible global influence (Fig. S7); the wide
confidence intervals at this screening sample size mean the parameter *ranking*, not the
precise indices, is the robust claim. Relatedly, the confinement-versus-immunosuppression trade-off that decides
whether stroma should be preserved or reduced is governed by parameters (stromal
confinement strength, drug shielding, and immune exclusion) that are mechanistically
motivated but not data-fitted; the phase boundary is therefore qualitative, and where
real PDAC falls on it is unresolved. Third, our spatial analyses are of 2D sections; true tissue architecture is
three-dimensional, and the myCAF-proximity positive control failed on our targeted-
panel Xenium annotation, succeeding only with full-panel author annotations (§3.2)—so
real-data spatial claims rest on the CosMx cohort and remain associations in a modest
sample (n = 15). Fourth, the mechanistic parameter sweeps used a single synthetic
architecture (chosen so one variable could be varied at a time); we confirmed the
pipeline runs on real CosMx patient tissue as its initial condition (§3.9, Fig. S6), but
have not yet calibrated its dynamics to patient-specific data, and the toxicity scale
remains in arbitrary units that are not clinically calibrated. Fifth, even
low modeled-exposure natural combinations can select for resistance under aggressive
suppression (competitive release, seen in aggressive-suppression regimes). Finally, the novelty survey is PubMed-
restricted, and the molecular grounding for most food-medicine-homology compounds rests
on pathway-level mechanistic inference rather than solved compound-bound structures
(§3.8). Relatedly, several entries are still encoded at the food rather than the
molecule level; only where a standardized active ingredient is available (e.g., SAC for
garlic) have we begun resolving a food to a defined, bioavailable species, and even
there the PDAC-myCAF action is hypothetical (§3.8).

Beyond these, several structural assumptions bound the biology the model can speak to,
and we state them plainly. **(i) Confinement is largely encoded, not discovered.** The
physical rules by which local myCAF density blocks daughter-cell placement and lowers
carrying capacity are model inputs; the phase behavior we report (an intermediate
stromal optimum where confinement dominates) is therefore a consequence of those rules,
not independent evidence about human tissue. The defensible claim is conditional—*if*
myCAF exerts measurable physical confinement, an intermediate stromal optimum can
emerge—whereas the stronger claim, that human PDAC myCAF actually contains the tumor to
the predicted degree, is untested and would require relating myCAF/ECM density to
invasion fronts, tumor budding, ductal spread, perineural invasion, and dissemination in
real tissue. **(ii) The pro-tumor assumption is now tested, not assumed away.** Our baseline set the
direct CAF pro-proliferative term to zero, which risks biasing the model toward
containment; real CAFs also support tumors through paracrine cytokines and growth
factors, metabolic support, ECM remodeling, invasion facilitation, and therapy-induced
survival signaling. The very CosMx cohort we use (Shiau et al. [20]) reports
CAF–malignant interaction changes and IL-6-family signaling linked to therapy resistance,
and shikonin reversal of CAF-induced gemcitabine resistance is documented [18]. We
therefore added two axes—a modest paracrine proliferation boost (`caf_protumor`) and a
signaling-based drug-tolerance term (`caf_survival`, in which local myCAF attenuates the
drug's effect on sensitive cells independently of physical penetration)—and recomputed
the confinement × immune-exclusion phase map with this biology switched on (§3.4,
Fig. S11). Restoring it shrinks the "stroma-as-resource" regime (from 7/16 to 6/16
sampled cells) and lowers the average benefit of keeping stroma, but does not erase the
regime: preserving/modulating stroma still aids control where physical confinement is
strong and immune exclusion is weak. The central claim is thus qualitatively robust to
the previously omitted biology, though the favorable region is narrower; full calibration
of these axes against the IL-6/JAK–STAT data remains a priority.
**(iii) CAF states are fixed rather than plastic.** myCAF/iCAF/apCAF are modeled as
stable types, whereas real CAFs interconvert (TGF-β vs IL-1/JAK–STAT balance,
therapy-induced transitions, LRRC15⁺ and senescent/inflammatory substates); a claim of
CAF *reprogramming* rather than depletion ultimately needs a subtype-transition model.
**(iv) Immunity is reduced to CD8 T cells.** Macrophages, MDSCs, Tregs, NK and dendritic
cells, and the vasculature are present only as static context; Fig. 3 shows treatment-
associated macrophage shifts that the control mechanism does not yet act on, a gap
between spatial observation and simulated mechanism. **(v) Resistance is not
PDAC-specific.** The resistance fitness cost is borrowed from NSCLC and the resistant
phenotype is a fixed binary decoupled from CAF state, whereas gemcitabine resistance is
mixed genetic/epigenetic/metabolic/microenvironmental and is, in our own cited biology,
CAF-induced; a continuous, reversible drug-tolerance state coupled to CAF signaling
(with resensitization during holidays) would be more faithful. **(vi) The controller is
partly idealized, and outcomes are local.** Our primary adaptive arm switches on the true
simulated burden; to test clinical realizability we added a CA19-9-based observation model
with interval sampling, measurement noise, a biomarker-to-burden lag, non-secretors,
minimum on/off durations, and a confirm-before-stop safety rule (§3.5, Fig. S12). The
adaptive advantage survived but was roughly halved by observation delay and noise,
degraded with longer sampling intervals, and vanished in non-secretors—so the benefit is
real but monitoring-limited, and a deployable controller would still need a validated,
sufficiently frequent readout (imaging, CA19-9, or ctDNA) and formal control-theoretic
analysis. Moreover, the model controls a
single primary focus, and local burden is not a validated surrogate for survival:
clinically meaningful endpoints (local and distant progression-free survival, biliary-
obstruction-free and pain-free survival, conversion to resectability, overall survival,
and patient-reported quality of life) depend on systemic disease the model does not
represent. **(vii) Stochastic and calibration limits.** Most exploratory results use three-seed
means; for the central compound comparison we therefore re-ran a 30-seed analysis with
medians, 95% intervals, and rank-stability box plots, summarized as a
progression-versus-exposure Pareto frontier (§3.6, Fig. S10). That re-analysis both
supports a narrow robust claim (adaptive combinations control at a fraction of
gemcitabine's exposure) and overturns several single-agent rankings that did not survive
replication—so single-seed orderings elsewhere in the paper should be read as
provisional. Crucially, that analysis addresses only *stochastic* (architecture/seed)
uncertainty; a Monte Carlo over the *epistemic* uncertainty in the compound assumptions
themselves—effect coefficients, exposure weights, oral bioavailability, and synergy
(§3.6, Fig. S14)—shows the natural-compound control benefit is *not* robust: under sampled
assumptions those regimens control in 0–2% of draws versus 96% for gemcitabine, the
collapse being driven by bioavailability uncertainty. The compound rankings are therefore
hypotheses contingent on unmeasured pharmacology, not dependable predictions. The Sobol
analysis remains screening-level (§S7). Calibration and validation
are still not separated (no held-out patient cohort, posterior parameter distributions,
or identifiability analysis), and for the same reason we describe the patient-tissue runs
as *patient-geometry initialized* rather than patient-calibrated (§3.9).

### 4.6 Future directions
Several steps follow directly. Experimentally, the ranked regimens define a low-cost
validation ladder—2D PDAC cell-line dose–response and combination/synergy assays, then
patient-derived organoid screening—prioritizing the anti-CAF-backbone combinations and
the competitive-release predictions. Computationally, priorities are: parameter
calibration against experimental rates (especially proliferation and immune-barrier
terms, which the global Sobol analysis flags as dominant, Fig. S7); a larger Sobol
sample with multi-seed replication to tighten the indices beyond the present
screening-level estimate; a 3D registration pipeline over serial sections to replace 2D
metrics; annotation transfer from full-transcriptome references to make CAF subtyping
robust; molecular-dynamics or larger-scale virtual screening to strengthen the
molecular tier; and systematically resolving the remaining food-level entities to
standardized active ingredients with defined pharmacokinetics, as begun for garlic→SAC,
so that exposure and toxicity are grounded in real compound doses rather than
assumed weights. Broadening the organ context (e.g., the HCC "co-opted stroma" regime)
and extending the survey to Scopus/Web of Science under a PRISMA workflow would further
generalize and harden the framework.

### 4.7 Conclusion
Within a control—rather than eradication—objective, this work's central contribution is
a principle: by implementing the myCAF barrier as a physical structure with an explicit
confinement-versus-immunosuppression trade-off, it reframes the therapeutic question
from "deplete or preserve the stroma?" to "what is the optimal stromal state, and when?"
and maps the conditions that decide the answer. The food-medicine-homology compounds
serve as a case study demonstrating how the framework converts this principle into
concrete, ranked, testable hypotheses—candidate stromal targets, compound combinations,
and schedules—rather than as a claim that any compound treats PDAC. The result is a
reproducible path from real spatial data to condition-dependent stromal-control
principles and a prioritization engine for the experimental work that must follow.

---

## Figures

**Figure 1. Novelty landscape.** (a) Per-dimension PubMed hit counts from the targeted
survey (twelve queries: ten dimensional, two integrated); individual axes (blue) are well
populated while integrated, three-axis queries (red) return only off-target hits. (b)
Three-camp schematic: static network pharmacology, wet-lab herbal–CAF, and spatial/ABM
modeling are each populated but no directly matching study occupies their intersection in
our targeted search. *(assets/fig1_novelty.png)*

**Figure 2. The containment metric is sound; annotation is the limiting factor.**
(a) Synthetic contained (myCAF ring; barrier z = +22.3) vs diffuse (scattered; +0.9).
(b) On author-annotated CosMx data the myCAF-proximal positive control passes in 15/16
tumors (points left of 0 = correct). (c) The same metric fails on marker-annotated
Xenium data (1/5), localizing the failure to annotation. *(assets/fig2_validation.png)*

**Figure 3. Patient-level peritumoral composition, untreated (n=9) vs CRT (n=6) (CosMx
author annotations).** Per-cell-type rim (30 µm) enrichment z with the patient as the unit
of analysis; markers show group mean ± 95% CI (patient bootstrap) and small points show
individual patients. myCAF is tumor-adjacent and CD8⁺ T cells are excluded in *both*
groups. Untreated-versus-CRT differences are compared by Mann–Whitney U with
Benjamini–Hochberg correction across the ten cell types; none is significant (all adjusted
p > 0.19; myCAF p = 0.53), so the treatment contrasts are suggestive, non-significant
associations in this small, underpowered cross-sectional cohort (different patients, not
paired). *(assets/rim_scotia_stats.png)*

**Figure 4. The model predicts a non-zero optimal stromal level only under specified
conditions (6×6 grid, 5-seed phase uncertainty).** Under fixed sub-maximal treatment the
myCAF level was swept and the tumor-minimizing (optimal) level located, over a 6×6 grid of
confinement strength (caf_pressure) × immune-exclusion strength (cd8_barrier_alpha), with
5 seeds per cell. (a) Probability across seeds that keeping stroma is optimal,
P(keep-stroma)—a probabilistic phase boundary rather than a hard line: keeping stroma is
robust (P ≥ 0.8) in the low-immune-exclusion region, uncertain (P ≈ 0.4–0.6) in the
transition band, and disfavored (P ≤ 0.2) where immune exclusion is strong. (b) Mean
benefit of keeping stroma = tumor(no stroma) − tumor(optimal), annotated with half the
10–90% seed interval. The keep-stroma ("resource") regime occupies the low-immune-exclusion
region and holds in 20 of 36 sampled cells at P ≥ 0.6. *(assets/phase_dense.png)*

**Figure 5. Exposure-schedule comparison (within-model).** Tumor burden, resistant
fraction, and cumulative modeled exposure for untreated, continuous, and adaptive
dosing. Evaluated under the progression-constrained objective (§2.7): among
progression-free strategies, adaptive achieves simulated coexistence at ~1/5 the modeled
exposure of continuous, while the untreated arm is excluded for failing the progression
constraint. "Continuous full-intensity" is a within-model reference arm, not clinical
standard-of-care. *(assets/control_strategies.png)*

**Figure 6. Multiscale molecular grounding** with evidence tiers (experimental vs
docking vs mechanistic). *(assets/drug_structures_3d.png)*

**Figure S3. Sensitivity analysis.** (a) resistance_cost sweep—adaptive low-exposure
advantage is robust across the range; resistance_cost governs resistant-clone
suppression. (b) One-at-a-time ±50% tornado—outcome sensitivity concentrates in
tumor-immune balance parameters. *(assets/sensitivity.png)*

**Figure S4. The model predicts stromal gating of immune priming (illustrative).** In a
strongly immune-excluding regime, an illustrative immune-priming agent (encoded as immune
priming) was tested alone and combined with a barrier-opening anti-fibrotic agent. (a)
Tumor trajectories; (b) final tumor burden per arm. In the model, the agent alone is
nearly futile (dense stroma excludes the boosted T cells); opening the barrier first
unlocks it, and sustained co-administration outperforms a pulsed open-then-rest schedule.
The clinical motivation for this illustration is the telomerase (hTERT) peptide vaccine
GV1001, which failed to improve survival in unselected advanced PDAC (TeloVac [21]) but
did so in an eotaxin-high subgroup (KG4/2015 [22])—consistent with immune priming being
effective only where the stromal barrier does not exclude T cells. This is a mechanistic
illustration of barrier-gated immunotherapy within the model, not a clinical prediction
for GV1001. *(assets/gv1001_sequential.png)*

**Figure S5. Fair 2×2 comparison—agent × schedule (synergy off).** {Gemcitabine, lead
natural combination} × {continuous, adaptive}, five seeds. (a) Final tumor burden; (b)
cumulative exposure. All arms control the tumor; adaptive scheduling cuts exposure
~6–7-fold for both agents, and the natural agent's lower exposure at matched schedule
reflects its assigned exposure weight rather than superior tumor control.
*(assets/fair_2x2.png)*

**Figure S6. Patient-geometry initialized simulation—the framework runs on real SCOTIA
CosMx tissue.** A 1500-µm tumor-centered window from a treatment-naive patient (U7-a) and
a chemoradiation-treated patient (T4-a); each cell's measured position and type seeds the
ABM. (Left) patient tissue as the initial condition, colored by core cell type with
counts; the naive tissue is tumor/myCAF-dominated while the treated tissue is
iCAF/macrophage-shifted. (Right) untreated versus adaptive-natural-combo trajectories on
each real tissue: the adaptive regimen holds tumor burden below the untreated trajectory
in both, while the quantitative dynamics differ by architecture. In silico.
*(assets/patient_grounded.png)*

**Figure S7. Global (Sobol) sensitivity analysis complements the one-at-a-time sweep.**
Variance-based first-order ($S_1$, main effect) and total-order ($S_T$, including
interactions) Sobol indices for eight parameters spanning tumor kinetics, immunity,
resistance, and the two barrier axes (immune exclusion, physical confinement), computed
on the adaptive-therapy arm (Saltelli sampling, $D=8$, $N=32$; error bars are bootstrap
95% confidence intervals; runaway-growth samples were capped at a fixed population for
tractability, which only bounds already-progressing regimes and does not change their
progression classification). (a) Control score and (b) final tumor burden. Tumor
proliferation dominates both outcomes; the immune-exclusion barrier ($\alpha$) is the
next most influential for control and acts predominantly through interactions
($S_T\gg S_1$), whereas physical confinement and the resistance parameters have
negligible global influence. Given the wide intervals at this screening sample size, the
robust reading is the parameter ranking rather than the exact index values. In silico.
*(assets/sobol.png)*

**Figure S8. Single-seed exploratory regimen ranking (superseded).** Single agents and
combinations under adaptive scheduling on one seed, versus a continuous-gemcitabine
reference. **This single-seed ranking is exploratory only and is not used as a result:**
it is superseded by the multi-seed, progression-constrained Pareto analysis (Fig. S10),
which shows that several agents ranked highly here (e.g., wild ginseng, curcumin) do not
control the tumor on replication. Shown for transparency of the analysis history; in
silico, under assigned exposure weights. *(assets/natural_adaptive_optim.png)*

**Figure S9. Dose × drug-holiday optimization for the lead combination (illustrative).**
Seed-averaged grid search over simulated dose fraction and off-threshold under adaptive
scheduling; the optimum near a 40% simulated dose fraction is a model-parameter value
under the tested parameterization, not a clinical dose. *(assets/dose_band_optimization.png)*

**Figure S10. Multi-seed re-analysis (n = 30): progression-vs-exposure Pareto frontier
and ranking stability.** Each regimen was re-run over 30 seeds varying both tissue
architecture and simulation stochasticity, with combination synergy off. (a) Median final
tumor burden versus cumulative modeled exposure with 95% intervals on both axes; the
dashed line is the Pareto frontier (lower-left is better), and untreated (exposure 0, high
burden) sits on it trivially—illustrating why a single composite score that rewards low
exposure is misleading. (b) Distribution of per-seed rank (rank 1 = controlled at lowest
exposure) as box plots. Only the adaptive combinations (and continuous gemcitabine, at
~20× the exposure) control the tumor in every seed; several single agents that looked
strong on one seed (wild ginseng, curcumin, mugwort) are not robust, and the anti-fibrotic
pair fails across seeds. In silico, under assigned exposure weights.
*(assets/pareto_seeds.png)*

**Figure S11. The condition-dependent stromal-control result is robust to restoring CAF
pro-tumor biology.** The confinement (`caf_pressure`) × immune-exclusion
(`cd8_barrier_alpha`) phase map, recomputed with the tumor-supporting roles of CAFs (a)
switched off (baseline) and (b) switched on (a paracrine proliferation boost plus a
signaling-based drug-tolerance term, `caf_survival`, in which local myCAF attenuates the
drug's effect on sensitive cells—an IL-6/JAK–STAT-type therapy-induced resistance).
Color = benefit of keeping stroma (tumor with no stroma − tumor at the optimal stromal
level); cells labeled "keep" are the stroma-as-resource regime. Restoring pro-tumor
biology shrinks the resource regime (7/16 → 6/16 sampled cells) and lowers its average
benefit, but preservation still wins in the strong-confinement, weak-immune-exclusion
corner. Three seeds per cell; in silico. *(assets/phase_protumor.png)*

**Figure S12. The adaptive advantage under a realistic CA19-9 observation model.** The
adaptive controller re-run with feedback from a serum biomarker rather than true burden
(20 seeds). (a) Observation-model mechanics for one representative patient: true tumor
burden, the latent CA19-9 signal lagging it, interval-sampled noisy measurements
(points), and drug-on periods (shading); the on/off decision uses only the measurements,
subject to minimum durations and a confirm-before-stop rule. (b) Cumulative exposure by
arm: continuous, ideal adaptive (true-burden feedback), observed CA19-9 adaptive, and
observed with a non-secretor patient; control (fraction of seeds progression-free) is
annotated. The adaptive saving survives realistic observation (median exposure 67 vs 120
continuous) but is roughly halved relative to the ideal controller (21), and non-secretors
default to continuous therapy. (c) Measurement-interval sensitivity: exposure rises toward
continuous as the CA19-9 interval lengthens (14 → 28 → 56 days). In silico.
*(assets/obs_model.png)*

**Figure S13. Schedule effect within a single encoded agent (gemcitabine).** Varying only
the schedule of one agent at a fixed efficacy coefficient (20 seeds): continuous
full-intensity (a model reference, not the clinical standard), intermittent
clinical-approximation (three weeks on/one week off, q28), ideal adaptive, and observed
CA19-9 adaptive. (a) Cumulative modeled exposure and (b) final tumor burden. All schedules
maintain control (final ~0.00×, progression-free in every seed); they differ only in
exposure, which falls with less-continuous scheduling (128 → 98 → 71 → 18). The same
encoded agent under adaptive scheduling maintains control at substantially lower modeled
exposure than continuous full-intensity dosing—a schedule effect independent of the choice
of agent. In silico. *(assets/gem_schedule.png)*

**Figure S14. Epistemic (compound-assumption) uncertainty by Monte Carlo.** Instead of
architecture seeds alone, the compound *assumptions* are sampled: for each of 100 draws,
every API's effect coefficient (log-normal), exposure weight (log-normal), oral
bioavailability (API-specific Beta priors grounded in Table S2, used to scale the effect),
and combination synergy (Uniform[0, 0.3]) are drawn, and all API-resolved regimens are
re-ranked under those shared assumptions. (a) Per-draw rank distribution (rank 1 =
controlled at lowest exposure). (b) Fraction of draws in which each regimen controls the
tumor (progression-free). Under sampled uncertainty the natural-API regimens control in
only 0–2% of draws—despite controlling at nominal coefficients—whereas gemcitabine, with
defined pharmacology and a near-certain bioavailability prior, controls in 96%. The
collapse is driven chiefly by bioavailability uncertainty, showing that the natural
regimens' benefit is not robust to epistemic uncertainty in their pharmacology. In silico.
*(assets/mc_uncertainty.png)*

**Figure S15. Spatial geometry, not just abundance, changes the model's predictions.**
Contained (myCAF in peritumoral rings) versus diffuse (identical cell counts, scattered)
architectures, five seeds. (a) Final tumor burden by treatment modality at matched
sub-maximal dose: every modality is weaker against the confining geometry, and immune
therapy is ~12× weaker (residual 0.73× contained vs 0.06× diffuse) because the ring
barrier excludes CD8; the best combination (cytotoxic+immune) is robust, but the value of
immunotherapy is geometry-dependent at fixed abundance. (b) Tumor spatial spread
(90th-percentile radius) untreated versus after anti-fibrotic CAF reduction: reducing CAF
widens the invasion front more where the geometry confines the tumor (Δ +78 µm contained
vs +51 µm diffuse). In silico. *(assets/predict_geometry.png)*

**Figure S16. Treatment order in barrier-gated immunotherapy.** In a strongly
immune-excluding contained regime (five seeds), immune-alone and anti-fibrotic-alone both
fail; combining the two mechanisms roughly halves residual burden regardless of order, and
opening the barrier first (anti-fibrotic lead-in, then anti-fibrotic+immune) is marginally
better than immune-first. The sequence effect is modest relative to the benefit of using
both mechanisms; no fixed schedule fully controls this aggressive regime. In silico.
*(assets/predict_sequence.png)*

**Figure S17. A measurable spatial biomarker—peritumoral myCAF density, not abundance—
distinguishes the immune regime.** A panel of architectures (contained tissue with myCAF
progressively thinned, plus diffuse), five seeds; each tissue's spatial biomarkers are
measured and a fixed immune therapy is run. (a) Tumor burden under immune therapy versus
total myCAF abundance: abundance predicts weakly (|r| = 0.56) and fails on the decisive
contrast—contained and diffuse tissues share the same abundance but have opposite
outcomes. (b) The same versus peritumoral myCAF density (within 30 µm of tumor, per tumor
cell): a strong predictor (|r| = 0.92). (c) Biomarker map (peritumoral density × abundance)
colored by immune outcome: response tracks the rim-density axis, not the abundance axis
(green = immune-responsive, red = excluded). This yields a patient-stratification
hypothesis: a peritumoral myCAF-density readout should flag barrier-limited,
immune-excluded tumors. In silico. *(assets/predict_biomarker.png)*

**Figure S18. Why the spatial model is necessary: comparison to a well-mixed (mean-field)
model.** A mean-field model with the same mechanisms (proliferation, immune killing,
resistance, drug) but with myCAF acting only through its aggregate fraction M (immune
attenuation ∝ exp(−α·M)) and no geometry or local confinement. (a) Geometry-blindness:
because contained and diffuse tissues have identical myCAF abundance, the mean-field model
predicts identical immune-therapy outcomes, whereas the spatial ABM predicts a ~7-fold
difference. (b) Strategy inversion: sweeping the aggregate stromal level M under immune
therapy, the mean-field final burden increases monotonically with M at every
immune-exclusion strength α, so its optimum is always full stromal reduction (M = 0)—it
never produces the conditional keep-stroma regime that the spatial ABM finds (Fig. 4).
Spatial structure changes the selected strategy, not just the numbers. In silico.
*(assets/wellmixed_compare.png)*

---

## Declarations

**Data availability.** All datasets analyzed are from public repositories: Xenium
PDAC (GEO GSE274673), CosMx SMI PDAC (Mendeley Data doi:10.17632/kx6b69n3cb.1). All
analysis code and derived result tables are available at
https://github.com/kusi81/pdac-coexistence-control.

**Ethics.** This study used only publicly available, de-identified datasets; no new
human or animal experiments were performed. No ethical approval was required.

**Funding.** The author received no specific funding for this work.

**Competing interests.** The author declares no competing interests.

**Author contributions (CRediT).** Seung-Il Kim: conceptualization, methodology,
software, formal analysis, investigation, data curation, visualization,
writing – original draft, writing – review & editing.

**Acknowledgments.** [optional]

---

## References

1. Siegel RL, Kratzer TB, et al. Cancer statistics, 2025. *CA Cancer J Clin* 2025;75(1):10-45. PMID: 39817679. doi:10.3322/caac.21871
2. Hosein AN, Brekken RA, Maitra A. Pancreatic cancer stroma: an update on therapeutic targeting strategies. *Nat Rev Gastroenterol Hepatol* 2020;17(8):487-505. PMID: 32393771. doi:10.1038/s41575-020-0300-1
3. Öhlund D, Handly-Santana A, Biffi G, et al. Distinct populations of inflammatory fibroblasts and myofibroblasts in pancreatic cancer. *J Exp Med* 2017;214(3):579-596. PMID: 28232471. doi:10.1084/jem.20162024
4. Elyada E, Bolisetty M, et al. Cross-species single-cell analysis of pancreatic ductal adenocarcinoma reveals antigen-presenting cancer-associated fibroblasts. *Cancer Discov* 2019;9(8):1102-1123. PMID: 31197017. doi:10.1158/2159-8290.CD-19-0094
5. Biffi G, Oni TE, et al. IL1-induced JAK/STAT signaling is antagonized by TGFβ to shape CAF heterogeneity in pancreatic ductal adenocarcinoma. *Cancer Discov* 2019;9(2):282-301. PMID: 30366930. doi:10.1158/2159-8290.CD-18-0710
6. Rhim AD, Oberstein PE, et al. Stromal elements act to restrain, rather than support, pancreatic ductal adenocarcinoma. *Cancer Cell* 2014;25(6):735-747. PMID: 24856585. doi:10.1016/j.ccr.2014.04.021
7. Özdemir BC, Pentcheva-Hoang T, Carstens JL, et al. Depletion of carcinoma-associated fibroblasts and fibrosis induces immunosuppression and accelerates pancreas cancer with reduced survival. *Cancer Cell* 2014;25(6):719-734. PMID: 24856586. doi:10.1016/j.ccr.2014.04.005
8. Qu C, Wang Q, et al. Cancer-associated fibroblasts in pancreatic cancer: should they be deleted or reeducated? *Integr Cancer Ther* 2018;17(4):1016-1019. PMID: 30136592. doi:10.1177/1534735418794884
9. Gatenby RA, Silva AS, Gillies RJ, Frieden BR. Adaptive therapy. *Cancer Res* 2009;69(11):4894-4903. PMID: 19487300. doi:10.1158/0008-5472.CAN-08-3658
10. Zhang J, Cunningham JJ, Brown JS, Gatenby RA. Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. *Nat Commun* 2017;8(1):1816. PMID: 29180633. doi:10.1038/s41467-017-01968-5
11. West J, You L, et al. Towards multidrug adaptive therapy. *Cancer Res* 2020;80(7):1578-1589. PMID: 31948939. doi:10.1158/0008-5472.CAN-19-2669
12. McGehee C, Mori Y. A mathematical framework for comparison of intermittent versus continuous adaptive chemotherapy dosing in cancer. *NPJ Syst Biol Appl* 2024;10(1):140. PMID: 39614108. doi:10.1038/s41540-024-00461-2
13. Thomsen FJ, Dubbeldam JLA. Controllability in a class of cancer therapy models with co-evolving resistance. *Discrete Contin Dyn Syst B* 2026;32:40-64. doi:10.3934/dcdsb.2025121
14. Zheng J, Wu M, et al. Network pharmacology to unveil the biological basis of health-strengthening herbal medicine in cancer treatment. *Cancers (Basel)* 2018;10(11):461. PMID: 30469422. doi:10.3390/cancers10110461
15. Zhou Z, Nan Y, et al. Hawthorn with "homology of medicine and food": a review of anticancer effects and mechanisms. *Front Pharmacol* 2024;15:1384189. PMID: 38915462. doi:10.3389/fphar.2024.1384189
16. Zhao X, Xiu J, et al. Network pharmacology and bioinformatics study of six medicinal food homologous plants against colorectal cancer. *Int J Mol Sci* 2025;26(3):930. PMID: 39940699. doi:10.3390/ijms26030930
17. Chen L, Qu C, et al. Chinese herbal medicine (Qingyihuaji) suppresses invasion-promoting capacity of cancer-associated fibroblasts in pancreatic cancer. *PLoS One* 2014;9(4):e96177. PMID: 24781445. doi:10.1371/journal.pone.0096177
18. Hu X, Peng X, et al. Shikonin reverses cancer-associated fibroblast-induced gemcitabine resistance in pancreatic cancer cells by suppressing MCT4-mediated reverse Warburg effect. *Phytomedicine* 2024;123:155214. PMID: 38134861. doi:10.1016/j.phymed.2023.155214
19. Norton KA, Gong C, et al. Multiscale agent-based and hybrid modeling of the tumor immune microenvironment. *Processes (Basel)* 2019;7(1):37. PMID: 30701168. doi:10.3390/pr7010037
20. Shiau C, Cao J, et al. Spatially resolved analysis of pancreatic cancer identifies therapy-associated remodeling of the tumor microenvironment. *Nat Genet* 2024;56(11):2466-2478. PMID: 39227743. doi:10.1038/s41588-024-01890-9
21. Middleton G, Silcocks P, et al. Gemcitabine and capecitabine with or without telomerase peptide vaccine GV1001 in patients with locally advanced or metastatic pancreatic cancer (TeloVac): an open-label, randomised, phase 3 trial. *Lancet Oncol* 2014;15(8):829-840. PMID: 24954781. doi:10.1016/S1470-2045(14)70236-0
22. Jo JH, Kim YT, et al. Efficacy of GV1001 with gemcitabine/capecitabine in previously untreated patients with advanced pancreatic ductal adenocarcinoma having high serum eotaxin levels (KG4/2015): an open-label, randomised, phase 3 trial. *Br J Cancer* 2024;130(1):43-52. PMID: 37903909. doi:10.1038/s41416-023-02474-w
23. Gong Z, Ye H, et al. S-allyl-cysteine attenuates carbon tetrachloride-induced liver fibrosis in rats by targeting STAT3/SMAD3 pathway. *Am J Transl Res* 2018;10(5):1337-1346. PMID: 29887949
24. Tsukioka T, Takemura S, et al. Attenuation of bleomycin-induced pulmonary fibrosis in rats with S-allyl cysteine. *Molecules* 2017;22(4):543. PMID: 28353632. doi:10.3390/molecules22040543
25. Kawasaki H, Nussbaum G. Therapeutic potential of garlic, aged garlic extract and garlic-derived compounds on pancreatic cancer (Review). *Biomed Rep* 2025;22(3):54. PMID: 39926043. doi:10.3892/br.2025.1932
26. Hu J, Liu Y, et al. Eupatilin ameliorates hepatic fibrosis and hepatic stellate cell activation by suppressing β-catenin/PAI-1 pathway. *Int J Mol Sci* 2023;24(6):5933. PMID: 36983006. doi:10.3390/ijms24065933
27. Bae SH, Zheng YF, et al. Stereoselective determination of ginsenosides Rg3 and Rh2 epimers in rat plasma by LC-MS/MS: application to a pharmacokinetic study. *J Sep Sci* 2013;36(11):1904-1912. PMID: 23559579. doi:10.1002/jssc.201300107
28. Wu C, Li P, et al. Efficacy and safety of Shenyi capsule (ginsenoside Rg3) as adjuvant therapy for cancer: an overview of systematic reviews and meta-analyses. *Integr Cancer Ther* 2025;24:15347354251396519. PMID: 41342201. doi:10.1177/15347354251396519
29. Anand P, Kunnumakkara AB, et al. Bioavailability of curcumin: problems and promises. *Mol Pharm* 2007;4(6):807-818. PMID: 17999464. doi:10.1021/mp700113r

---

## Supplementary
- **Table S1** — full ABM parameter list with baseline values, units, perturbable flag, and literature grounding. See `S1_parameters.md` (24 parameters across global/tumor/myCAF/CD8/resistance groups + organ-context overrides).
- **Table S2** — from food-medicine entities to defined drug products: each resolved active ingredient (SAC, eupatilin, 20(S)-Rg3, curcumin) scored on nine pharmaceutical-development attributes (active ingredient, extraction, standardized content, purity, stereochemistry, formulation, oral bioavailability, tissue distribution, batch consistency), the mapped model axis, and defined/partial/undefined status. See `S2_drug_products.md`.
- **ODD protocol** — a full Overview–Design-concepts–Details description of the ABM (Grimm et al.), covering entities and state variables, the fixed per-step process schedule and update semantics, design concepts (emergence, sensing, stochasticity, observation), initialization, and every submodel with its update equations. See `S3_ODD_protocol.md`.
- **S3 Sensitivity analysis** — see `S3_sensitivity_draft.md` (resistance_cost sweep + OAT tornado; Fig. S3).
- **Supplementary Data** — systematic survey queries, counts, and retrieved records (`docs/literature_search/`).

## Master editorial checklist (pre-submission)
- [x] Ref [13]: published version confirmed & inserted (Thomsen & Dubbeldam, DCDS-B 2026;32:40-64)
- [x] Ref [G3]: resolved (Farrokhian et al., Sci Adv 2022;8:eabm7212)
- [x] Build Supplementary Table S1 (parameters + citations) — S1_parameters.md
- [x] Verify all in-text citation numbers against reference list — all 20 cited, all valid, no orphans
- [x] Author/affiliation block, corresponding author, ORCID, repository URL (§2.9), license (MIT) — complete
- [ ] Decide Abstract numeric inclusion (3.4–3.6 values) per target journal word limit  *(journal-dependent)*
- [ ] Target journal selection → convert citation style; fix figure numbering to house style  *(journal-dependent)*
- [x] Reframe to condition-dependent stromal control; compounds demoted to case study (title, abstract, §1, §3.6–3.7, §4.7)
- [x] Terminology pass: toxicity→modeled exposure; ground-truth→author-provided reference; patient-grounded→patient-geometry initialized; validated→positive-control support
- [x] Demote compound-ranking (→Fig S8) and dose-optimization (→Fig S9) to Supplement; promote molecular grounding to Fig 6; GV1001 caption softened
- [x] Expand §4.5 limitations: confinement-encoded, caf_protumor=0/IL-6, CAF plasticity, CD8-only immunity, non-PDAC resistance, no observation model, local≠survival, seeds/calibration
- [x] Delete control_score from paper; replace with progression-constrained multi-objective/Pareto evaluation (§2.7, §3.5, Fig 5 caption); reviewer #3
- [x] §3.6: remove single-seed ranking from main text (→ Fig S8, marked exploratory/superseded); lead with 30-seed progression-constrained Pareto; reviewer #4
- [x] CRT spatial analysis rewritten to patient-level statistics (Mann–Whitney U + Cliff's delta + bootstrap 95% CI + BH correction; n=9 vs 6, none significant); causal language removed; new Fig 3 (rim_scotia_stats.png); reviewer #6
- [x] Well-mixed (mean-field) vs spatial ABM comparison proving spatial necessity — geometry-blindness + strategy inversion (Fig. S18, §3.10); reviewer #2
- [x] Fig 3 glyph fixed (ASCII minus); Fig 5 shows no control_score; Fig S8 marked superseded
- [x] Recompute Figure 4 on a denser 6×6 grid, 5 seeds/cell, with probabilistic phase-boundary uncertainty P(keep-stroma) + benefit CI (phase_dense.png; §3.4); reviewer #1
- [x] Spatial-necessity ablation (Major Comment 1) answered by mapping to existing evidence in §3.10: well-mixed = no-confinement (Fig S18, inverts strategy to always-reduce), randomized geometry at matched counts (Fig S15, 12× outcome shift), immune-exclusion contribution = phase-map α-axis (Fig 4a flips keep→reduce). Central question ("does geometry change the selected strategy?") answered yes.
- [x] Recompute key rankings with 30 seeds + progression-vs-exposure Pareto frontier and rank-stability (Fig. S10; overturns several single-agent rankings)
- [x] Restore CAF pro-tumor axis (caf_protumor + new caf_survival drug-tolerance); re-map phase map with it on (Fig. S11; resource regime 7/16→6/16, persists)
- [x] Add CA19-9 observation model (interval, noise, lag, non-secretor, min-duration, safety); adaptive advantage survives but attenuates (Fig. S12)
- [x] Add intermittent (clinical-approx q28) gemcitabine arm + same-agent schedule-only comparison (Fig. S13); reframe continuous as within-model reference
- [x] Generalize GV1001 to "illustrative immune-priming agent" in main text; keep the GV1001 name + trial refs only in Fig. S4 caption
- [x] Resolve food entities to defined APIs (garlic→SAC, mugwort→eupatilin, ginseng→20(S)-Rg3) + drug-product attribute table (Table S2, refs [26-29])
- [x] Monte Carlo epistemic uncertainty over compound assumptions (effect, exposure weight, bioavailability, synergy) on API regimens (Fig. S14); natural control collapses to 0-2% vs gem 96% — not robust to pharmacology uncertainty
- [x] Non-trivial spatial predictions (§3.10): geometry changes immune-therapy efficacy ~12× at fixed abundance (Fig. S15a), CAF reduction widens invasion front more when confined (Fig. S15b), open-first beats immune-first modestly (Fig. S16), and peritumoral myCAF density (|r|=0.92) — not abundance (|r|=0.56) — predicts the immune regime (Fig. S17)
- [x] CRTL defined inline (§2.2, chemoradiation + additional line; excluded from untreated-vs-CRT comparison)
- [x] Editorial pass: cover-page placeholder removed; CRT causal→cross-sectional language (§3.3, Figs 2b/3); novelty "empty/confirmed"→"no directly matching study" and query count unified to 12 (Fig 1); Fig 4 caption "resource"→"non-zero optimum"; HCC decoupled from results; per-experiment seed table (§2.8); archived-release note (§2.9)
- [x] ODD-protocol ABM specification (Grimm et al.): entities/state variables, fixed per-step schedule + update semantics, design concepts, initialization, all submodel equations, boundary/RNG conventions (S3_ODD_protocol.md; referenced §2.4)
