# A spatially grounded agent-based framework for coexistence-control of pancreatic cancer using food-medicine-homology compounds

*Draft manuscript — assembled from section drafts (2026-07-21). In-silico, hypothesis-generating study. Placeholder items marked [ ].*

**Author:** Seung-Il Kim (김승일 / 金昇日), Independent Researcher [+ collaborating PI to be added]
**Corresponding author:** Seung-Il Kim · kusi81kim@gmail.com · **ORCID:** [0009-0007-5965-9212](https://orcid.org/0009-0007-5965-9212)

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
agent-based modeling of the tumor microenvironment. A systematic PubMed survey
confirmed that, while each line is individually well populated, their integration
is essentially absent. Our framework builds a spatial agent-based model of the
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
stromal state—not depletion or preservation per se—the control target. Within this
framework, adaptively scheduled, low-exposure food-medicine-homology regimens
controlled tumor burden under specified conditions, though any regimen's advantage
depends on assumed compound weights and the stromal regime. We stress that these
outputs are testable hypotheses—candidate stromal targets, combinations, and
schedules—intended to focus, not replace, experimental validation. The framework offers a reproducible route from real spatial data to
prioritized, low-toxicity control strategies for PDAC coexistence.

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
sits. Into this simulation we added compounds that are eaten as food and used as
traditional medicine in Northeast Asia. The model lets us test many combinations,
timings, and doses in silico and rank the most promising, low-toxicity strategies.
These predictions are hypotheses meant to guide future laboratory experiments, not
to replace them.

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
and remedy (e.g., garlic, ginseng, Platycodon, mugwort, hawthorn)—offer
multi-target activity with favorable safety and accessibility, making them
attractive backbones for chronic control regimens rather than acute cytotoxic
bursts [14,15]. However, the computational study of these compounds in oncology has
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
explicit benefit–cost trade-off. Rather than optimizing for eradication, we define a
*coexistence-control* objective with explicit resistance dynamics and a resistance
fitness cost, and we first map, over the confinement-versus-immunosuppression
trade-off, the conditions under which preserving/modulating the stroma outperforms
depleting it; we then prioritize compound combinations, sequences, and doses under
adaptive scheduling within that landscape. We emphasize that the framework is a
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
tumors (9 untreated, 6 CRT, 1 CRTL) with a 1,009-gene panel and **author-provided**
major-type and CAF-subtype annotations (717,493 cells); global pixel coordinates
were converted to microns (0.12028 µm/px). MIBI-TOF (colorectal, squidpy) served as
a negative-architecture control. The two PDAC datasets play complementary roles: the
author-annotated CosMx cohort provides the annotation ground truth for metric
validation (§3.2), while the Xenium cohort illustrates the limits of marker-based
typing.

### 2.3 Barrier and containment metrics
All spatial metrics operate on (coordinates, cell-type labels) and use permutation
nulls that preserve tissue geometry. (i) **Barrier score**: the fraction of
straight-line tumor→immune paths that pass within a corridor (default 30 µm) of a
barrier-cell (myCAF), compared against a matched-null that resamples barrier
positions at fixed density; reported as a z-score, so the metric tests
*interposition geometry* rather than abundance. (ii) **Rim enrichment**: for a
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
cells (with iCAF/macrophage as context). Update rules per step:

- **Tumor proliferation.** Each tumor cell divides with probability `k_prolif`
  (baseline 0.11/day; doubling ~6 days), gated by local contact inhibition
  (≤`tumor_density_cap` neighbors within the kill radius); baseline apoptosis
  `k_tumor_apoptosis`.
- **myCAF barrier (activation/turnover).** Peri-tumoral stroma activates to myCAF with
  probability `k_caf_activate` within `caf_ring_um` (150 µm) of tumor, up to
  `caf_cap_per_tumor`, and reverts (turnover) at `k_caf_death`; anti-fibrotic
  perturbation lowers activation so turnover reduces barrier mass. An optional
  `caf_protumor` term boosts local proliferation (0 for PDAC; large for the HCC context,
  reflecting cirrhosis-as-pro-tumor soil).
- **myCAF physical containment.** Local myCAF density ρ (within the kill radius) exerts
  three effects on the tumor. (i) *Confinement*: a dividing cell's daughter is blocked
  from being placed into a stroma-dense location with probability
  `caf_confine`·min(ρ/`caf_confine_ref`, 1), so the tumor cannot expand through the
  stromal wall. (ii) *Mechanical pressure*: the local carrying capacity is scaled by
  exp(−`caf_pressure`·ρ/`caf_confine_ref`), so the tumor cannot pack where stroma is
  dense. (iii) *Impaired drug delivery*: the drug's anti-proliferative effect on
  sensitive cells is attenuated by exp(−`caf_drug_block`·ρ/`caf_confine_ref`). Setting
  these to zero recovers a model in which myCAF affects only immunity.
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
and two costs: immune exclusion (via `cd8_barrier_alpha`) and impaired drug delivery
(effect iii). Their balance determines whether preserving or reducing stroma aids
control, which we map in §3.4. Parameters are literature-grounded rather than fit to
data; full values and citations are in Supplementary Table S1. Two organ contexts
(PDAC, HCC) share the stellate/TGF-β axis and differ only in `caf_protumor`.

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
simulated trajectory we compute: time-to-progression (first time tumor reaches 1.5×
baseline; if never, the horizon is used and the run is *censored*), peak and final
burden, final resistant fraction, cumulative toxicity, burden AUC, and a composite
`control_score` = time-to-progression / (cumulative toxicity + 1). Because
`control_score` is artifactually inflated for zero-toxicity (untreated) runs,
strategies are ranked by progression control first (censored vs progressed) and then
by toxicity. Adaptive scheduling toggles treatment on above an upper band (`adapt_on`)
and off below a lower band (`adapt_off`) around the tumor-burden reference.

### 2.8 In-silico experimental design
Simulations were run on synthetic tissues generated with matched cell counts in two
architectures—`contained` (myCAF rings around tumor islets, CD8 excluded) and
`diffuse` (identical counts, no architecture)—so that abundance is held constant and
only geometry differs. Experiments comprised: (E1) untreated vs continuous
maximum-dose vs adaptive scheduling; (E2) single agents vs combinations vs timed
sequential cycles; (E3) a seed-averaged dose × drug-holiday grid search; and a
sensitivity analysis combining a `resistance_cost` sweep with one-at-a-time ±50%
perturbation of nine key parameters. Stochastic results are averaged over three seeds
(42, 7, 123).

### 2.9 Implementation and availability
The framework is implemented in Python (3.13) using scanpy/anndata, squidpy, pyarrow,
NumPy/SciPy, and Matplotlib, with an interactive Streamlit dashboard exposing the
spatial, simulation, optimization, and molecular-viewer modes. Spatial data: Xenium
GSE274673; CosMx from Mendeley Data (doi:10.17632/kx6b69n3cb.1) [20]. All analysis
code is available at https://github.com/kusi81/pdac-coexistence-control.

---

## 3. Results

### 3.1 The components of our approach are individually well-established, but their integration is unoccupied
To position the framework, we ran a systematic PubMed survey (ten Boolean queries
spanning the intersecting dimensions of our approach; Methods §2.1) and quantified
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
food-medicine-homology compounds to modulate the myCAF barrier—was empty in our
survey, defining the gap this work occupies. We note the limits of this analysis: it
queries PubMed only and is a targeted novelty search rather than a PRISMA systematic
review; nonetheless, the near-absence of integrated hits amid abundant single-axis
literature supports that the contribution lies in the integration.

### 3.2 The containment metric is validated by an author-annotated positive control, and localizes annotation—not metric—as the limiting factor
On synthetic tissues, the barrier score cleanly separated an architecture in which
stroma is interposed between tumor and immune cells (contained; z ≈ +22.3) from one in
which the same cell counts are randomly intermixed (diffuse; z ≈ +0.9), confirming
that the matched-null construction responds to spatial *geometry* rather than
abundance (Fig. 2a). We next tested a known biological positive control—myCAFs lie
closer to malignant cells than iCAFs [3,20]—on the author-annotated CosMx dataset
(717,493 cells across 16 tumors). myCAF was significantly closer to malignant cells
than iCAF in 15 of 16 tumors (the single exception showed a 0.2 µm difference,
p = 0.279, a tie), and the gap widened after therapy (−2.6 µm untreated vs −18.8 µm
CRT) (Fig. 2b). Given faithful annotations, our metric recovers the expected
architecture.

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

### 3.3 myCAFs form a tumor-adjacent barrier that tightens under therapy while cytotoxic T cells remain excluded
We characterized the peritumoral rim (30 µm shell) composition and its change with
neoadjuvant chemoradiotherapy in the CosMx cohort (9 untreated, 6 CRT; Fig. 3). myCAFs
were the dominant tumor-adjacent population, strongly enriched at the rim in both
groups and further enriched after therapy (mean rim z = +9.1 → +15.7), whereas iCAFs
were excluded and displaced farther out (z = −3.5 → −9.2). Pericytes shifted from
neutral to tumor-adjacent after therapy (z = −1.6 → +3.3), consistent with vascular
remodeling. The immune compartment showed persistent exclusion: cytotoxic CD8⁺ T cells
remained depleted from the rim (z = −5.6 → −4.0), with CD4⁺ T and regulatory T cells
similarly excluded and plasma cells most strongly excluded (z = −11.3 → −4.1);
macrophages were the only rim-enriched immune lineage, declining after therapy
(z = +6.8 → +2.4).

These observations refine the intuitive expectation that therapy relieves immune
exclusion. In this cohort, CRT did *not* recruit cytotoxic T cells to the tumor rim:
CD8⁺ T cells remained excluded while the myCAF barrier grew tighter. Rather than
contradicting the containment hypothesis, this supports it—the myCAF barrier persists
and is reinforced under treatment, maintaining spatial separation between malignant
cells and cytotoxic effectors, precisely the configuration our control framework seeks
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

The optimal stromal level was regime-dependent. Where confinement was strong and immune
exclusion weak, the optimum was interior (a non-zero stromal level), and preserving that
level reduced tumor burden by up to ~0.26-fold relative to full stromal depletion—the
myCAF-as-controllable-resource regime. Where immune exclusion dominated, the optimum
collapsed to zero (full stromal reduction was best), because the immune- and
drug-shielding costs of stroma outweighed its confinement benefit. Across the sampled
grid, the resource regime occupied 7 of 16 conditions (Fig. 4a), with the benefit of
keeping stroma concentrated in the low-immune-exclusion, high-confinement corner
(Fig. 4b).

This reframes the central question from "deplete or preserve the stroma?" to "what is the
optimal stromal state, and when?" It also formalizes the CAF-biology debate: the same
peritumoral myCAF barrier can be net-protective or net-harmful depending on the local
balance, so both the failure of stromal ablation [6,7] and reports that CAF reduction
improves immunotherapy emerge as different regimes of one trade-off. We stress that the
phase boundary depends on parameters that are mechanistically motivated but not
data-fitted; where real PDAC lies on it is an empirical question this framework poses
rather than answers (Discussion §4.5).

### 3.5 Adaptive scheduling achieves tumor coexistence at a fraction of the toxicity of continuous dosing
On a controllable-but-not-eradicable tumor, we compared no treatment, continuous
maximum-dose therapy, and adaptive on/off dosing (Fig. 5). Continuous dosing drove the
sensitive population to extinction (final burden 0.00× baseline) but at the highest
cumulative toxicity (120 units) and left a purely resistant residue. Adaptive dosing
held the tumor at a stable 0.82× baseline—coexistence rather than eradication—at
roughly one-fifth the toxicity (25 units), keeping the resistant fraction low (0.03);
both treated arms remained below the progression threshold over 150 days, whereas the
untreated tumor progressed at day 103. Because our `control_score` is artifactually
high for the zero-toxicity untreated arm, we rank by progression control first and then
by toxicity, under which adaptive scheduling is preferred (control_score 5.7 vs 1.2).

### 3.6 Food-medicine-homology regimens control tumor burden at low predicted toxicity, with anti-CAF agents minimizing resistance
We evaluated single agents and combinations under adaptive scheduling, benchmarked
against continuous gemcitabine (Fig. 6; in-silico predictions). Every natural regimen
maintained control at markedly lower predicted toxicity than gemcitabine (5–38 vs 128
units). The most favorable control-per-toxicity profiles were single anti-CAF/
immunomodulatory agents—garlic (toxicity 5, 1.0× burden, resistant 0.05, control_score
25.0) and wild ginseng (toxicity 5, 0.9×, 24.0)—stabilizing the tumor near baseline at
~1/25 the toxicity of gemcitabine. Two patterns stand out. First, deeper suppression
did not equate to better control: an aggressive combination (garlic + mugwort) reduced
burden to 0.4× but at a sharply elevated resistant fraction (0.34), an in-silico
illustration of competitive release. The balanced combination curcumin + garlic +
ginsenoside-Rg3 achieved moderate suppression (0.8×) with the lowest resistance (0.01)
at modest toxicity (11), and was carried forward for dose optimization (§3.7). Second,
anti-fibrotic activity alone was insufficient: a purely anti-fibrotic pairing (danshen
+ astragaloside) barely held the tumor (1.2×) at the highest natural-arm toxicity (38).
These results support a control-oriented design principle—pair a low-toxicity anti-CAF
backbone with a resistance-sparing anti-proliferative partner rather than maximizing
cytotoxic pressure.

### 3.7 Dose and drug-holiday optimization further lowers predicted toxicity
For the lead combination (curcumin + garlic + ginsenoside-Rg3), a seed-averaged grid
search over dose intensity and off-threshold under adaptive scheduling (200 days;
Fig. 7) found an optimum at 80% dose with an off-threshold of 0.5, achieving control at
cumulative toxicity 27 (resistant 0.02, final 0.8× baseline). This modestly but
consistently outperformed full-dose adaptive scheduling (toxicity 32 at comparable
control), indicating that a sub-maximal, adaptively timed dose preserves control while
further trimming toxicity—reinforcing that, under a control objective, less drug
delivered adaptively can be better.

### 3.8 Multiscale molecular grounding links compounds to model parameters under a transparent evidence hierarchy
Each compound is annotated with its principal target(s) and an explicit evidence tier—
experimental co-crystal, docking prediction, or pathway-level mechanistic inference
(Fig. 8). Only the two conventional agents are supported by experimental co-crystals in
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
mechanistic-tier compounds that dominate our low-toxicity regimens are grounded in
pathway pharmacology, and their molecular engagement remains a hypothesis for
experimental confirmation.

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
metric reproduced the myCAF-proximal/iCAF-distal architecture in 15 of 16 tumors and
showed the myCAF barrier tightening under chemoradiation while cytotoxic T cells
remained excluded (§3.3). We are careful in interpreting this: a denser peritumoral
myCAF rim with persistent T-cell exclusion is equally consistent with an
*immunosuppressive* barrier as with a tumor-*containing* one—the observation is a
necessary, not sufficient, condition for containment. This ambiguity is precisely
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

### 4.5 Limitations
Our study has important limitations that bound its claims. First and foremost, it is
entirely in silico: outputs are ranked, testable hypotheses, not evidence of efficacy,
and require experimental validation. Second, the ABM is phenomenological with
literature-inferred parameters rather than data-fitted; although the adaptive advantage
was robust to ±50% parameter variation (§S3), outcomes were most sensitive to tumor
proliferation and immune-barrier parameters, which therefore require empirical
calibration, and control was not guaranteed under the most aggressive proliferation
regimes. Relatedly, the confinement-versus-immunosuppression trade-off that decides
whether stroma should be preserved or reduced is governed by parameters (stromal
confinement strength, drug shielding, and immune exclusion) that are mechanistically
motivated but not data-fitted; the phase boundary is therefore qualitative, and where
real PDAC falls on it is unresolved. Third, our spatial analyses are of 2D sections; true tissue architecture is
three-dimensional, and the myCAF-proximity positive control failed on our targeted-
panel Xenium annotation, succeeding only with full-panel author annotations (§3.2)—so
real-data spatial claims rest on the CosMx cohort and remain associations in a modest
sample (n = 15). Fourth, simulations used a single synthetic architecture and an
arbitrary-unit toxicity scale that is not clinically calibrated. Fifth, even
low-toxicity natural combinations can select for resistance under aggressive
suppression (competitive release; §3.6). Finally, the novelty survey is PubMed-
restricted, and the molecular grounding for most food-medicine-homology compounds rests
on pathway-level mechanistic inference rather than solved compound-bound structures
(§3.8).

### 4.6 Future directions
Several steps follow directly. Experimentally, the ranked regimens define a low-cost
validation ladder—2D PDAC cell-line dose–response and combination/synergy assays, then
patient-derived organoid screening—prioritizing the anti-CAF-backbone combinations and
the competitive-release predictions. Computationally, priorities are: parameter
calibration against experimental rates (especially proliferation and immune-barrier
terms); global sensitivity analysis (e.g., Sobol) to capture interactions beyond the
one-at-a-time sweep; a 3D registration pipeline over serial sections to replace 2D
metrics; annotation transfer from full-transcriptome references to make CAF subtyping
robust; and molecular-dynamics or larger-scale virtual screening to strengthen the
molecular tier. Broadening the organ context (e.g., the HCC "co-opted stroma" regime)
and extending the survey to Scopus/Web of Science under a PRISMA workflow would further
generalize and harden the framework.

### 4.7 Conclusion
Within a control—rather than eradication—objective, this work provides a reproducible
computational path from real spatial data to condition-dependent stromal targets and
prioritized, low-exposure food-medicine-homology regimens for PDAC coexistence. By
implementing the myCAF barrier as a physical structure with an explicit
confinement-versus-immunosuppression trade-off, it reframes the therapeutic question
from "deplete or preserve the stroma?" to "what is the optimal stromal state, and when?",
and it yields concrete, testable hypotheses—stromal targets, compound combinations, and
schedules—to guide the experimental work that must follow.

---

## Figures

**Figure 1. Novelty landscape.** (a) Per-dimension PubMed hit counts from the
systematic survey; individual axes (blue) are well populated while integrated,
three-axis queries (red) return only off-target hits. (b) Three-camp schematic: static
network pharmacology, wet-lab herbal–CAF, and spatial/ABM modeling are each populated
but their intersection—this work—is empty. *(assets/fig1_novelty.png)*

**Figure 2. The containment metric is sound; annotation is the limiting factor.**
(a) Synthetic contained (myCAF ring; barrier z = +22.3) vs diffuse (scattered; +0.9).
(b) On author-annotated CosMx data the myCAF-proximal positive control passes in 15/16
tumors (points left of 0 = correct). (c) The same metric fails on marker-annotated
Xenium data (1/5), localizing the failure to annotation. *(assets/fig2_validation.png)*

**Figure 3. Peritumoral composition, untreated vs CRT (CosMx author annotations).**
Rim (30 µm) enrichment z by cell type; myCAF is tumor-adjacent and tightens with CRT
while CD8⁺ T cells remain excluded. *(assets/rim_scotia.png)*

**Figure 4. The myCAF barrier is a controllable resource only under specified
conditions.** Under fixed sub-maximal treatment the myCAF level was swept and the
tumor-minimizing (optimal) level located, over a grid of confinement strength
(caf_pressure) × immune-exclusion strength (cd8_barrier_alpha). (a) Optimal myCAF level
per regime: >0 = preserve/modulate stroma (resource); 0 = reduce stroma. (b) Benefit of
keeping stroma = tumor(no stroma) − tumor(optimal stroma); >0 = stroma aids control. The
resource regime occupies the low-immune-exclusion, high-confinement region.
*(assets/phase_map.png)*

**Figure 5. Control strategies.** Tumor burden, resistant fraction, and cumulative
toxicity for untreated, continuous, and adaptive dosing; adaptive achieves coexistence
at ~1/5 the toxicity of continuous. *(assets/control_strategies.png)*

**Figure 6. Food-medicine-homology regimens** ranked by control, toxicity, and
resistance vs continuous gemcitabine. *(assets/natural_adaptive_optim.png)*

**Figure 7. Dose × drug-holiday optimization** for the lead combination; optimum at 80%
dose / 0.5 off-threshold. *(assets/dose_band_optimization.png)*

**Figure 8. Multiscale molecular grounding** with evidence tiers (experimental vs
docking vs mechanistic). *(assets/drug_structures_3d.png)*

**Figure S3. Sensitivity analysis.** (a) resistance_cost sweep—adaptive low-toxicity
advantage is robust across the range; resistance_cost governs resistant-clone
suppression. (b) One-at-a-time ±50% tornado—outcome sensitivity concentrates in
tumor-immune balance parameters. *(assets/sensitivity.png)*

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

---

## Supplementary
- **Table S1** — full ABM parameter list with baseline values, units, perturbable flag, and literature grounding. See `S1_parameters.md` (24 parameters across global/tumor/myCAF/CD8/resistance groups + organ-context overrides).
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
- [ ] Define "arbitrary units" toxicity; cross-ref Methods §2.5/2.7
- [ ] Fig 3 / Fig 5 / Fig 6 minor glyph pass (unicode minus; replot from saved CSVs — cosmetic)
- [ ] CRTL definition footnote (§2.2) — confirm meaning from SCOTIA metadata (excluded from analysis)
