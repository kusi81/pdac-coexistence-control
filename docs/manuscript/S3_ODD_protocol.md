# Supplementary — ODD protocol for the PDAC tumor–stroma–immune agent-based model

This document describes the agent-based model (ABM) following the **ODD (Overview, Design
concepts, Details)** protocol of Grimm et al. (2006; update 2010, 2020), so the model can
be understood and re-implemented independently of the source code. All quantities refer to
the reference implementation in `pipeline/abm.py`; parameter symbols and baseline values
are in Supplementary Table S1, and the food-medicine→API/drug-product mapping is in
Supplementary Table S2.

---

## 1. Overview

### 1.1 Purpose and patterns
The model is a **hypothesis-generating** in-silico sandbox for studying *condition-dependent
stromal control* of pancreatic ductal adenocarcinoma (PDAC): whether, and under what
spatial conditions, the myofibroblastic cancer-associated fibroblast (myCAF) barrier should
be preserved/modulated versus depleted to restrain—not eradicate—the tumor, and how
low-exposure, adaptively scheduled compound regimens perform under an explicit
resistance-and-coexistence objective. It is **not** a validated predictor of patient
response; parameters are grounded to literature direction and order of magnitude, not
fitted to specific patients or drugs. The model is considered to reproduce qualitative
*patterns* rather than quantitative outcomes: (i) a peritumoral myCAF barrier that excludes
cytotoxic T cells; (ii) competitive release of resistance under continuous maximum-intensity
dosing that adaptive scheduling mitigates; (iii) a regime-dependent optimal stromal level.

### 1.2 Entities, state variables, and scales
**Entities.** Individual cells (off-lattice point agents) of five types: `Tumor`, `myCAF`,
`iCAF`, `CD8_T`, `Macrophage`. There are no super-individuals; the "barrier" and "tumor
mass" are emergent aggregates, not explicit agents.

**Agent state variables.**
- Position `(x, y)` in micrometres (continuous, off-lattice).
- Cell-type label ∈ {Tumor, myCAF, iCAF, CD8_T, Macrophage}.
- Resistance flag (Boolean; tumor cells only; non-tumor cells are always `False`).

**Global (model) state variables.**
- Simulation time `t` (days); baseline count `n0` (initial tumor number).
- Cumulative modeled exposure `cum_toxicity`.
- Treatment state `drug_on` (Boolean) and current phase label.
- Observation-model state: latent biomarker `biomarker` (ratio to baseline), last measured
  value `obs_ratio`, `last_obs_t`, `last_switch_t`, consecutive-low-read counter `low_reads`.
- `crit_time` (first time tumor reaches the progression threshold).

**Scales.** Space: a square field of side `field_um` = 1500 µm (baseline). Time: discrete
steps of `dt_days` = 0.5 day; a run of *D* days is `round(D / dt_days)` steps. One tumor
"interaction radius" `kill_radius_um` = 20 µm defines local neighbourhoods; the myCAF
"corridor" width for immune gating is `corridor_um` = 30 µm; the peritumoral activation
band is `caf_ring_um` = 150 µm.

### 1.3 Process overview and scheduling
Each time step executes the following processes **in fixed order** (method `step()`):

1. **Observation update** — advance the latent biomarker toward the true burden ratio
   (first-order lag; §3.3.10).
2. **Regimen application** — decide `drug_on` (fixed, scheduled, or adaptive/observed),
   set the effective (drug-modified) parameter vector `eff`, and accrue exposure
   `cum_toxicity += Σ toxicity(active compounds) · dt` while dosing.
3. **Tumor proliferation** — density-capped division with myCAF containment, mechanical
   pressure, drug-penetration block, signaling drug-tolerance, and optional pro-tumor boost
   (§3.3.1).
4. **myCAF activation** — a fraction of peritumoral iCAF/Macrophage cells convert to myCAF,
   up to a cap (§3.3.4).
5. **myCAF turnover** — a fraction of myCAF revert to iCAF (§3.3.5).
6. **CD8 migration and killing** — CD8 step toward the nearest tumor with speed attenuated
   by corridor myCAF density, then kill tumor within `kill_radius_um` (resistant clones
   evade) (§3.3.6–3.3.7).
7. **Tumor apoptosis** — baseline plus drug-induced apoptosis (drug-induced on sensitive
   cells only) (§3.3.8).
8. **CD8 turnover** — random CD8 removal (§3.3.9).
9. **Apply structural changes** — remove killed/apoptotic tumor and dead CD8; append
   division daughters; recruit new CD8 at the tissue margin (Poisson).
10. **Boundary clip** — clip all positions into `[0, field_um]`.
11. **Record** — append the step's summary to the history and, at requested progress
    fractions, store a full `(coords, labels)` snapshot.

**Update semantics.** Within a step, proliferation decisions use the tumor k-d tree built
at the *start* of the step, so daughters created this step do not influence divisions in the
same step (effectively synchronous for division). Cell removals and additions are collected
and applied together at step end. CD8 positions are updated in place before killing is
evaluated, so killing uses post-move positions. Processes are otherwise sequential in the
order above.

---

## 2. Design concepts

- **Basic principles.** (i) The myCAF barrier is a double-edged structure: it physically
  confines tumor expansion but simultaneously excludes CD8 T cells and impairs drug
  delivery—so its net value is conditional. (ii) Adaptive therapy (Gatenby): a fitness cost
  on resistance lets a maintained sensitive population competitively suppress resistant
  clones during treatment holidays. (iii) Compounds act as multiplicative perturbations of
  mechanistic parameters (network-pharmacology mapping).
- **Emergence.** The regime-dependent optimal stromal level, the ~order-of-magnitude
  dependence of immunotherapy efficacy on myCAF *geometry* at fixed abundance, the
  invasion-front response to CAF reduction, and competitive release all *emerge* from local
  rules; none is directly imposed.
- **Adaptation & objectives.** Cells do not optimize. The *controller* (treatment) adapts:
  it toggles dosing to keep the observed tumor burden within a band, pursuing the objective
  of long control at low cumulative exposure (evaluated post hoc by `control_metrics`).
- **Prediction.** The adaptive controller acts on either the true burden (idealized) or a
  lagged, noisy, interval-sampled CA19-9-type biomarker (observation model).
- **Sensing.** Cells sense only *local* information via k-d-tree radius queries: tumor cells
  sense neighbour counts within `kill_radius_um` (and myCAF counts for containment); CD8
  cells sense the direction/distance to the nearest tumor and the myCAF density along that
  path. The controller senses aggregate burden (or its biomarker proxy).
- **Interaction.** Local and density-mediated: contact inhibition (division cap), myCAF
  corridor attenuation of CD8 motility, CD8 killing within radius, tumor-proximal CAF
  activation.
- **Stochasticity.** All events are Bernoulli/Poisson draws from a single seeded generator
  (`numpy.random.default_rng(seed)`): division, CAF activation, turnover, killing,
  apoptosis, CD8 recruitment (Poisson), daughter placement angle/distance, resistance
  phenotype switching, and biomarker measurement noise. Replication uses independent seeds
  (Table in §2.8 of the main text).
- **Collectives.** The peritumoral myCAF "barrier" and the tumor "mass" are emergent
  collectives with no explicit representation.
- **Observation.** Per step, the model records tumor/resistant/myCAF/CD8 counts, CD8
  infiltration (fraction within 25 µm of tumor), exposure, dosing state, and the biomarker;
  full spatial snapshots are stored at requested progress points and read by the spatial
  metrics (barrier score, rim enrichment). `control_metrics` derives time-to-progression,
  peak/final burden, resistant fraction, cumulative exposure, and AUC.

---

## 3. Details

### 3.1 Initialization
Tissue is initialized in one of two ways:
- **Synthetic** (`synthetic.make_tissue`): a `contained` architecture (tumor islets ringed
  by myCAF, iCAF distal, CD8 excluded outside the ring, macrophages in stroma) or a
  `diffuse` architecture with **identical cell counts** but scattered myCAF/immune cells.
  Both are seeded, so architecture varies with seed.
- **Patient-grounded** (`patient_grounded.py`): each cell's measured position and mapped
  type from a real CosMx window are used directly as the initial condition.

At initialization, each tumor cell is independently marked resistant with probability
`init_resistant_frac`. All state counters (`t=0`, `cum_toxicity=0`, `biomarker=1`,
`drug_on=True`, `low_reads=0`) are set, and the initial state is recorded.

### 3.2 Input data
The model uses no time-varying external environmental drivers. The only exogenous inputs
are (a) the spatial initial condition (synthetic or patient tissue) and (b) the treatment
program (fixed dose, sequential `schedule`, or adaptive/observed controller). Compound
definitions (parameter multipliers, toxicity weights) are static inputs (Table S2,
`SUBSTANCES`/`TOXICITY`).

### 3.3 Submodels
Notation: `dt` = `dt_days`; `ρ` = local myCAF count within `kill_radius_um` of a focal
tumor cell; `ρ_ref` = `caf_confine_ref`; `eff[·]` = drug-modified parameter (baseline value
times the product of active-compound multipliers, dose-interpolated); `rng` = the seeded
generator.

**3.3.1 Tumor proliferation.** For each existing tumor cell, let `local` be the tumor count
within `kill_radius_um` and let the effective carrying capacity be
`eff_cap = tumor_density_cap · (1 + soil) · exp(−caf_pressure · ρ / ρ_ref)`, where
`soil = caf_protumor · min(ρ₂/caf_boost_ref, 1)` uses the myCAF count `ρ₂` within
`2·kill_radius_um` (pro-tumor boost; baseline `caf_protumor = 0`). Division is attempted
only if `local < eff_cap`. The per-step division probability is `kp · boost · dt` with
`boost = 1 + soil` and:
- resistant cell: `kp = k_prolif · (1 − resistance_cost)` (drug-independent);
- sensitive cell: `kp = k_prolif − drug_reduction · pen · surv`, where
  `drug_reduction = max(0, k_prolif − eff[k_prolif])` is the drug's anti-proliferative
  effect, `pen = exp(−caf_drug_block · ρ / ρ_ref)` is physical drug-penetration attenuation,
  and `surv = exp(−caf_survival · ρ / ρ_ref)` is signaling-based drug tolerance
  (IL-6/JAK–STAT-type; baseline `caf_survival = 0`).
On division, a daughter is placed at a random angle and distance `Uniform(4, kill_radius_um)`
from the parent. **Containment:** if the daughter's target site has myCAF density `ρ_d`, the
placement is blocked (division aborted) with probability
`caf_confine · min(ρ_d / ρ_ref, 1)`. A daughter inherits resistance from a resistant parent,
or switches to resistant with probability `mutation_rate` (phenotypic switching).

**3.3.2 Runtime population cap (analysis only).** An optional `max_tumor` cap halts new
divisions once the tumor count reaches it. This is a tractability bound used in some
compute-heavy analyses; it only affects already-progressing regimes and does not change
their progression classification (default off).

**3.3.3 Drug effect composition.** A regimen is a list of `(compound, dose)`; each compound
contributes multipliers on perturbable parameters (`k_prolif`, `k_caf_activate`, `k_kill`,
`cd8_recruit`, `cd8_speed_um`, `k_tumor_apoptosis`), linearly interpolated by dose
(`1 + (m−1)·dose`) and multiplied across compounds; an optional `synergy` term further
lowers `k_prolif`/`k_caf_activate` for ≥2-compound regimens. The product gives `eff[·]`.

**3.3.4 myCAF activation.** If the current myCAF count is below `caf_cap_per_tumor · (#tumor)`,
stromal cells (iCAF or Macrophage) whose distance to the nearest tumor is in
`(40, caf_ring_um)` µm convert to myCAF each with probability `eff[k_caf_activate] · dt`,
up to the cap. This builds/thickens the peritumoral ring; anti-fibrotic compounds lower
`eff[k_caf_activate]`.

**3.3.5 myCAF turnover.** Each myCAF reverts to iCAF with probability `k_caf_death · dt`.
Activation and turnover set an equilibrium barrier mass; lowering activation (anti-fibrotic)
reduces the equilibrium myCAF, i.e. thins the barrier over time.

**3.3.6 CD8 migration.** For each CD8 cell, find the nearest tumor. Compute corridor myCAF
density `corr` by sampling points along the CD8→tumor segment (at fractional positions
0.15–0.85, spacing `corridor_um/2`) and averaging the myCAF count within `corridor_um` of
each sample. Motility `mob = exp(−cd8_barrier_alpha · corr) ∈ (0,1]` attenuates the step,
which is `min(eff[cd8_speed_um]·dt·mob, distance-to-target)` toward the tumor. A dense
corridor thus corrals CD8 in the periphery.

**3.3.7 CD8 killing.** After moving, any CD8 within `kill_radius_um` of its nearest tumor
kills it with probability `k_kill_eff · dt`, where `k_kill_eff = eff[k_kill]·(1 −
resistant_immune_evasion)` for a resistant target and `eff[k_kill]` for a sensitive target.

**3.3.8 Tumor apoptosis.** Each tumor cell dies with probability
`apop · dt`, where `apop = k_tumor_apoptosis` for resistant cells and
`eff[k_tumor_apoptosis]` for sensitive cells (drug-induced apoptosis applies to sensitive
cells only).

**3.3.9 CD8 turnover and recruitment.** Each CD8 is removed with probability `k_cd8_death·dt`.
New CD8 are recruited as `Poisson(eff[cd8_recruit]·dt)` cells placed uniformly within
`cd8_entry_margin_um` of the field edge.

**3.3.10 Observation model (optional).** When `obs_model` is on, the latent biomarker
follows the true burden ratio `r = n_tumor/n0` with first-order lag:
`biomarker += (r − biomarker)·(dt / max(obs_lag_days, dt))`. Every `obs_interval` days the
controller reads `obs_ratio = biomarker · η`, with multiplicative log-normal noise `η`
(median 1, CV `obs_noise_cv`). A non-secretor patient (`nonsecretor=True`) yields an
uninformative marker → the controller defaults to continuous dosing. Otherwise, the on/off
decision applies: force ON if `obs_ratio ≥ obs_safety_mult`; else de-escalate to OFF only
after `obs_confirm` consecutive reads ≤ `adapt_off` **and** a minimum on-duration
`min_on_days`; re-escalate to ON when `obs_ratio ≥ adapt_on` **and** a minimum off-duration
`min_off_days` has elapsed. With `obs_model` off, the adaptive controller uses the true
burden directly (idealized).

**3.3.11 Adaptive control band.** In adaptive mode, dosing turns ON when burden ≥
`adapt_on · n0` and OFF when burden ≤ `adapt_off · n0`, with hysteresis between the bands.

### 3.4 Boundary and numerical conventions
- **Boundary:** positions are clipped into `[0, field_um]` each step (absorbing at the
  wall; no periodicity, no explicit collision resolution—off-lattice overlaps are permitted
  but limited by the local density cap).
- **Neighbourhoods:** all "within radius" queries use `scipy.spatial.cKDTree`
  ball/nearest queries; the tumor tree is rebuilt once per step (pre-division), the myCAF
  tree twice (pre-division and pre-migration).
- **Randomness:** a single `numpy.random.default_rng(seed)` drives all draws; a fixed seed
  fully determines a run (deterministic given identical inputs and library versions).
- **Time step:** `dt_days = 0.5`. Outcomes were checked to be qualitatively stable to the
  time step and field size in exploratory runs; formal step/size-convergence and
  identifiability analyses remain future work.

### 3.5 Reproducibility
The model is implemented in Python 3.13 (NumPy/SciPy, Matplotlib; scanpy/anndata/squidpy
for spatial I/O). Every figure is regenerated by a dedicated script in `pipeline/` from the
seeds listed in main-text §2.8; a version-tagged, DOI-archived release will accompany
publication so that results are reproducible from the exact commit. Code:
https://github.com/kusi81/pdac-coexistence-control.
