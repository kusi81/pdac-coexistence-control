# Supplementary Table S1 — Agent-based model parameters

Baseline values, units, perturbability, and grounding for the ABM (Methods §2.4).
Source: `pipeline/abm.py` (`DEFAULT_PARAMS`, `PERTURBABLE`, `CONTEXTS`). Parameters
are literature-grounded or model-calibration choices, **not fit to data**; the
sensitivity analysis (§S3, Fig. S3) quantifies robustness to ±50% variation.
"Perturb." = whether a compound regimen can scale the parameter (drug action).

## Global

| Parameter | Baseline | Unit | Perturb. | Description / grounding |
|---|---|---|:--:|---|
| `field_um` | 1500 | µm | — | Square tissue field side length. |
| `dt_days` | 0.5 | day | — | Simulation time step. |
| `seed` | 0 | — | — | RNG seed (results averaged over seeds 42/7/123). |

## Tumor

| Parameter | Baseline | Unit | Perturb. | Description / grounding |
|---|---|---|:--:|---|
| `k_prolif` | 0.11 | /day | ✅ | Division probability per day; set for ~6-day doubling so immune killing can compete. Anti-proliferative compounds lower it. |
| `tumor_density_cap` | 10 | cells | — | Max tumor neighbors within kill radius (contact inhibition). |
| `k_tumor_apoptosis` | 0.01 | /day | ✅ | Baseline spontaneous apoptosis. |

## myCAF barrier

| Parameter | Baseline | Unit | Perturb. | Description / grounding |
|---|---|---|:--:|---|
| `k_caf_activate` | 0.18 | /day | ✅ | Stroma→myCAF activation probability near tumor; anti-fibrotic compounds lower it so turnover shrinks barrier mass. |
| `caf_ring_um` | 150 | µm | — | Distance band from tumor within which activation occurs. |
| `caf_cap_per_tumor` | 0.9 | — | — | myCAF saturation cap per tumor cell. |
| `caf_protumor` | 0.0 (PDAC) | — | — | Local proliferation boost from myCAF. 0 = pure barrier (PDAC); large for HCC "co-opted stroma" context [7; PMC7419619]. |
| `caf_boost_ref` | 18 | cells | — | Local myCAF count giving full pro-tumor boost (higher → more linear/less saturating). |
| `k_caf_death` | 0.07 | /day | — | myCAF deactivation/turnover; the route by which reduced activation lowers barrier mass. |
| `caf_confine` | 0.8 | — | — | Physical containment: probability a tumor daughter is blocked from a stroma-dense location (scaled by local myCAF). 0 recovers a no-containment model. |
| `caf_confine_ref` | 3 | cells | — | Local (20 µm) myCAF count for full confinement/pressure/drug-block effect; set to the observed peritumoral myCAF density. |
| `caf_pressure` | 1.2 | — | — | Mechanical exclusion: local carrying capacity scaled by exp(−`caf_pressure`·ρ/ref), limiting tumor packing/spread where stroma is dense. |
| `caf_drug_block` | 0.6 | — | — | Impaired drug delivery: drug's anti-proliferative effect attenuated by exp(−`caf_drug_block`·ρ/ref) in stroma-dense regions. |
| `caf_survival` | 0.0 (baseline) | — | — | CAF-induced drug tolerance (signaling): local myCAF further attenuates the drug's effect on sensitive cells by exp(−`caf_survival`·ρ/ref), multiplicative with `caf_drug_block` but representing IL-6/JAK–STAT-type therapy-induced resistance [18,20] rather than physical exclusion. Switched on (0.4) in the pro-tumor robustness check (Fig. S11). |

## CD8⁺ T cell

| Parameter | Baseline | Unit | Perturb. | Description / grounding |
|---|---|---|:--:|---|
| `cd8_speed_um` | 120 | µm/day | ✅ | Migration speed toward nearest tumor (unobstructed). Immunomotility compounds raise it. |
| `cd8_barrier_alpha` | 0.9 | — | — | Attenuation of CD8 migration by corridor myCAF density — the key immune-exclusion gate. |
| `corridor_um` | 30 | µm | — | Corridor width over which barrier density is evaluated. |
| `kill_radius_um` | 20 | µm | — | Contact distance at which CD8 can kill tumor. |
| `k_kill` | 1.0 | /day | ✅ | Tumor death probability per day on contact (~1 kill/day). |
| `cd8_recruit` | 30 | cells/day | ✅ | New CD8 influx at tissue margin. Immunorecruiting compounds raise it. |
| `k_cd8_death` | 0.03 | /day | — | CD8 turnover. |
| `cd8_entry_margin_um` | 60 | µm | — | Marginal band width where new CD8 enter. |

## Resistance (adaptive-therapy core)

| Parameter | Baseline | Unit | Perturb. | Description / grounding |
|---|---|---|:--:|---|
| `init_resistant_frac` | 0.01 | — | — | Pre-existing resistant fraction (small phenotypic pre-pool). PDAC gemcitabine resistance is predominantly non-genetic/adaptive (epigenetic, microenvironmental), modeled as pre-pool + phenotype switching rather than point mutation. |
| `mutation_rate` | 0.001 | /division | — | Sensitive→resistant phenotypic switching per division. |
| `resistant_immune_evasion` | 0.45 | — | — | Probability a resistant clone evades CD8 killing (PD-L1↑ / antigen–MHC-I(B2M) loss). Grounded in [G1]. |
| `resistance_cost` | 0.24 | — | — | Resistant fitness cost: resistant cells grow at ~76% of sensitive rate drug-free (NSCLC-derived) [G3]; enables competitive suppression of resistance during holidays. Whether a cost is required depends on turnover [G2]. |

## Organ contexts (overrides)

| Context | `caf_protumor` | `k_prolif` | `k_caf_activate` | Rationale |
|---|:--:|:--:|:--:|---|
| PDAC | 0.0 | 0.11 | 0.18 | Desmoplasia acts as a barrier (immune exclusion dominant); stromal ablation worsened outcomes [6,7]. |
| HCC | 0.9 | 0.11 | 0.22 | Cirrhosis (activated HSC) acts as pro-tumor soil — fibrosis promotes proliferation [PMC7419619]. |

## Perturbable parameters (drug action targets)
`k_prolif` (anti-proliferative), `k_caf_activate` (anti-fibrotic), `k_kill` &
`cd8_recruit` & `cd8_speed_um` (immunomodulatory), `k_tumor_apoptosis`
(pro-apoptotic). A regimen applies per-compound multiplicative factors to these,
scaled by dose; drug action affects sensitive tumor cells only.

## Grounding references (supplementary)
- **[G1]** Juneja VR, McGuire KA, et al. PD-L1 on tumor cells is sufficient for immune evasion in immunogenic tumors and inhibits CD8 T cell cytotoxicity. *J Exp Med* 2017;214:895-904. PMID: 28302645. doi:10.1084/jem.20160801
- **[G2]** Strobl MAR, West J, et al. Turnover modulates the need for a cost of resistance in adaptive therapy. *Cancer Res* 2021;81:1135-1147. PMID: 33172930. doi:10.1158/0008-5472.CAN-20-0806
- **[G3]** Farrokhian N, Maltas J, Dinh M, et al. Measuring competitive exclusion in non-small cell lung cancer. *Sci Adv* 2022;8:eabm7212. doi:10.1126/sciadv.abm7212 (내성 적합도 비용 ~76% 근거).
- [6,7] = main-text refs (Rhim 2014 PMID 24856585; Özdemir 2014 PMID 24856586).
- PMC7419619 (HCC stroma pro-tumor), gemcitabine non-genetic resistance review — Methods 각주로 정식화 예정.

## 편집 체크리스트
- [ ] [G3] Science Adv abm7212 저자·제목 확정(또는 대체 인용)
- [ ] PMC7419619 → 정식 PMID·서지로 변환
- [ ] gemcitabine 비유전적 내성 리뷰 인용 1편 확정(§2.4 각주)
- [ ] 독성 가중치(물질별)·synergy 항은 별도 Table S2로 분리할지 결정
