# Supplementary Table S2 — From food-medicine entities to defined drug products

Each food-level entity used in the regimens is resolved, where a standardized active
ingredient exists, onto a defined molecular species. This table makes explicit the
pharmaceutical-development attributes on which a *food* and a *drug product* differ—axes
the ABM's single parameter perturbation does not itself capture—and marks where each is
**defined**, **partial**, or **undefined**. This is a scope statement, not a claim that
any entry is a validated PDAC therapeutic: the modeled effect presumes an exposure and a
target engagement that remain to be established (§3.8, §4.6).

| Attribute | S-allylcysteine (garlic) | Eupatilin (mugwort, *A. argyi*) | 20(S)-ginsenoside Rg3 (ginseng) | Curcumin (turmeric) |
|---|---|---|---|---|
| **Active ingredient** | S-allyl-L-cysteine (defined) | Eupatilin, 5,7-dihydroxy-3′,4′,6-trimethoxyflavone (defined) | 20(S)-ginsenoside Rg3 (defined) | Curcumin / diferuloylmethane (defined) |
| **Source / extraction** | Aged garlic extract; forms during aqueous aging | *Artemisia argyi* leaf extract | Steam-processed (red) ginseng; heat converts Rb1→Rg3 | Turmeric rhizome; solvent extraction of curcuminoids |
| **Standardized content** | AGE standardized to SAC content | Standardized active ingredient of a marketed gastroprotective product (partial) | Standardized capsule product (Shenyi/Rg3) [28] | 95%-curcuminoid extracts are commonplace (partial) |
| **Purity attainable** | ≥98% achievable | High (pharma-grade precedent) | Epimer-resolved single stereoisomer achievable | 95% curcuminoids; residual demethoxy-curcuminoids |
| **Stereochemistry** | L-cysteine backbone; no critical epimer | No critical stereocentre | **20(S) vs 20(R) epimers differ in activity and oral PK** [27] | Keto–enol tautomerism; no chiral centre |
| **Formulation** | Oral, chemically stable, odorless [25] | Oral tablet precedent | Poorly soluble → capsule/nanocarrier dependent [26 review context] | **Poor solubility → phospholipid/nano/piperine formulations** [29] |
| **Oral bioavailability** | High (~98% in rodents) [25] | Moderate; not fully characterized (partial) | Low; formulation-dependent [27] | **Very low; rapid metabolism/glucuronidation** [29] |
| **Tissue distribution** | Crosses membranes; systemic (partial) | Not established for PDAC tissue (undefined) | Limited; carrier-enhanced (partial) | Low tissue levels without formulation (partial) [29] |
| **Batch consistency** | Assayable by SAC content | Pharma-standardized precedent | Rg3-content assay | Curcuminoid assay; source-variable |
| **Mapped model axis** | Anti-fibrotic → CAF activation (`k_caf_activate`↓) | Anti-fibrotic → CAF activation (`k_caf_activate`↓) | Immune/anti-proliferative (`k_kill`↑, `k_prolif`↓) | Anti-proliferative (`k_prolif`↓) |
| **Direct PDAC-myCAF evidence** | Untested hypothesis (liver/lung data) | Untested hypothesis (liver HSC data) [26] | Indirect (PD-L1, PI3K/Akt) | Indirect |
| **Overall status** | Defined molecule; PK strong, PDAC-myCAF action hypothetical | Partially defined; PDAC action hypothetical | Defined stereoisomer; PK/formulation-limited | Defined molecule; bioavailability-limited |

**How to read this table.** A green-light "defined" on active ingredient does not imply a
green light on the whole row. The clearest cases (SAC, eupatilin) are *anti-fibrotic*
molecules mapped to the CAF-modulating axis that the phase map (§3.4) shows is beneficial
only in the immune-exclusion-dominant regime; the immune/anti-proliferative APIs
(20(S)-Rg3, curcumin) are chemically defined but bioavailability- or formulation-limited,
so their encoded effects presume an exposure a conventional preparation may not reach. The
right development question the framework raises is therefore not "which food," but "which
molecule, at which achievable tissue exposure, engaging which axis"—the subject of the
pharmacokinetic/CMC work in §4.6.

*Evidence keys refer to the main-text reference list: [25] SAC/AGE review; [26] eupatilin
anti-fibrotic (β-catenin/PAI-1); [27] 20(S)/20(R)-Rg3 stereoselective PK; [28] Shenyi
(Rg3) standardized clinical product; [29] curcumin bioavailability.*
