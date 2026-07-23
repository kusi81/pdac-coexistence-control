# Cover Letter

**To:** The Editors, *PLOS Computational Biology*
**Re:** Submission of an original Research Article
**Date:** 23 July 2026

**Manuscript title:** *Condition-dependent stromal control of pancreatic cancer in a spatially grounded agent-based model* (subtitle: *food-medicine-homology compounds as a hypothesis-generating application*)

**Corresponding author:** Seung-Il Kim (Independent Researcher) · kusi81kim@gmail.com · ORCID 0009-0007-5965-9212

---

Dear Editors,

I am pleased to submit the above manuscript for consideration as a Research Article in *PLOS Computational Biology*. The work develops and analyzes a spatial agent-based model (ABM) of the pancreatic ductal adenocarcinoma (PDAC) tumor–stroma–immune ecosystem. Its central contribution is a *principle of condition-dependent stromal control*: by implementing the myofibroblastic CAF (myCAF) barrier as a physical structure with an explicit trade-off between confining the tumor and excluding immune cells / impairing drug delivery, the model shows that the optimal target is the stromal *state*—regime-dependent—rather than depletion or preservation per se. To my knowledge it is also the first study to connect three research lines that are individually well developed but, as a systematic PubMed survey in the paper confirms, essentially never integrated: network pharmacology of food-medicine-homology compounds, wet-laboratory studies of herbal modulation of PDAC CAFs, and agent-based modeling of the tumor microenvironment. Food-medicine-homology compounds enter as a *hypothesis-generating case study* that demonstrates how the framework prioritizes candidates—not as a claim of therapeutic efficacy.

**The central question and contribution.** Attempts to ablate the dense myofibroblastic CAF (myCAF) stroma of PDAC have paradoxically accelerated disease, implying the barrier also restrains the tumor. Rather than treating this as a paradox to be resolved by more depletion, we ask whether the myCAF barrier can be managed as a *controllable resource*. The model's key methodological feature is that it represents the myCAF barrier **physically**: local stroma confines tumor expansion but simultaneously excludes cytotoxic T cells and impairs drug delivery. This built-in trade-off means preserving stroma is not universally beneficial, and a variance-based phase analysis shows that the optimal target is the *stromal state itself*—regime-dependent—rather than depletion or preservation per se. Under an explicit coexistence-control objective with resistance dynamics, adaptively scheduled low-exposure regimens control tumor burden under specified conditions.

**Why this fits *PLOS Computational Biology*.** The paper is method-driven, fully reproducible, and explicitly hypothesis-generating rather than clinical: every output is framed as a candidate target, combination, or schedule to focus—not replace—experimental validation. It couples real single-cell spatial transcriptomics (Xenium and CosMx) to a mechanistic multiscale model, and it foregrounds methodological rigor and honest limitation-mapping over positive-result framing. Specifically, the model is grounded on real patient tissue as initial conditions, compares agents and schedules in a fair factorial design, resolves food-level entities toward standardized active ingredients (e.g., garlic → S-allylcysteine), and reports both one-at-a-time and global (Sobol) sensitivity analyses that identify tumor proliferation and the immune-exclusion barrier as the dominant controls. Where results did not support the initial hypothesis, we say so and re-frame accordingly.

**Openness.** All analyzed datasets are from public repositories (GEO GSE274673; Mendeley Data doi:10.17632/kx6b69n3cb.1), and all analysis and modeling code is openly available at https://github.com/kusi81/pdac-coexistence-control. A preprint of this work is being posted to bioRxiv concurrently with submission.

**Declarations.** This manuscript is original, has not been published elsewhere, and is not under consideration by another journal. The author received no specific funding for this work and declares no competing interests. As an independent, unaffiliated researcher, I am actively seeking a collaborating laboratory to pursue the experimental validation ladder the paper lays out; this does not affect the integrity or availability of the computational work presented here.

I believe this framework offers your readership a reproducible route from real spatial data to prioritized, low-toxicity control strategies for PDAC, together with a transparent account of where such an approach holds and where it breaks down. Thank you for considering the manuscript; I would be glad to respond to any questions.

Sincerely,

Seung-Il Kim
Independent Researcher
kusi81kim@gmail.com · ORCID 0009-0007-5965-9212

---

*Suggested reviewers (optional): [to be completed — 3–5 names with expertise in agent-based tumor modeling, PDAC CAF biology, and adaptive/evolutionary therapy; avoid conflicts].*
*Opposed reviewers (optional): [none / to be completed].*
