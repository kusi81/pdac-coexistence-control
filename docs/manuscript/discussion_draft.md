# Discussion (draft v1)

> **작성 메모:** IMRaD 초고 완성 절. 전부 in-silico/hypothesis 프레이밍 유지. 한계 정직
> (검증부재·2D·현상학적·주석). [n]=Introduction 참조. §=Results 절.

---

## 4. Discussion

### 4.1 A framework that bridges three disconnected research lines
We present an in-silico framework that connects three individually mature but
mutually disconnected bodies of work—static network pharmacology of food-medicine-
homology compounds, wet-laboratory studies of herbal modulation of PDAC CAFs, and
spatial/agent-based modeling of the tumor microenvironment. Our systematic survey
quantified this gap: single-axis literature is abundant, yet queries requiring all
three axes returned only off-target hits (§3.1, Fig. 1). The contribution is
therefore integrative—a spatially grounded, adaptive-control model in which
food-medicine-homology compounds are encoded as mechanistic perturbations and
optimized against a coexistence objective that treats the myCAF barrier as a
resource rather than an obstacle.

### 4.2 Biological plausibility
Several results align with independent biology. On author-annotated CosMx data,
our metric reproduced the myCAF-proximal/iCAF-distal architecture in 15 of 16
tumors and showed the myCAF barrier tightening under chemoradiation while
cytotoxic T cells remained excluded (§3.3)—consistent with the paradoxical
worsening of PDAC after stromal ablation [6,7] and with the reframing of myCAF as
a containment structure rather than a purely pro-tumor element. In simulation,
adaptive scheduling achieved durable coexistence at a fraction of the toxicity of
continuous dosing (§3.4), and our sensitivity analysis showed that the resistance
fitness cost governs competitive suppression of resistant clones (§S3), mirroring
the adaptive-therapy rationale [9–12].

### 4.3 Positioning relative to prior art
The framework advances each of the three neighboring camps on the axis it lacks:

| Prior camp | What it provides | What it lacks | Our addition |
|---|---|---|---|
| Static network pharmacology of food-medicine-homology compounds [14,16] | compound–target–pathway maps | time, space, scheduling | dynamic, spatial, adaptive-control simulation |
| Wet-lab herbal–CAF studies in PDAC [17,18] | real CAF-modulating biology | prediction/optimization | in-silico ranking of combinations, sequences, doses |
| Spatial/ABM models of PDAC TME [19,20] | tumor–stroma–immune dynamics | natural compounds; coexistence objective | food-medicine-homology perturbations + control objective |

Crucially, we do not claim novelty in *identifying* molecular targets of these
compounds—that space is saturated (§3.1)—but in embedding that knowledge in a
dynamic, spatial, control-theoretic framework.

### 4.4 Translational implications
The design principle that emerges is control-oriented: pair a low-toxicity
anti-CAF backbone with a resistance-sparing anti-proliferative partner rather than
maximizing cytotoxic pressure (§3.5). In-silico, food-medicine-homology regimens
maintained tumor control at a small fraction of the toxicity of continuous
gemcitabine, and a sub-maximal, adaptively timed dose further reduced burden
(§3.6). For a disease in which conventional therapy imposes severe toxicity for
limited benefit, a strategy oriented toward slowing proliferation, exploiting the
myCAF barrier, and preserving quality of life—rather than eradication—may be a
more realistic near-term goal, potentially extending functional survival and
easing symptom burden. We emphasize that these are model-generated hypotheses, and
that accessibility and low toxicity are advantages only if efficacy is confirmed
experimentally and clinically.

### 4.5 Limitations
Our study has important limitations that bound its claims. First and foremost, it
is entirely in silico: outputs are ranked, testable hypotheses, not evidence of
efficacy, and require experimental validation. Second, the ABM is phenomenological
with literature-inferred parameters rather than data-fitted; although the adaptive
advantage was robust to ±50% parameter variation (§S3), outcomes were most
sensitive to tumor proliferation and immune-barrier parameters, which therefore
require empirical calibration, and control was not guaranteed under the most
aggressive proliferation regimes. Third, our spatial analyses are of 2D sections;
true tissue architecture is three-dimensional, and the myCAF-proximity positive
control failed on our targeted-panel Xenium annotation, succeeding only with
full-panel author annotations (§3.2)—so real-data spatial claims rest on the CosMx
cohort and remain associations in a modest sample (n = 15). Fourth, simulations
used a single synthetic architecture and an arbitrary-unit toxicity scale that is
not clinically calibrated. Fifth, even low-toxicity natural combinations can select
for resistance under aggressive suppression (competitive release; §3.5). Finally,
the novelty survey is PubMed-restricted, and the molecular grounding for most
food-medicine-homology compounds rests on pathway-level mechanistic inference
rather than solved compound-bound structures (§3.7).

### 4.6 Future directions
Several steps follow directly. Experimentally, the ranked regimens define a
low-cost validation ladder—2D PDAC cell-line dose–response and combination/synergy
assays, then patient-derived organoid screening—prioritizing the anti-CAF-backbone
combinations and the competitive-release predictions. Computationally, priorities
are: parameter calibration against experimental rates (especially proliferation
and immune-barrier terms); global sensitivity analysis (e.g., Sobol) to capture
interactions beyond the one-at-a-time sweep; a 3D registration pipeline over serial
sections to replace 2D metrics; annotation transfer from full-transcriptome
references (or continued use of large-panel author-annotated data) to make CAF
subtyping robust; and molecular-dynamics or larger-scale virtual screening to
strengthen the molecular tier. Broadening the organ context (e.g., the HCC "co-
opted stroma" regime) and extending the survey to Scopus/Web of Science under a
PRISMA workflow would further generalize and harden the framework.

### 4.7 Conclusion
Within a control—rather than eradication—objective, this work provides a
reproducible computational path from real spatial data to prioritized,
low-toxicity, food-medicine-homology regimens for PDAC coexistence. By treating
the myCAF barrier as a resource and resistance as an evolutionary constraint, it
reframes the therapeutic goal toward durable tumor containment, and it yields
concrete, testable hypotheses to guide the experimental work that must follow.

---

## 편집 체크리스트
- [ ] 4.3 표를 Introduction 3진영 문단·§3.1과 용어 일치
- [ ] 임상함의(4.4)에서 과대주장 없는지 재검(전부 model-generated 명시됨)
- [ ] 한계(4.5) 항목이 Results 각 절과 1:1 대응하는지 확인
- [ ] future(4.6) wet-lab 사다리 = 별도 로드맵 문서와 정합
- [ ] 인용번호 최종 대조
