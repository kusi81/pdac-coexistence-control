# Abstract (draft v1)

> **작성 메모 (한글):** Introduction v2와 용어 통일(myCAF, food-medicine homology,
> coexistence-control, resistance fitness cost, Xenium, hypothesis-generating).
> 범위 정직: in-silico 가설 *생성·우선순위화* — 실험 검증본 아님. 정량 수치는 모델
> *예측*으로 명시(과대주장 금지). 저널별로 (A) 일반형 ~250단어 / (B) 구조형 중 택1.
> PLOS Comput Biol 투고 시 Author Summary 필수 → (C) 포함.

---

## (A) 일반형 초록 (unstructured, ~250 words)

Pancreatic ductal adenocarcinoma (PDAC) combines dismal survival with a dense,
myofibroblastic cancer-associated fibroblast (myCAF) stroma that both supports
and physically constrains the tumor. Because attempts to ablate this stroma have
paradoxically accelerated disease, we ask whether the myCAF barrier can instead
be treated as a *controllable resource* for restraining, rather than eradicating,
the tumor. We present an in-silico framework that integrates three previously
disconnected research lines—static network pharmacology of "food-medicine
homology" compounds, wet-laboratory studies of herbal modulation of PDAC CAFs,
and agent-based modeling of the tumor microenvironment. A systematic PubMed
survey confirmed that, while each line is individually well populated, their
integration is essentially absent. Our framework builds a spatial agent-based
model of the PDAC tumor–myCAF–immune ecosystem, grounded in single-cell spatial
transcriptomics (Xenium), and encodes food-medicine-homology compounds (e.g.,
garlic, ginseng, Platycodon, mugwort) as mechanistic parameter perturbations
across a molecular-to-cellular-to-tissue hierarchy. Rather than optimizing for
eradication, we define a coexistence-control objective that leverages the myCAF
barrier and explicitly represents resistance dynamics with a resistance fitness
cost. In simulation, adaptive schedules achieved tumor containment comparable to
continuous maximum-dose regimens at substantially lower predicted toxicity, and
compound combinations with sequential myCAF-modulating and anti-proliferative
phases outperformed single agents, with anti-CAF components predicted to minimize
resistance emergence. We stress that these outputs are testable hypotheses—ranked
candidate combinations and schedules—intended to focus, not replace, experimental
validation. The framework offers a reproducible route from real spatial data to
prioritized, low-toxicity control strategies for PDAC coexistence.

---

## (B) 구조형 초록 (structured — Cancers 등)

**Background:** Pancreatic ductal adenocarcinoma (PDAC) is defined by a
myofibroblastic cancer-associated fibroblast (myCAF) stroma that both supports and
physically constrains the tumor; stromal ablation has paradoxically worsened
outcomes, motivating strategies that exploit rather than remove the barrier.

**Methods:** We developed an in-silico framework combining (i) a spatial
agent-based model of the PDAC tumor–myCAF–immune ecosystem grounded in
single-cell spatial transcriptomics (Xenium), (ii) barrier/containment metrics
via matched-null spatial tests, and (iii) mechanistic encoding of "food-medicine
homology" compounds as multiscale parameter perturbations. A coexistence-control
objective with explicit resistance dynamics and a resistance fitness cost was
used to rank combinations, sequences, and adaptive schedules. Novelty was
assessed by a systematic PubMed survey.

**Results:** The systematic survey confirmed that static network pharmacology of
food-medicine-homology compounds, wet-lab herbal–CAF studies, and microenvironment
modeling are each well populated but never integrated. In simulation, adaptive
scheduling matched continuous maximum-dose containment at markedly lower predicted
toxicity; sequential myCAF-modulating/anti-proliferative combinations outperformed
single agents; and anti-CAF components minimized predicted resistance.

**Conclusions:** The framework converts real spatial data into prioritized,
testable, low-toxicity control hypotheses for PDAC coexistence. Outputs are
hypothesis-generating and require experimental validation.

---

## (C) Author Summary (PLOS Comput Biol 양식, ~150–200 words, 비전문가 대상)

Pancreatic cancer is one of the deadliest cancers, and patients often lose
healthy time very quickly. Most treatments aim to destroy the tumor outright, but
this frequently fails and carries heavy toxicity. The tumor is wrapped in a dense
"scar-like" tissue built by cells called myofibroblastic cancer-associated
fibroblasts (myCAFs). Surprisingly, experiments that tried to remove this tissue
made the cancer worse—suggesting the barrier also holds the tumor back. We asked a
different question: instead of destroying the tumor, can we *manage* it—slowing its
growth and using the myCAF barrier to keep it contained, so patients feel better
and live longer? To explore this safely and cheaply, we built a computer
simulation of pancreatic tumor tissue based on real spatial measurements of where
each cell sits. Into this simulation we added compounds that are eaten as food and
used as traditional medicine in East Asia. The model lets us test many
combinations, timings, and doses in silico and rank the most promising, low-toxicity
strategies. These predictions are hypotheses meant to guide future laboratory
experiments, not to replace them.

---

## Keywords
pancreatic ductal adenocarcinoma; cancer-associated fibroblasts; myCAF; adaptive
therapy; tumor coexistence; agent-based model; spatial transcriptomics;
food-medicine homology; network pharmacology; drug resistance

---

## 편집 체크리스트
- [ ] 타깃 저널 확정 → (A) vs (B) 선택, 단어수 상한 맞추기(예: PLOS 250, Cancers 200 구조형)
- [ ] 정량 수치를 넣을지 결정: 넣는다면 Results의 확정 시뮬 수치(독성 비·TTP)로 교체하고
      "predicted/in-silico" 수식어 유지. 지금은 방향성만 서술(과대주장 회피).
- [ ] Xenium 시료(GSE274673) 명기 위치는 Methods(초록엔 플랫폼명만)
- [ ] Introduction/Discussion과 핵심 용어 최종 통일
