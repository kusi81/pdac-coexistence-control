# Introduction (draft v2 — 인용 확정본)

> **작성 메모 (한글):** 국제 게재 목표 → 영문. 범위는 정직하게 "실공간데이터 접지
> in-silico *가설 생성·우선순위화* 프레임워크"로 한정(검증 논문 아님). 체계검색(§novelty_assessment)
> 에서 확인된 실제 선행연구를 인용해 3진영 차별화를 5번째 문단에 명시 → 리뷰어 반박 선제 차단.
> **v2 변경:** 전 인용 PMID·DOI를 PubMed E-utilities로 확정(pipeline/verify_citations*.py),
> 본문 번호 재정렬, ⚠️ 미확인 표시 제거. [13]만 arXiv 프리프린트(게재본 확인 여지).
> 타깃 저널 후보: *PLOS Computational Biology*, *npj Systems Biology and Applications*,
> *Cancers*, *Frontiers in Pharmacology/Oncology*.

---

## Introduction

Pancreatic ductal adenocarcinoma (PDAC) remains among the most lethal human
malignancies, with a five-year survival of roughly 13% and a median survival
measured in months for advanced disease [1]. A defining feature of PDAC is its
extensive desmoplastic stroma, in which cancer-associated fibroblasts (CAFs) and
a dense extracellular matrix can constitute the majority of the tumor mass [2].
For patients, the clinical reality is not only high mortality but a steep, often
rapid decline in remaining quality-adjusted life: the pace at which PDAC erodes
patient time and comfort is itself a therapeutic target. Interventions that
meaningfully *slow* tumor expansion—even without eradication—could translate
into reduced symptom burden and extended functional survival.

CAFs in PDAC are not a monolithic population. Single-cell and spatial studies
have resolved at least three functional states—myofibroblastic CAFs (myCAF),
inflammatory CAFs (iCAF), and antigen-presenting CAFs (apCAF)—that occupy
distinct spatial niches and are, importantly, plastic and interconvertible
[3,4,5]. myCAFs, enriched immediately adjacent to tumor epithelium, deposit the
fibrotic matrix that both supports and physically constrains the tumor. This
duality has a critical therapeutic corollary: attempts to *ablate* the stroma in
preclinical PDAC paradoxically accelerated disease and worsened survival [6,7],
reframing the field's question from "should CAFs be deleted?" to "should they be
reeducated?" [8]. We take this further and treat the myCAF-derived barrier not as
an obstacle to be removed but as a *controllable resource*—a containment
structure that, if strategically preserved and modulated, may be leveraged to
restrain rather than merely attack the tumor.

This perspective aligns naturally with a shift in therapeutic philosophy from
maximum-tolerated-dose (MTD) eradication toward evolutionary *control*. Adaptive
therapy, which maintains a population of treatment-sensitive cells to
competitively suppress resistant subclones through intermittent or modulated
dosing, has shown that stabilizing tumor burden can outperform attempts to
minimize it, provided resistant cells carry a fitness cost [9,10,11,12]. Related
control-theoretic analyses now frame tumor management explicitly in terms of
*controllability* under co-evolving resistance [13]. In PDAC—where MTD regimens
impose severe toxicity for limited benefit—an approach oriented toward durable
coexistence, low toxicity, and preserved quality of life is especially
compelling, yet the spatial dimension of CAF-mediated containment has not been
integrated into such control frameworks.

A complementary opportunity lies in the compounds available for long-term,
low-toxicity control. "Food-medicine homology" (medicine-food homology)
substances—plant and animal materials used in East Asian traditions as both food
and remedy (e.g., garlic, ginseng, Platycodon, mugwort, hawthorn)—offer
multi-target activity with favorable safety and accessibility, making them
attractive backbones for chronic control regimens rather than acute cytotoxic
bursts [14,15]. However, the computational study of these compounds in oncology
has been dominated by *static* network pharmacology: compound–target–pathway
mapping and molecular docking that predict mechanism but do not simulate
dynamics, space, or scheduling [14,15,16]. Our own systematic PubMed survey
(Methods) confirmed this imbalance quantitatively—food-medicine-homology
network-pharmacology cancer studies are abundant, whereas their integration with
dynamic, spatial, or adaptive-control modeling is essentially absent.

Indeed, the relevant prior art separates into three camps that, to our
knowledge, have not been bridged. First, static network-pharmacology studies of
food-medicine-homology compounds against cancer establish molecular targets but
remain time- and space-agnostic [16]. Second, wet-laboratory studies show that
specific herbal compounds and formulas modulate PDAC CAFs—suppressing CAF
activation [17] and reversing CAF-induced gemcitabine resistance [18]—providing
biological grounding but no predictive or optimization framework. Third,
agent-based and multiscale models of the tumor microenvironment, increasingly
informed by spatial transcriptomics, capture tumor–stroma–immune dynamics but do
not incorporate natural-compound perturbations or a coexistence objective
[19,20]. No existing work, to our knowledge, connects the compound knowledge of
the first camp to the dynamic spatial models of the third in service of the
CAF-containment biology of the second.

Here we present an in-silico framework that occupies precisely this gap. We build
a spatial agent-based model of the PDAC tumor–myCAF–immune ecosystem, grounded in
real single-cell spatial transcriptomics (Xenium) [20], in which
food-medicine-homology compounds are encoded as mechanistic parameter
perturbations acting across a molecular-to-cellular-to-tissue hierarchy. Rather
than optimizing for eradication, we define a *coexistence-control* objective that
exploits the myCAF barrier and explicitly models resistance dynamics and a
resistance fitness cost, and we use it to prioritize combinations, sequences, and
doses under adaptive scheduling. We emphasize that the framework is a
hypothesis-*generating* and hypothesis-*prioritizing* engine: its outputs are
testable predictions—candidate compound combinations and schedules ranked by
predicted control and toxicity—intended to focus, not replace, subsequent
experimental validation.

The remainder of this paper is organized as follows. We first describe the
spatial-transcriptomic grounding and the barrier/containment metrics used to
characterize CAF-mediated tumor constraint. We then detail the agent-based
control model, its resistance dynamics, and the encoding of food-medicine-
homology compounds. Finally, we present in-silico experiments comparing
continuous versus adaptive schedules and single-agent versus combination and
sequential regimens, and we discuss the resulting testable hypotheses together
with the framework's current limitations, including the 2D-section constraint and
the need for experimental parameter calibration.

---

## References (PMID·DOI 확정 — PubMed E-utilities 검증)

1. Siegel RL, Kratzer TB, et al. Cancer statistics, 2025. *CA Cancer J Clin*
   2025;75(1):10-45. PMID: 39817679. doi:10.3322/caac.21871
2. Hosein AN, Brekken RA, Maitra A. Pancreatic cancer stroma: an update on
   therapeutic targeting strategies. *Nat Rev Gastroenterol Hepatol*
   2020;17(8):487-505. PMID: 32393771. doi:10.1038/s41575-020-0300-1
3. Öhlund D, Handly-Santana A, Biffi G, et al. Distinct populations of
   inflammatory fibroblasts and myofibroblasts in pancreatic cancer. *J Exp Med*
   2017;214(3):579-596. PMID: 28232471. doi:10.1084/jem.20162024
4. Elyada E, Bolisetty M, et al. Cross-species single-cell analysis of pancreatic
   ductal adenocarcinoma reveals antigen-presenting cancer-associated
   fibroblasts. *Cancer Discov* 2019;9(8):1102-1123. PMID: 31197017.
   doi:10.1158/2159-8290.CD-19-0094
5. Biffi G, Oni TE, et al. IL1-induced JAK/STAT signaling is antagonized by TGFβ
   to shape CAF heterogeneity in pancreatic ductal adenocarcinoma. *Cancer Discov*
   2019;9(2):282-301. PMID: 30366930. doi:10.1158/2159-8290.CD-18-0710
6. Rhim AD, Oberstein PE, et al. Stromal elements act to restrain, rather than
   support, pancreatic ductal adenocarcinoma. *Cancer Cell* 2014;25(6):735-747.
   PMID: 24856585. doi:10.1016/j.ccr.2014.04.021
7. Özdemir BC, Pentcheva-Hoang T, Carstens JL, et al. Depletion of
   carcinoma-associated fibroblasts and fibrosis induces immunosuppression and
   accelerates pancreas cancer with reduced survival. *Cancer Cell*
   2014;25(6):719-734. PMID: 24856586. doi:10.1016/j.ccr.2014.04.005
8. Qu C, Wang Q, et al. Cancer-associated fibroblasts in pancreatic cancer:
   should they be deleted or reeducated? *Integr Cancer Ther* 2018;17(4):1016-1019.
   PMID: 30136592. doi:10.1177/1534735418794884
9. Gatenby RA, Silva AS, Gillies RJ, Frieden BR. Adaptive therapy. *Cancer Res*
   2009;69(11):4894-4903. PMID: 19487300. doi:10.1158/0008-5472.CAN-08-3658
10. Zhang J, Cunningham JJ, Brown JS, Gatenby RA. Integrating evolutionary
    dynamics into treatment of metastatic castrate-resistant prostate cancer.
    *Nat Commun* 2017;8(1):1816. PMID: 29180633. doi:10.1038/s41467-017-01968-5
11. West J, You L, et al. Towards multidrug adaptive therapy. *Cancer Res*
    2020;80(7):1578-1589. PMID: 31948939. doi:10.1158/0008-5472.CAN-19-2669
12. McGehee C, Mori Y. A mathematical framework for comparison of intermittent
    versus continuous adaptive chemotherapy dosing in cancer. *NPJ Syst Biol Appl*
    2024;10(1):140. PMID: 39614108. doi:10.1038/s41540-024-00461-2
13. [Preprint] Controllability in a class of cancer therapy models with
    co-evolving resistance. *arXiv*:2503.07533, 2025. (게재본 확인 필요)
14. Zheng J, Wu M, et al. Network pharmacology to unveil the biological basis of
    health-strengthening herbal medicine in cancer treatment. *Cancers (Basel)*
    2018;10(11):461. PMID: 30469422. doi:10.3390/cancers10110461
15. Zhou Z, Nan Y, et al. Hawthorn with "homology of medicine and food": a review
    of anticancer effects and mechanisms. *Front Pharmacol* 2024;15:1384189.
    PMID: 38915462. doi:10.3389/fphar.2024.1384189
16. Zhao X, Xiu J, et al. Network pharmacology and bioinformatics study of six
    medicinal food homologous plants against colorectal cancer. *Int J Mol Sci*
    2025;26(3):930. PMID: 39940699. doi:10.3390/ijms26030930
17. Chen L, Qu C, et al. Chinese herbal medicine (Qingyihuaji) suppresses
    invasion-promoting capacity of cancer-associated fibroblasts in pancreatic
    cancer. *PLoS One* 2014;9(4):e96177. PMID: 24781445.
    doi:10.1371/journal.pone.0096177
18. Hu X, Peng X, et al. Shikonin reverses cancer-associated fibroblast-induced
    gemcitabine resistance in pancreatic cancer cells by suppressing MCT4-mediated
    reverse Warburg effect. *Phytomedicine* 2024;123:155214. PMID: 38134861.
    doi:10.1016/j.phymed.2023.155214
19. Norton KA, Gong C, et al. Multiscale agent-based and hybrid modeling of the
    tumor immune microenvironment. *Processes (Basel)* 2019;7(1):37.
    PMID: 30701168. doi:10.3390/pr7010037
20. Shiau C, Cao J, et al. Spatially resolved analysis of pancreatic cancer
    identifies therapy-associated remodeling of the tumor microenvironment.
    *Nat Genet* 2024;56(11):2466-2478. PMID: 39227743.
    doi:10.1038/s41588-024-01890-9

---

## 편집 체크리스트 (v2)
- [x] 전 인용 PMID·DOI 확정 (PubMed 검증) — [13] arXiv만 프리프린트
- [ ] [13] controllability 게재본(저널) 여부 최종 확인 → 있으면 교체
- [ ] 타깃 저널 확정 후 인용 스타일(Vancouver/저널양식) 변환
- [ ] Abstract·Methods와 용어 통일(myCAF, food-medicine homology, coexistence-control)
- [ ] bioRxiv/arXiv 인접연구 추가 인용 여부 검토(적응요법 threshold 등)
- [ ] Xenium 데이터 출처(GSE274673) Methods에 명기, 본문 [20]과 구분
