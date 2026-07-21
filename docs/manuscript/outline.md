# PDAC 논문 전체 아웃라인 (manuscript skeleton)

> **범위:** in-silico 가설 생성·우선순위화 프레임워크 (검증본 아님).
> **상태 태그:** ✅자료있음(코드·결과 존재) · 🔬실행필요(스크립트는 있으나 재실행/집계) ·
> ✍️작성필요(글) · ⚠️정직고지(한계).
> 연결: [Introduction v2](introduction_draft.md) · [Abstract](abstract_draft.md) ·
> [Novelty](../literature_search/novelty_assessment.md)

---

## Title (후보)
- **주:** "A spatially grounded agent-based framework for coexistence-control of
  pancreatic cancer using food-medicine-homology compounds"
- 대안1: "Modeling the myCAF barrier as a controllable resource: an in-silico
  framework for low-toxicity coexistence therapy in PDAC"
- 대안2: "From static network pharmacology to dynamic spatial control:
  food-medicine-homology compounds for PDAC tumor containment in silico"

## Running title / Authors / Affiliations
- ✍️ 저자: 본인(Independent Researcher) + 협업 PI 확보 시 추가 → [[pdac-paper-funding-track]] 전략
- Corresponding author, ORCID, 연락처

---

## Abstract + Keywords
- ✅ [abstract_draft.md](abstract_draft.md) (A 일반형 / B 구조형 / C Author Summary)

---

## 1. Introduction
- ✅ [introduction_draft.md](introduction_draft.md) v2 (인용 20개 확정)
- 골격: PDAC 치명성·myCAF 이중성 → 장벽=통제자원 → 박멸→통제 패러다임 →
  약식동원 저독성 백본 → **3진영 미연결(체계검색)** → 본 기여 → 논문 구성

---

## 2. Methods

### 2.1 Systematic novelty survey
- ✅ [lit_search.py](../../pipeline/lit_search.py), [verify_citations*.py](../../pipeline/)
- 내용: PubMed E-utilities, 10 불린 쿼리(교집합 차원), 132편, 통합쿼리 A=2·B=3(off-target)
- ⚠️ 한계: PubMed only; Scopus/WoS·PRISMA는 §4.5/future

### 2.2 Spatial data & cell annotation
- ✅ [data_loader.py](../../pipeline/data_loader.py): load_xenium_bundle(µm coords),
  annotate_pdac (2단계 module-score: PDAC_COARSE_MODULES + CAF_SUBTYPE_MODULES, FertigLab)
- 시료: **GSE274673 Xenium PDAC 6종** (naive 446/449/450, CRT 447/448/451) — XENIUM_TREATMENT
- 보조: MIBI-TOF(squidpy), Zhou Visium(살베이지 복구) — 대조/한계 예시
- ✍️ 전처리 파라미터·QC 기준 표

### 2.3 Barrier / containment metrics
- ✅ [spatial_core.py](../../pipeline/spatial_core.py):
  barrier_score(matched-null corridor density, tumor_shell_um),
  rim_enrichment(peritumoral shell, n_perms), radial_profile, proximity_test,
  median_nn_distance
- ✍️ 수식·귀무모형(matched-null) 정의, 유의성 검정(permutation) 서술

### 2.4 Agent-based control model
- ✅ [abm.py](../../pipeline/abm.py): TumorABM (tumor/CAF/immune agents)
- 핵심 규칙: 증식 kp, CAF 장벽 효과(k_caf_death, caf_boost_ref), 면역 사멸
- **내성 동역학:** self.res 배열, init_resistant_frac=0.01, mutation_rate=0.001,
  resistant_immune_evasion=0.45, **resistance_cost=0.24**(문헌 접지),
  kp_res=k_prolif·(1−resistance_cost); 약물은 sensitive만, 면역사멸은 evasion 감쇠
- 맥락: CONTEXTS(pdac caf_protumor=0 / hcc 0.9) — context_params()
- ✍️ 규칙 의사코드(pseudocode) 박스, 상태전이 다이어그램

### 2.5 Food-medicine-homology compound encoding
- ✅ abm.py SUBSTANCES(17종: curcumin, ginsenoside_rg3, garlic, mugwort,
  wild_ginseng, platycodon, sea_cucumber, deer_antler, danshen, milk_thistle,
  astragaloside, gemcitabine, erlotinib, entecavir[가설] + generics)
- 각 물질: 표적축(anti-prolif / anti-CAF / immune), EVIDENCE_BADGE 근거수준, TOXICITY
- compose_regimen, classify_substances, build_cycle_schedule
- ✍️ 물질-파라미터 매핑 표(부록 S1), 근거수준 기준 명시

### 2.6 Multiscale molecular grounding
- ✅ [mol_page.py](../../pipeline/mol_page.py): PDB 근거 3계층(experimental/docking/mechanistic),
  검증 PDB(1M17 EGFR+에를로티닙, 1P62 dCK+젬시타빈, 5IKR COX-2, 1BG1 STAT3 …)
- ✅ [render_compare_structures.py](../../pipeline/render_compare_structures.py)
- ⚠️ 정직고지: docking/mechanistic은 참조리간드/예측잔기 — "초록=우리물질"은 experimental만
- ✍️ 근거계층이 결론 신뢰도를 어떻게 가르는지 서술

### 2.7 Control objective & metrics
- ✅ abm.py control_metrics(): ttp_days, progression_censored, final_resistant_frac,
  cum_toxicity, control_score; adaptive on/off(adapt_on/adapt_off)
- ✍️ coexistence-control 목적함수 정의(박멸 아님), TTP·독성 정의, 순위기준
  (censored-first → lowest toxicity; frontier: TTP vs toxicity)

### 2.8 In-silico experimental design
- ✅ [sim_page.py](../../pipeline/sim_page.py)(동시/순차사이클/적응형),
  [opt_page.py](../../pipeline/opt_page.py) + optimize_dose_band / optimize_band_aggr(SEEDS=[42,7,123])
- 실험: (E1) continuous vs adaptive; (E2) single vs combination vs sequential;
  (E3) dose×band grid; (E4) context pdac vs hcc(부록)
- ✍️ 시드·반복·파라미터 스윕 명세

### 2.9 Implementation & reproducibility
- ✅ Python 3.13, squidpy 1.8.3/scanpy/anndata, Streamlit 대시보드(app.py, 4 모드)
- ✍️ Code availability(리포 URL), 데이터 출처(GSE274673), 환경 고정(requirements.txt)

---

## 3. Results

### 3.1 The integration gap is real (novelty landscape)
- ✅ 표: 차원별 count(개별 축 포화 D5=33/D7=31 vs 통합 A=2·B=3 off-target)
- **Figure 1**: 차원별 히트 막대 + 3진영 벤다이어그램(빈 교집합 강조)

### 3.2 Barrier metric validation — 지표는 정상, 주석이 관건 (2-데이터셋)
- ✅ 합성 검증: contained z≈+22.3 vs diffuse z≈+0.9 (특이도)
- ✅ **양성대조 재현(SCOTIA 저자주석 CosMx, scotia_posctrl.py):** myCAF vs iCAF→Malignant
  proximity **15/16 통과** = SCOTIA(ref[20]) 재현 → **우리 지표 정당성 확립**.
- ⚠️ **정직한 대조 — Xenium 주석 실패:** 동일 지표가 Xenium(GSE274673) module-score
  주석에선 양성대조 실패(diag_mycaf.py/refine_annotate.py). 진단: 480 타깃패널 +
  argmax 타이핑 한계(골드마커 POSTN/LRRC15 넣어도 4/5 실패) → **지표 아닌 주석 문제**.
  MIBI-TOF 대장=장벽없음(위양성 회피), Zhou Visium=spot해상도 한계.
- **핵심 메시지:** 좋은 주석엔 지표가 생물학을 복원, 나쁜 주석엔 실패 → robust CAF
  공간분석은 전체패널/저자주석 필요. (프레임워크 "가설생성" 규정과 일관)
- **Figure 2**: 합성 contained/diffuse + SCOTIA 양성대조(15/16) + Xenium 실패 대조

### 3.3 Real PDAC spatial characterization — Untreated vs CRT (SCOTIA 완성 F3)
- ✅ **완성(scotia_rim.py, assets/rim_scotia.png):** SCOTIA 저자주석 15시료
  (Untreated9/CRT6) 종양-rim(30µm) 세포농축.
- 핵심 결과: **myCAF 강한 종양인접(+9.1→+15.7), CRT로 더 조여짐** = 장벽 논지 실증;
  iCAF 원위(−3.5→−9.2); **CD8+ T는 양쪽 다 배제(−5.6→−4.0) — CRT로도 rim 진입 못함
  (면역배제 지속)**; Macrophage 인접(+6.8); Plasma 강배제(−11.3); Pericyte CRT후 인접(+3.3).
- 서사: myCAF 장벽 + CD8 배제가 치료 하에서 유지·강화 → "장벽=통제자원" 직접 지지.
- **Figure 3**: rim_scotia.png (양성대조 내장). ⚠️출판전 사소한 글리프(−) 수정.
- ⚠️ n=15 소규모·저자주석 재사용(순환 아님=양성대조 목적) 명시.

### 3.4 Adaptive scheduling controls at lower toxicity
- ✅ continuous(박멸, 고독성) vs adaptive(공존, 저독성); resistance_cost가 가능케 함
- **Figure 4**: 종양·내성 궤적(continuous vs adaptive) + 독성 누적 비교

### 3.5 Combination & sequential regimens
- ✅ myCAF 약화(ginsenoside/garlic/mugwort) ↔ anti-prolif(erlotinib) 순차사이클;
  anti-CAF가 내성 최소; anti-fibrotic 단독은 불충분(anti-prolif 파트너 필요)
- **Figure 5**: 단일 vs 조합 vs 순차 — 통제점수/독성/내성분율

### 3.6 Dose×band optimization
- ✅ 최적 용량 80% / off 0.5 (독성 27 < 전용량 32); frontier(TTP vs toxicity)
- **Figure 6**: dose×band 히트맵(독성·최종종양) + 최적점

### 3.7 Multiscale molecular hypotheses
- ✅ 표적 결합 근거계층표; 에를로티닙/커큐민/젬시타빈 구조 비교
- **Figure 7**(부록 가능): 구조 비교 + 근거계층 캡션

---

## 4. Discussion

### 4.1 Principal contribution — the bridge
- ✍️ 3진영을 잇는 최초의 동적·공간·적응통제 프레임(정적 netpharm→동적 확장)

### 4.2 Biological plausibility
- ✍️ myCAF 이중성·기질제거 역설[6,7]과 정합; adaptive therapy 경쟁·적합도비용[9–12]

### 4.3 Positioning vs prior art
- ✍️ [16](정적 약식동원 netpharm)·[17,18](한약-CAF wet-lab)·[19,20](모델)과 1:1 차별화표

### 4.4 Translational implications
- ✍️ 공존·삶의질·저독성·약식동원 접근성; 잔여수명·고통경감(연구동기 반영)

### 4.5 Limitations (정직)
- ⚠️ 2D 단면(3D 아님); 파라미터 다수 문헌유추(민감도분석 필요);
  근접성 양성대조 실패; PubMed-only 서베이; in-silico only(검증 부재);
  표본 소규모(n=6); docking/mechanistic 구조의 해석 한계

### 4.6 Future work
- ✍️ 연속절편 3D 정합 파이프라인; wet-lab 검증 사다리(2D 세포주→오가노이드,
  [[pdac-paper-funding-track]]); MD 시뮬/GPU 확장; Scopus/WoS + PRISMA; 민감도분석

---

## 5. Conclusion
- ✍️ 실공간데이터 → 우선순위화된 저독성 통제 가설. 박멸 아닌 공존의 계산적 경로.

---

## Back matter
- Data availability: GSE274673, MIBI-TOF(squidpy), 코드 리포
- Code availability: pipeline/ + Streamlit app
- Author contributions / Funding(현재 없음 → 향후 국가R&D) / Conflicts(none)
- ✍️ Supplementary: S1 물질-파라미터표, S2 전체 쿼리·결과, S3 민감도분석, S4 대시보드 스냅샷

---

## Figures 요약 (7개)
| # | 내용 | 상태 |
|---|---|---|
| F1 | 신규성 landscape(막대+벤) | ✅**조립완료** assets/fig1_novelty.png (fig1.py) |
| F2 | 지표 검증(합성/SCOTIA통과/Xenium실패) | ✅**조립완료** assets/fig2_validation.png (fig2.py) |
| F3 | Untreated vs CRT rim 조성 | ✅**완료** assets/rim_scotia.png (scotia_rim.py) |
| F4 | adaptive vs continuous 궤적·독성 | ✅ |
| F5 | 단일/조합/순차 비교 | ✅ |
| F6 | dose×band 최적화 히트맵 | ✅ |
| F7 | 분자 구조·근거계층(부록) | ✅ |

## 진행 로그
- [x] F3(SCOTIA) 완료 → Fig 3 (rim_scotia.png)
- [x] Fig 2 조립 (fig2_validation.png)
- [x] Results §3.2–3.3 본문 (results_3_2_3_3_draft.md)
- [x] **민감도분석 S3 완료** → Fig S3 (sensitivity.png), S3_sensitivity_draft.md.
      핵심: 적응형 저독성 이점 강건(cost 비의존), resistance_cost=내성억제 역할,
      민감도는 종양-면역 균형(k_prolif·barrier·k_kill)에 집중. 한계: 극단 체제서 통제 저하.
- [x] **Results §3.4–3.7 본문** (results_3_4_3_7_draft.md). 실행값 확정: 3.4 적응형
      0.82x/독성25 vs 연속 0.00x/독성120; 3.5 자연물 독성5–38 vs 젬128, 마늘/산삼 CS25/24,
      마늘+쑥 0.4x이나 내성0.34(경쟁적방출), 커큐민+마늘+Rg3 내성0.01, 항섬유 단독 부족;
      3.6 최적 용량0.8/off0.5 독성27; 3.7 근거3계층(실험 1P62/1M17만, 나머지 도킹/기전).

- [x] **Fig 1 + §3.1 완료** (fig1_novelty.png, results_3_1_draft.md) → **Results §3.1–3.7 전체 완결**
- [x] **Methods §2.1–2.9 완료** (methods_draft.md)
- [x] **Discussion §4.1–4.7 완료** (discussion_draft.md, 3진영 차별화표 포함) → **IMRaD 초고 본문 완결**

## 남은 마감 작업
- [x] **단일 manuscript.md 병합** (IMRaD 전체+그림캡션+참조20+마스터체크리스트)
- [x] **Supp Table S1** (S1_parameters.md, 24 파라미터+맥락 오버라이드+접지인용 G1-G3)
1. ✍️ 타깃 저널 선정 + Abstract 수치 삽입(3.4–3.6) + 인용스타일 변환
2. ⚠️ [G3] Science Adv abm7212 서지 확정 · [13] arXiv 게재본 확인 · repo URL · 저자블록
3. ⚠️ 인용번호 전수 대조 · 독성가중치 Table S2 분리 여부
