# 투고 실행 가이드 — bioRxiv 선공개 → PLOS Computational Biology (waiver)

> 원고 콘텐츠는 완결(`manuscript.md`). 이 문서는 **당신이 실행할 단계**와 준비물 정리.
> 실제 업로드·제출·계정 로그인은 당신이 직접(공개 게시 행위).

## 전략 요약
1. **bioRxiv 프리프린트 선공개** (무료) → 우선권·타임스탬프·협업 지렛대 즉시 확보.
2. **PLOS Computational Biology 정식 투고** (waiver 신청 → 게재비 $0 목표).
- bioRxiv는 PLOS 투고를 **방해하지 않음** — PLOS는 프리프린트를 명시적으로 환영하고,
  bioRxiv에는 저널 직접 전송 기능도 있음.

---

## 준비물 체크 (모두 준비됨)
- [x] 원고 본문 `manuscript.md` (IMRaD + Declarations + 참조 20)
- [x] 그림 8개 `assets/` (fig1_novelty, fig2_validation, rim_scotia, control_strategies,
      natural_adaptive_optim, dose_band_optimization, drug_structures_3d, sensitivity)
- [x] Supp Table S1 `S1_parameters.md`, S3 본문 `S3_sensitivity_draft.md`
- [x] 저자: Seung-Il Kim, ORCID 0009-0007-5965-9212, kusi81kim@gmail.com
- [x] 코드 리포: https://github.com/kusi81/pdac-coexistence-control
- [ ] **제출용 단일 파일**(PDF 또는 Word) — 아래 "파일 변환" 참고
- [ ] Fig 3/5/6 글리프 마감(cosmetic) — 카메라레디 전 권장

---

## STAGE 1 — bioRxiv (당신이 실행)

**1. 계정** — biorxiv.org 가입(무료). ORCID 연동 권장.

**2. 새 제출** → 메타데이터 입력:
- **Title:** A spatially grounded agent-based framework for coexistence-control of
  pancreatic cancer using food-medicine-homology compounds
- **Authors:** Seung-Il Kim (Independent Researcher), ORCID 0009-0007-5965-9212, 교신
- **Abstract:** `manuscript.md`의 Abstract 붙여넣기
- **Category(주제):** **Cancer Biology** (1순위) 또는 Systems Biology / Bioinformatics
- **License:** **CC-BY 4.0** 권장 (PLOS도 CC-BY라 호환; 재사용·인용 최대)
- **Type:** New Results

**3. 파일 업로드:** 제출용 PDF(본문+그림) + 그림 원본(선택) + 보충자료(S1/S3).

**4. 제출 → 스크리닝**(보통 1–2 영업일) → 공개되면 **DOI 부여**. 이 DOI를 협업
제안·이력서·리포 About의 website 칸에 사용.

---

## STAGE 2 — PLOS Computational Biology (bioRxiv 후, 당신이 실행)

**1. 형식 확인** — journals.plos.org/ploscompbiol 제출 가이드.
- Author Summary **필수** (이미 보유 ✅), 비구조 초록(보유 ✅).
- Data Availability / Funding / Competing Interests / Ethics 문장 필요 (이미 Declarations에 ✅).
- 참조 = Vancouver 스타일(현재와 유사) — 제출 시스템이 대부분 처리.

**2. ⭐ Waiver 신청 (핵심):**
- 제출 과정 중 **"Publication fee assistance"** 단계에서 개인 waiver 신청.
- 근거: **무소속·무연구비 개인 연구자**(Independent Researcher, no funding).
- Funding 문장이 "no specific funding"이라 waiver 논리와 일관 ✅.
- 편집부가 사안별 심사 → 승인 시 게재비 $0.

**3. Cover letter** (제출 시) — 초안:
> We submit "…coexistence-control of pancreatic cancer using food-medicine-homology
> compounds," an in-silico, hypothesis-generating framework that bridges static
> network pharmacology, wet-lab CAF biology, and spatial agent-based modeling. A
> systematic survey (Fig 1) shows this integration is unoccupied. We validate our
> containment metric against an author-annotated positive control (15/16, Fig 2) and
> report testable, low-toxicity regimens for tumor coexistence. As an unaffiliated
> researcher without funding, I request consideration for a publication-fee waiver.
> A preprint is available on bioRxiv [DOI].

**4. 제출 후** — 리뷰(수 주~개월). 리뷰 코멘트 대응은 준비된 한계 섹션(§4.5)·민감도
(§S3)로 상당 부분 방어 가능.

---

## 파일 변환 (다음 작업 후보)
`manuscript.md` → 제출용 단일 파일:
- **옵션 A (권장): Word .docx** — 편집·공동저자 추가·저널 시스템 호환 최선.
- **옵션 B: PDF** — bioRxiv 제출에 무난.
→ 원하면 다음 단계로 생성.

## 남은 cosmetic (카메라레디 전)
- Fig 3/5/6 유니코드 마이너스 글리프 정리(CSV replot)
- "arbitrary units" 독성 정의 한 문장(§2.5/2.7)
- CRTL 각주(§2.2)
