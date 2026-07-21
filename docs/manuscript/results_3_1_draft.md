# Results §3.1 (draft v1)

> **작성 메모:** 신규성 landscape. 데이터: lit_search.py 체계검색(10 쿼리, 132편),
> Fig 1=fig1_novelty.png, 상세=novelty_assessment.md. 정직: PubMed-only, 표적검색
> (PRISMA 아님). [n]=Introduction 참조. Methods §2.1과 상호참조.

---

## 3.1 The components of our approach are individually well-established, but their integration is unoccupied

To position the framework, we ran a systematic PubMed survey (ten Boolean
queries spanning the intersecting dimensions of our approach; Methods §2.1) and
quantified how much prior literature occupies each axis versus their integration
(Fig. 1a). The individual axes are well populated: coexistence/control with
natural compounds (56 records), network-pharmacology of food-medicine-homology
compounds against cancer (33), and spatial/agent-based modeling of PDAC (31) each
returned substantial literature. In contrast, queries requiring the *simultaneous*
combination of dynamic/computational modeling, natural compounds, and CAF biology
returned only 2–3 records, and on inspection all were off-target—reviews or
wet-laboratory studies in other tumor types (e.g., an ovarian TCM–chemotherapy
review; an NSCLC phytochemical study), none of which couples a spatial/dynamic
model to natural-compound perturbation of PDAC CAF biology.

The relevant prior art therefore resolves into three camps that our survey shows
are individually mature but mutually disconnected (Fig. 1b): (i) static network
pharmacology of food-medicine-homology compounds, which maps compound–target
relationships but is time- and space-agnostic [14,16]; (ii) wet-laboratory studies
of herbal modulation of PDAC CAFs, which provide biology but no predictive or
optimization framework [17,18]; and (iii) agent-based and spatial models of the
PDAC microenvironment, which capture tumor–stroma–immune dynamics but incorporate
neither natural compounds nor a coexistence-control objective [19,20]. The triple
intersection—a dynamic, spatially grounded, adaptive-control framework that uses
food-medicine-homology compounds to modulate the myCAF barrier—was empty in our
survey, defining the gap this work occupies.

We note the limits of this analysis: it queries PubMed only (excluding, e.g.,
bioRxiv/arXiv preprints and some non-English venues) and is a targeted novelty
search rather than a PRISMA systematic review. Nonetheless, the near-absence of
integrated hits, combined with the abundance of single-axis literature, provides
quantitative support that the contribution of this work lies in the integration
rather than in any single component—consistent with our framing of the framework
as connecting, rather than replacing, these established research lines.

---

## 편집 체크리스트
- [ ] Fig 1 = fig1_novelty.png (a 막대 / b 벤)
- [ ] 인용 [14][16][17][18][19][20] Introduction과 일치
- [ ] off-target 예시 PMID 각주(41372788 난소 리뷰, 39851099 NSCLC) 넣을지 결정
- [ ] Methods §2.1(쿼리 전문·PRISMA 향후)과 상호참조
- [ ] "약 7편"(camp ②=D3+D4 합) 계산근거 각주
