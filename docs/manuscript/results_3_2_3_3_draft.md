# Results §3.2–3.3 (draft v1)

> **작성 메모 (한글):** Introduction v2/Abstract 용어 통일. 수치는 실제 실행값
> (scotia_posctrl.py, scotia_rim.py, diag_mycaf.py). 전부 방향성 + "consistent with"
> hedging. 데이터: SCOTIA=CosMx SMI 717,493셀×1,009유전자 16시료(Untreated 9 / CRT 6 /
> CRTL 1) [20]; Xenium=GSE274673 5시료(자체 module-score 주석). Fig 2=검증, Fig 3=완성 F3.
> 저널 확정 후 인용스타일·figure 번호 조정.

---

## 3.2 The containment metric is validated by an author-annotated positive control, and localizes annotation—not metric—as the limiting factor

We first established that our spatial barrier/containment metrics behave correctly
on ground-truth-like inputs. On synthetic tissues, the barrier score cleanly
separated an architecture in which stroma is interposed between tumor and immune
cells (contained; z ≈ +22.3) from one in which the same cell counts are randomly
intermixed (diffuse; z ≈ +0.9), confirming that the matched-null construction
responds to spatial *geometry* rather than abundance (Fig. 2a).

We next asked whether the metric reproduces a known biological positive control on
real data: myofibroblastic CAFs (myCAF) are reported to lie closer to malignant
cells than inflammatory CAFs (iCAF) [3,20]. Using the author-annotated CosMx SMI
dataset of Shiau et al. [20] (717,493 cells across 16 tumors), we applied our
proximity test (median distance of each CAF subtype to the nearest malignant cell,
with label-permutation nulls) per sample. myCAF was significantly closer to
malignant cells than iCAF in 15 of 16 tumors (the single exception, U2-b, showed a
0.2 µm difference, p = 0.279, i.e., a statistical tie rather than a reversal). The
effect was consistent across treatment groups and, notably, larger after therapy:
the myCAF–iCAF proximity gap widened from −2.6 µm (untreated; myCAF 20.6 µm vs iCAF
23.1 µm) to −18.8 µm (CRT; myCAF 42.1 µm vs iCAF 61.0 µm) (Fig. 2b). Thus, given
faithful cell annotations, our metric recovers the expected myCAF-proximal /
iCAF-distal architecture.

Critically, the *same* metric failed this positive control on our targeted-panel
Xenium data (GSE274673) when cells were typed with marker-based module scores: myCAF
was scored as marginally *farther* from tumor than iCAF (e.g., +1.5 to +6.0 µm;
positive control failed in 4–5 of 5 tumors). We traced this to the annotation, not
the metric. The 480-gene Xenium panel contained only 4 of 8 canonical myCAF
subtype markers, and the two present pan-myofibroblast markers (ACTA2, TAGLN) are
shared with perivascular/smooth-muscle cells. A refined annotation that added the
gold-standard myCAF markers available in the panel (POSTN, LRRC15, CTHRC1) and
introduced a dedicated perivascular class did separate perivascular cells but still
failed the positive control (4 of 5 tumors), indicating that marker-based
winner-take-all typing on a targeted panel cannot resolve the fine myCAF/iCAF
spatial gradient (Fig. 2c). Together these results establish two points that frame
the rest of the study: (i) our containment metrics are sound, as they reproduce the
literature positive control on full-panel author-annotated data; and (ii) robust
CAF-subtype spatial analysis requires full-transcriptome or author-level
annotations rather than targeted-panel module scoring—consistent with our framing
of the framework as hypothesis-generating rather than annotation-confirmatory.

---

## 3.3 myCAFs form a tumor-adjacent barrier that tightens under therapy while cytotoxic T cells remain excluded

Having validated the metric, we characterized which cells occupy the immediate
peritumoral rim (30 µm shell) and how this changes with neoadjuvant
chemoradiotherapy, using the author-annotated CosMx cohort (9 untreated, 6 CRT;
one CRTL sample excluded) (Fig. 3).

myCAFs were the dominant tumor-adjacent population, strongly enriched at the tumor
rim in both groups and further enriched after therapy (mean rim enrichment z = +9.1
untreated → +15.7 CRT), whereas iCAFs were excluded from the rim and displaced
farther out under treatment (z = −3.5 → −9.2). This treatment-associated tightening
of the myCAF barrier and displacement of iCAF mirrors the widening myCAF–iCAF
proximity gap seen in §3.2 and is consistent with CRT-driven remodeling toward a
matrix-rich, myCAF-encased tumor front [20]. Pericytes shifted from neutral to
tumor-adjacent after therapy (z = −1.6 → +3.3), consistent with vascular
remodeling, while endothelial cells remained near-neutral.

The immune compartment showed persistent exclusion. Cytotoxic CD8⁺ T cells were
depleted from the tumor rim in both groups and remained excluded after therapy
(z = −5.6 untreated → −4.0 CRT); CD4⁺ T cells and regulatory T cells were similarly,
if more mildly, excluded (CD4⁺ −2.6 → −1.3; Treg −3.6 → −1.6). Plasma cells were the
most strongly rim-excluded population (z = −11.3 → −4.1), consistent with a
predominantly distal, stromal localization. Macrophages were the only immune-lineage
population enriched at the tumor rim, and this enrichment declined after therapy
(z = +6.8 → +2.4).

These observations directly test, and refine, the intuitive expectation that
therapy relieves immune exclusion. In this cohort, CRT did *not* recruit cytotoxic
T cells to the tumor rim: CD8⁺ T cells remained excluded while the myCAF barrier
grew tighter. Rather than contradicting the containment hypothesis, this pattern
supports it—the myCAF-derived barrier persists and is reinforced under treatment,
maintaining a spatial separation between malignant cells and cytotoxic effectors.
This is precisely the configuration our control framework seeks to exploit: a
durable, treatment-stable myCAF structure that can, in principle, be modulated to
restrain rather than merely attack the tumor (§3.4–§3.6). We note that these are
associations in a modest cohort (n = 15) and that our use of the authors'
annotations makes the myCAF/iCAF positive control confirmatory by construction;
the peritumoral compositional findings for immune and vascular lineages, however,
are independent of that construction.

---

## 편집 체크리스트
- [ ] Fig 2 패널 배치 확정 (a 합성 / b SCOTIA 양성대조 / c Xenium 실패 대조)
- [ ] Fig 3 = rim_scotia.png (글리프 − 수정본)
- [ ] px→µm(0.12028), 30µm shell, permutation n=200–300을 Methods와 상호참조
- [ ] 인용번호 [3][20] Introduction과 일치 확인
- [ ] "CRTL" 정의(장기 CRT?) Methods에 각주
- [ ] 수치 유효숫자·부호 최종 대조 (scotia_posctrl.csv / scotia_rim.csv)
