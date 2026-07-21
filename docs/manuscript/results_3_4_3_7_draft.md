# Results §3.4–3.7 (draft v1)

> **작성 메모:** 수치는 실제 실행값(compare_control.py, optimize_natural.py,
> optimize_dose_band.py; 2026-07-21). 전부 in-silico *예측*으로 명시. artifact 정직
> 고지(무처치 control_score, 공격적 조합의 경쟁적 방출). Fig 4=control_strategies.png,
> Fig 5=natural_adaptive_optim.png, Fig 6=dose_band_optimization.png, Fig 7=분자구조.
> [n]=Introduction 참조번호. 저널 확정 후 인용스타일 조정.

---

## 3.4 Adaptive scheduling achieves tumor coexistence at a fraction of the toxicity of continuous dosing

Using the spatial ABM with an explicit resistant subclone, we compared three
strategies on a controllable-but-not-eradicable tumor: no treatment, continuous
maximum-dose therapy, and adaptive on/off dosing that toggles treatment around a
tumor-burden band (Fig. 4). Continuous dosing drove the sensitive population to
extinction (final burden 0.00× baseline) but required the highest cumulative
toxicity (120 arbitrary units) and, by removing the sensitive competitors, left a
purely resistant residue. Adaptive dosing instead held the tumor at a stable
0.82× baseline—coexistence rather than eradication—at roughly one-fifth the
cumulative toxicity (25 units), while keeping the resistant fraction low (0.03).
Both treated arms remained below the progression threshold over the full 150-day
horizon (censored), whereas the untreated tumor progressed at day 103.

We note that our composite control_score (progression time divided by cumulative
toxicity) assigns an artifactually high value to the untreated arm because its
toxicity denominator is zero; we therefore rank strategies by progression-control
first (censored vs progressed) and only then by toxicity, which excludes the
untreated arm despite its nominal score. Under this ranking, adaptive scheduling
is the preferred strategy: it achieves durable control at a small fraction of the
toxicity burden of continuous dosing (control_score 5.7 vs 1.2).

---

## 3.5 Food-medicine-homology regimens control tumor burden at low predicted toxicity, with anti-CAF agents minimizing resistance

We next encoded food-medicine-homology compounds as mechanistic perturbations and
evaluated single agents and combinations under adaptive scheduling, benchmarked
against conventional continuous gemcitabine (Fig. 5; all values in-silico
predictions). Every natural regimen tested maintained the tumor below the
progression threshold for the full horizon, but at markedly lower predicted
toxicity than gemcitabine (5–38 units vs 128). The most favorable
control-per-toxicity profiles were single anti-CAF/immunomodulatory agents—garlic
(toxicity 5, final 1.0×, resistant 0.05, control_score 25.0) and wild ginseng
(toxicity 5, final 0.9×, control_score 24.0)—which stabilized the tumor near
baseline at roughly one-twenty-fifth the toxicity of gemcitabine.

Two patterns are worth emphasizing. First, deeper tumor suppression did not equate
to better control. Aggressive natural combinations reduced burden further (e.g.,
garlic + mugwort, 0.4× baseline) but at the cost of a sharply elevated resistant
fraction (0.34), a direct in-silico illustration of competitive release: hard
suppression of sensitive cells selects for resistance. The balanced combination
curcumin + garlic + ginsenoside-Rg3 achieved moderate suppression (0.8×) while
holding resistance lowest (0.01) at modest toxicity (11), and was carried forward
for dose optimization (§3.6). Second, anti-fibrotic activity alone was
insufficient: a purely anti-fibrotic pairing (danshen + astragaloside) barely held
the tumor (1.2× baseline) at the highest natural-arm toxicity (38), indicating
that stromal modulation requires an anti-proliferative or immunomodulatory partner
to constrain the tumor. Together these results support a control-oriented design
principle—pair a low-toxicity anti-CAF backbone with a resistance-sparing
anti-proliferative component rather than maximizing cytotoxic pressure.

---

## 3.6 Dose and drug-holiday optimization further lowers predicted toxicity

For the lead balanced combination (curcumin + garlic + ginsenoside-Rg3), we
performed a seed-averaged grid search over dose intensity and off-threshold (band
lower bound) under adaptive scheduling over 200 days (Fig. 6). The optimum reduced
dose to 80% with an off-threshold of 0.5, achieving control at a cumulative
toxicity of 27 with a resistant fraction of 0.02 and final burden 0.8× baseline.
This modestly but consistently outperformed full-dose adaptive scheduling
(toxicity 32 at comparable control, 0.7× burden), indicating that a sub-maximal
dose with appropriately timed holidays preserves control while further trimming
toxicity. The result reinforces the central theme: within a control (rather than
eradication) objective, less drug—delivered adaptively—can be better.

---

## 3.7 Multiscale molecular grounding links compounds to model parameters under a transparent evidence hierarchy

To connect the cellular-scale perturbations to their molecular rationale, each
compound is annotated with its principal target(s) and an explicit evidence tier:
experimental co-crystal, computational docking prediction, or pathway-level
mechanistic inference without a solved structure (Fig. 7). This tiering is
deliberately conservative about what a 3D structure can claim. Only the two
conventional agents are supported by experimental co-crystals in which the
compound itself is resolved in the binding site—gemcitabine bound to
deoxycytidine kinase (PDB 1P62) and erlotinib bound to the EGFR kinase domain
(PDB 1M17). A minority of interactions rest on docking predictions (e.g.,
curcumin–STAT3, PDB 1BG1; the entecavir–KDM5B repurposing hypothesis, PDB 5A1F),
where the displayed ligand is a reference rather than our compound. The majority
of the food-medicine-homology anti-fibrotic compounds act through the TGF-β/Smad
and NF-κB axes for which no compound-bound structure is available and are
therefore represented at the mechanistic tier only.

Each molecular target maps to the ABM parameter class it perturbs—anti-
proliferative targets to the tumor proliferation rate, anti-CAF/anti-fibrotic
targets (predominantly TGF-β/Smad) to CAF activation, and immunomodulatory targets
to CD8 recruitment and killing—providing the molecule-to-cell-to-tissue link that
the framework operationalizes. We stress that this hierarchy is a transparency
device, not a claim of binding validation: the mechanistic-tier compounds that
dominate our low-toxicity regimens are grounded in pathway pharmacology, and their
molecular engagement remains a hypothesis for experimental confirmation.

---

## 편집 체크리스트
- [ ] Fig 4/5/6/7 파일 확정 + 캡션 (Fig 7=drug_structures_3d.png 재사용 여부)
- [ ] "arbitrary units" 독성 정의 Methods 상호참조
- [ ] control_score artifact 처리(censored-first)를 Methods §2.7과 일치
- [ ] 경쟁적 방출(garlic+mugwort 0.34) → Discussion에서 재언급(통제 vs 박멸)
- [ ] 수치 유효숫자 최종 대조(optimize_natural/dose_band 콘솔값)
- [ ] in-silico/predicted 수식어 전수 확인
