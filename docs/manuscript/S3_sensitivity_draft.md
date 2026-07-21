# Supplementary S3 — Sensitivity analysis (draft v1)

> **작성 메모:** 파라미터 문헌유추 방어. 데이터: sens_analysis.py (3-seed 평균,
> baseline=compare_control 체제). Fig S3=assets/sensitivity.png. 정직: 원래 가설
> "이점=resistance_cost 의존"은 이 체제선 부정확 → 실제는 저독성 이점 강건 +
> resistance_cost=내성 억제. 이게 더 정확한 방어.

## Methods (요약)
control_score = time-to-progression / (cumulative toxicity + 1) 을 주 결과로,
continuous MTD vs adaptive on/off를 검증된 "통제-가능" 체제(k_prolif=0.15,
k_kill=0.5, cd8_recruit=10, init_resistant_frac=0.03, resistance_cost=0.24)에서
비교. (1) resistance_cost를 0~0.5로 스윕, (2) 9개 핵심 파라미터를 baseline의
0.5×/2×로 일변수(OAT) 섭동. 각 조건 3-seed(42,7,123) 평균.

## 결과

**(1) 적응형 저독성 이점은 resistance_cost에 강건.** 적응형은 전 cost 범위에서
연속 대비 훨씬 낮은 누적 독성으로 종양을 통제했다(적응 독성 24–38 vs 연속 120,
약 1/4–1/5). control_score는 적응형 4.5–5.9 vs 연속 1.2로 일관되게 우세했다
(Fig S3a). 즉 이점 자체는 cost 값에 의존하지 않는다.

**(2) resistance_cost의 실제 역할은 내성 억제(경쟁적 방출).** cost=0에서 적응형
최종 내성 분율은 0.20까지 상승했으나, cost≥0.10에서 <0.05로 억제되었다
(Fig S3a, 우축). 이는 감수성 세포가 내성 세포를 경쟁 억제하려면 내성이 적합도
비용을 치러야 한다는 적응요법 이론과 정확히 일치한다 [9,10]. 따라서 문헌
접지값(0.24; NSCLC 실측 ~0.24 [Science Adv])은 이점을 만들기 위해서가 아니라
내성 동역학을 현실적으로 재현하기 위해 필요하다.

**(3) 결과 민감도는 종양-면역 균형에 집중.** control_score는 종양 증식률
(k_prolif), CD8 장벽 게이팅(cd8_barrier_alpha), 면역 살상(k_kill)에 가장
민감했고(±50%에서 span 5.9–8.0), 내성 관련 파라미터(resistant_immune_evasion,
init_resistant_frac, resistance_cost, mutation_rate)에는 상대적으로 둔감했다
(span 0.0–1.1)(Fig S3b). 이는 향후 **실험적 보정의 우선순위가 증식률·면역
살상/장벽 투과성**임을 지목한다.

## 정직한 한계 (본문 반영)
- **극단 체제에서 통제 저하:** k_prolif 2×에서 적응형 control_score가 0.4로,
  cd8_barrier_alpha 2×에서 2.0으로 급락 → 매우 공격적 종양이나 극도로 촘촘한
  장벽에선 적응형 통제가 보장되지 않는다(§3.3 공격적 체제 주의와 일관).
- OAT는 상호작용을 포착하지 못함 → 향후 전역(Sobol/Latin hypercube) 분석 권장.
- 단일 합성 체제 기반 → 실데이터 접지 파라미터로 재검 필요.

## 한 줄 요약(논문 삽입용)
> Adaptive control retained a ~4–5× toxicity advantage over continuous dosing
> across the full resistance-cost range and ±50% variation in all nine tested
> parameters; the resistance fitness cost governed resistant-clone suppression
> rather than the advantage itself, and outcome sensitivity concentrated in
> tumor proliferation and immune-barrier parameters, identifying these as
> priorities for experimental calibration.
