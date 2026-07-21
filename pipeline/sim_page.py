"""
교란 시뮬레이션 페이지 (경로 B) — 물질 조합을 넣고 종양 반응을 시간에 따라 시뮬레이션.

app.py에서 모드가 '교란 시뮬레이션'일 때 render()가 호출된다.
"""

import numpy as np
import pandas as pd
import altair as alt
import streamlit as st

from synthetic import make_tissue
from abm import (simulate, compose_regimen, SUBSTANCES, DEFAULT_PARAMS,
                 EVIDENCE_BADGE, CONTEXTS, context_params,
                 classify_substances, build_cycle_schedule, control_metrics)
import abm_plots


def _line_chart(long_df, ycol, color_domain, color_range, ytitle, height=320):
    """반응형 Altair 라인차트 (컨테이너 폭에 맞춰 동적 스케일)."""
    return alt.Chart(long_df).mark_line(point=True).encode(
        x=alt.X("t:Q", title="시간 (days)"),
        y=alt.Y(f"{ycol}:Q", title=ytitle,
                scale=alt.Scale(zero=False, nice=True)),
        color=alt.Color("arm:N", title="",
                        scale=alt.Scale(domain=color_domain, range=color_range)),
        tooltip=["t:Q", "arm:N", alt.Tooltip(f"{ycol}:Q", format=".0f")],
    ).properties(height=height).interactive()


def render():
    st.sidebar.markdown("---")
    st.sidebar.markdown("**장기 맥락**")
    context = st.sidebar.selectbox(
        "context", list(CONTEXTS.keys()),
        format_func=lambda c: CONTEXTS[c]["label"], label_visibility="collapsed")
    ctx = CONTEXTS[context]

    st.sidebar.markdown("**초기 조직**")
    tissue = st.sidebar.selectbox(
        "합성 조직", ["contained", "diffuse"],
        format_func=lambda m: {"contained": "contained — 섬유아세포 링 + CD8 배제",
                               "diffuse": "diffuse — 구조 없음"}[m])
    seed = st.sidebar.number_input("seed", value=42, step=1)
    days = st.sidebar.slider("시뮬레이션 기간 (days)", 8, 200, 20, 2,
                             help="적응형/저항성 통제는 긴 기간(100~200일) 권장 — "
                                  "경쟁적 방출·공존 역학이 그때 드러납니다.")

    st.sidebar.markdown("**물질 조합 (한약 성분 등)**")
    _ev_mark = {"strong": "🟢", "moderate": "🟡", "weak": "🟠",
                "reference": "⚪", "hypothesis": "🔴"}
    picked = st.sidebar.multiselect(
        "물질 선택", list(SUBSTANCES.keys()),
        format_func=lambda k: f"{_ev_mark.get(SUBSTANCES[k].get('evidence',''),'')} "
                              f"{SUBSTANCES[k]['label']}",
        default=["curcumin"])
    doses = {}
    for name in picked:
        doses[name] = st.sidebar.slider(
            f"용량 — {SUBSTANCES[name]['label']}", 0, 100, 100, 10,
            key=f"dose_{name}") / 100.0
    synergy = st.sidebar.slider("조합 시너지", 0.0, 1.0, 0.0, 0.1,
                                help="조합 시 증식·CAF 억제 방향의 약한 상호작용")

    st.sidebar.markdown("**투여 방식**")
    dosing = st.sidebar.radio(
        "dosing", ["동시 투여", "순차 사이클 (시차)", "적응형 통제 (공존)"],
        label_visibility="collapsed",
        help="순차: 약화→증식억제→기질리셋 반복. 적응형: 종양 부담 band로 투여 on/off "
             "(Gatenby) — 박멸이 아니라 최소 부담으로 오래 통제(공존).")
    phase_days = 5
    if dosing.startswith("순차"):
        phase_days = st.sidebar.slider("phase 길이 (days)", 2, 10, 5, 1)
    adaptive = dosing.startswith("적응형")
    adapt_on, adapt_off = 1.1, 0.7
    if adaptive:
        adapt_on = st.sidebar.slider("투여 시작 band (초기 대비 배수)", 1.05, 1.6, 1.1, 0.05,
                                     help="종양이 초기의 이 배수를 넘으면 투여 시작")
        adapt_off = st.sidebar.slider("휴약 band (초기 대비 배수)", 0.3, 0.95, 0.7, 0.05,
                                      help="종양이 이 배수 밑으로 떨어지면 휴약")

    # 저항성(내성) 역학 — adaptive therapy의 근거. 문헌 접지 기본값.
    st.sidebar.markdown("**저항성 모델** (adaptive therapy)")
    resist_on = st.sidebar.checkbox(
        "내성 아형 활성화", value=adaptive,
        help="내성 종양(약물 회피) + 적합도 비용(내성은 약물 없을 때 느리게 증식). "
             "경쟁적 방출·adaptive therapy를 재현. 문헌 접지값 사용.")
    res_cost = 0.24
    if resist_on:
        res_cost = st.sidebar.slider("내성 적합도 비용", 0.0, 0.5, 0.24, 0.02,
                                     help="내성이 약물 없을 때 감수성의 (1-cost)배 속도로 증식. "
                                          "NSCLC 실측 ~0.24 (Science Adv). 0=비용 없음.")

    go = st.sidebar.button("▶ 시뮬레이션 실행", type="primary",
                           use_container_width=True)

    # ---------- header ----------
    st.title("💊 교란 시뮬레이션 (In-silico)")
    st.markdown(
        "물질(한약 성분·화학요법)을 **파라미터 교란**으로 인코딩해 종양-기질-면역 "
        "동역학을 시간에 따라 돌립니다. 물질이 **섬유화 장벽을 풀어 면역 접근을 여는지**, "
        "**증식을 직접 막는지**를 종양 궤적으로 봅니다.")

    _cbadge = "#C0392B" if context == "hcc" else "#2E86AB"
    st.markdown(
        f"<span style='background:{_cbadge}22;border-left:4px solid {_cbadge};"
        f"padding:.25rem .6rem;border-radius:4px'><b>맥락: {ctx['label']}</b> · "
        f"기질세포 = {ctx['display']['stellate']} → {ctx['display']['barrier']}</span>",
        unsafe_allow_html=True)
    st.caption(("🔬 " + ctx["note"]) +
               ("  이 맥락에서는 **섬유화가 종양 증식을 촉진**(포섭된 촉진자)하므로 "
                "항섬유화 약물이 특히 유효합니다." if context == "hcc"
                else "  이 맥락에서는 섬유화가 **장벽(면역배제)** 로 작동합니다."))

    st.warning(
        "⚠️ **이것은 가설 탐색 샌드박스이지 실제 반응 예측이 아닙니다.** 파라미터는 "
        "문헌에서 기전의 방향·자릿수만 접지했을 뿐 특정 환자/약물에 피팅되지 않았습니다. "
        "물질→배수 매핑은 약리 문헌 기반 '가설'이며, 실측 용량-반응으로 검증해야 합니다.")

    with st.expander("물질별 기전 → 파라미터 매핑 (문헌 근거)", expanded=not go):
        rows = []
        for k, v in SUBSTANCES.items():
            eff = ", ".join(f"{p}×{m}" for p, m in v["effects"].items())
            rows.append({"물질": v["label"],
                         "근거 수준": EVIDENCE_BADGE.get(v.get("evidence", ""), ""),
                         "파라미터 효과": eff, "기전": v["rationale"]})
        st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
        st.caption("k_prolif=증식, k_caf_activate=장벽형성, k_kill=면역살상, "
                   "cd8_recruit=면역유입, k_tumor_apoptosis=세포사멸. "
                   "배수 1.0=효과없음, <1=억제, >1=증가.")
        st.warning("⚠️ **근거 수준 주의:** 녹용은 직접 항암 근거가 약하고(보조적), 대형 "
                   "배당체(도라지·해삼·산삼)는 대부분 in vitro/동물 연구입니다. 파라미터는 "
                   "기전의 방향만 접지한 **가설**이며 사람 대상 효과를 뜻하지 않습니다.")

    if not go:
        st.info("왼쪽에서 물질과 용량을 정한 뒤 **▶ 시뮬레이션 실행**을 누르세요. "
                "대조(무처치)와 처치군이 함께 계산됩니다.")
        return

    if not picked:
        st.error("물질을 하나 이상 선택하세요.")
        return

    # ---------- run ----------
    coords, labels, _ = make_tissue(mode=tissue, seed=int(seed))
    params = context_params(context)   # 맥락 프리셋(PDAC/HCC)을 대조·처치 모두에 적용
    # 저항성 모델: 켜면 문헌접지값+사용자 비용, 끄면 내성 발생 차단(구 동작)
    if resist_on:
        params = dict(params, resistance_cost=res_cost)
    else:
        params = dict(params, init_resistant_frac=0.0, mutation_rate=0.0)
    picked_subs = [(n, doses[n]) for n in picked]
    treated_label = " + ".join(SUBSTANCES[n]["label"].split(" ")[0] for n in picked)
    sequential = dosing.startswith("순차")

    # 순차 사이클 스케줄 구성 (약화 → 증식억제 → 기질리셋 반복)
    schedule = None
    anti, cyto = classify_substances(picked_subs)
    if sequential:
        schedule = build_cycle_schedule(anti, cyto, phase_days=phase_days)

    with st.spinner(f"{days}일 시뮬레이션 실행 중 (대조 + 처치)…"):
        snap_fracs = (0.0, 0.5, 1.0)
        hist_ctrl, _ = simulate(coords, labels, days=days, params=params,
                                snapshots=(0.0, 1.0))
        if adaptive:
            treated_label += " (적응형)"
            hist_tr, snaps_tr = simulate(coords, labels, days=days, params=params,
                                         regimen_subs=picked_subs, synergy=synergy,
                                         adaptive=True, adapt_on=adapt_on,
                                         adapt_off=adapt_off, snapshots=snap_fracs)
            # 비교용 연속(관행) 투여
            hist_simul, _ = simulate(coords, labels, days=days, params=params,
                                     regimen_subs=picked_subs, synergy=synergy,
                                     snapshots=(0.0, 1.0))
        elif sequential and schedule:
            hist_tr, snaps_tr = simulate(coords, labels, days=days, params=params,
                                         schedule=schedule, synergy=synergy,
                                         snapshots=snap_fracs)
            hist_simul, _ = simulate(coords, labels, days=days, params=params,
                                     regimen_subs=picked_subs, synergy=synergy,
                                     snapshots=(0.0, 1.0))
        else:
            hist_tr, snaps_tr = simulate(coords, labels, days=days, params=params,
                                         regimen_subs=picked_subs, synergy=synergy,
                                         snapshots=snap_fracs)
            hist_simul = None

    # 순차 사이클 구성 표시
    if sequential:
        if not schedule:
            st.error("순차 사이클을 만들 물질이 부족합니다. 물질을 더 선택하세요.")
            return
        _lab = lambda subs: ", ".join(SUBSTANCES[n]["label"].split(" ")[0]
                                      for n, _ in subs) or "—"
        phase_txt = "  →  ".join(f"**{ph['label']}**({ph['days']}일: {_lab(ph['substances'])})"
                                 for ph in schedule)
        st.info(f"🔄 순차 사이클 ({sum(ph['days'] for ph in schedule)}일/사이클, {days}일 반복): "
                + phase_txt)

    n0 = hist_ctrl[0]["n_tumor"]
    nc = hist_ctrl[-1]["n_tumor"]
    nt = hist_tr[-1]["n_tumor"]
    pct_ctrl = 100.0 * (nc - n0) / n0
    pct_tr = 100.0 * (nt - n0) / n0
    # 처치가 대조 대비 줄인 종양 (성장 억제율)
    inhibition = 100.0 * (nc - nt) / nc if nc else 0.0

    tox_tr = hist_tr[-1]["cum_toxicity"]
    efficiency = inhibition / tox_tr if tox_tr > 0 else 0.0

    # ---------- summary cards ----------
    st.subheader("결과 요약")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("대조 대비 억제율", f"{inhibition:.0f}%",
              help="(대조최종−처치최종)/대조최종. 양수 = 처치가 종양을 줄임.")
    c2.metric("처치 종양 (최종)", f"{nt:,}", f"{pct_tr:+.0f}%", delta_color="inverse")
    c3.metric("누적 독성(부담)", f"{tox_tr:.1f}",
              help="활성 물질 독성가중치×기간 적산. 낮을수록 환자 부담↓. 예시값이며 임상수치 아님.")
    c4.metric("효율 (억제%/독성)", f"{efficiency:.1f}",
              help="같은 억제를 더 낮은 독성으로 달성할수록 높음.")

    # ---------- 통제·공존 지표 (적응형/저항성) ----------
    if adaptive or resist_on:
        st.markdown("#### 🧬 통제·공존 지표 (adaptive therapy)")
        m_tr = control_metrics(hist_tr, n0=n0)
        d1, d2, d3, d4 = st.columns(4)
        d1.metric("진행까지 시간 TTP",
                  f"{m_tr['ttp_days']:.0f}d" + ("*" if m_tr['progression_censored'] else ""),
                  help="종양이 초기 1.5배 도달 시각. *=관측기간 내 진행 없음(통제 유지).")
        d2.metric("내성 분율(최종)", f"{m_tr['final_resistant_frac']:.2f}",
                  help="종료 시 약물내성 종양 비율. 낮을수록 통제 양호(경쟁적 방출 억제).")
        d3.metric("종양(초기 대비)", f"{m_tr['final_frac']:.2f}x")
        d4.metric("통제점수 TTP/부담", f"{m_tr['control_score']:.1f}",
                  help="적은 부담으로 오래 통제할수록 높음.")

        # 적응형 vs 연속(관행) 비교
        if adaptive and hist_simul is not None:
            m_co = control_metrics(hist_simul, n0=n0)
            cmp = pd.DataFrame([
                {"전략": "연속(관행)", "TTP": f"{m_co['ttp_days']:.0f}d",
                 "최종종양": f"{m_co['final_frac']:.2f}x",
                 "내성": f"{m_co['final_resistant_frac']:.2f}",
                 "누적독성": f"{m_co['cum_toxicity']:.0f}"},
                {"전략": "적응형(공존)", "TTP": f"{m_tr['ttp_days']:.0f}d",
                 "최종종양": f"{m_tr['final_frac']:.2f}x",
                 "내성": f"{m_tr['final_resistant_frac']:.2f}",
                 "누적독성": f"{m_tr['cum_toxicity']:.0f}"},
            ])
            st.dataframe(cmp, hide_index=True, use_container_width=True)
            if m_tr["cum_toxicity"] < m_co["cum_toxicity"] * 0.85:
                cut = 100 * (1 - m_tr["cum_toxicity"] / max(m_co["cum_toxicity"], 1))
                st.success(f"✅ 적응형이 종양을 통제하면서 **부담을 {cut:.0f}% 낮췄습니다** "
                           f"(독성 {m_tr['cum_toxicity']:.0f} vs 연속 {m_co['cum_toxicity']:.0f}). "
                           "'박멸이 아니라 최소 부담으로 공존' — Gatenby 적응요법.")

        # 내성 분율 궤적 (경쟁적 방출)
        rrows = [{"t": h["t"], "arm": treated_label, "resistant_frac": h["resistant_frac"]}
                 for h in hist_tr]
        rdom, rcol = [treated_label], ["#C0392B"]
        if adaptive and hist_simul is not None:
            rrows += [{"t": h["t"], "arm": "연속(관행)", "resistant_frac": h["resistant_frac"]}
                      for h in hist_simul]
            rdom, rcol = [treated_label, "연속(관행)"], ["#27AE60", "#C0392B"]
        st.altair_chart(
            _line_chart(pd.DataFrame(rrows), "resistant_frac", rdom, rcol, "내성 분율"),
            use_container_width=True)
        st.caption("내성 분율이 치솟으면 '경쟁적 방출'(감수성이 사라져 내성이 장악). 적응형은 "
                   "감수성을 남겨 내성을 억누릅니다. 내성 적합도 비용은 문헌 접지(NSCLC ~0.24).")

    # 순차 vs 동시 비교 (사용자 핵심 질문)
    if sequential and hist_simul is not None:
        ns = hist_simul[-1]["n_tumor"]
        inh_s = 100.0 * (nc - ns) / nc if nc else 0.0
        tox_s = hist_simul[-1]["cum_toxicity"]
        eff_s = inh_s / tox_s if tox_s > 0 else 0.0
        st.markdown("#### 순차 사이클 vs 동시 투여")
        cmp = pd.DataFrame([
            {"방식": "동시 투여", "억제율": f"{inh_s:.0f}%",
             "누적 독성": f"{tox_s:.1f}", "효율": f"{eff_s:.1f}"},
            {"방식": "순차 사이클", "억제율": f"{inhibition:.0f}%",
             "누적 독성": f"{tox_tr:.1f}", "효율": f"{efficiency:.1f}"},
        ])
        st.dataframe(cmp, hide_index=True, use_container_width=True)
        if efficiency > eff_s * 1.1 and inhibition > inh_s - 10:
            st.success(f"✅ 순차 사이클이 비슷한 억제({inhibition:.0f}% vs {inh_s:.0f}%)를 "
                       f"**더 낮은 독성**({tox_tr:.1f} vs {tox_s:.1f})으로 달성 → 효율 우위. "
                       "적응요법(adaptive therapy)의 논리와 부합합니다.")
        elif inh_s > inhibition + 10:
            st.info(f"이 설정에선 동시 투여가 억제는 더 크나({inh_s:.0f}% vs {inhibition:.0f}%) "
                    f"독성도 큽니다({tox_s:.1f} vs {tox_tr:.1f}). phase 길이·조합을 바꿔보세요.")
    else:
        if inhibition > 15:
            st.success(f"✅ 처치군이 대조 대비 종양을 **{inhibition:.0f}%** 줄였습니다.")
        elif inhibition > 3:
            st.info(f"처치군이 대조 대비 종양을 {inhibition:.0f}% 줄였습니다 (부분 효과).")
        else:
            st.warning("이 조합/용량으로는 뚜렷한 억제가 관찰되지 않았습니다.")

    # 대조/처치 궤적을 long-form으로 (Altair 반응형 차트용)
    def _long(key):
        rows = []
        for h in hist_ctrl:
            rows.append({"t": h["t"], "arm": "대조", key: h[key]})
        for h in hist_tr:
            rows.append({"t": h["t"], "arm": treated_label, key: h[key]})
        return pd.DataFrame(rows)

    dom = ["대조", treated_label]

    # ---------- trajectories ----------
    st.subheader("① 종양 성장 궤적")
    tcol1, tcol2 = st.columns(2)
    with tcol1:
        st.markdown("**종양 세포 수**")
        st.altair_chart(
            _line_chart(_long("n_tumor"), "n_tumor", dom, ["#7F8C8D", "#C0392B"],
                        "종양 세포 수"),
            use_container_width=True)
    with tcol2:
        # 누적 독성(부담) 궤적 — 순차면 동시 투여와 비교
        st.markdown("**누적 독성(환자 부담)**")
        trows = [{"t": h["t"], "arm": treated_label, "cum_toxicity": h["cum_toxicity"]}
                 for h in hist_tr]
        tdom, tcolors = [treated_label], ["#C0392B"]
        if sequential and hist_simul is not None:
            trows += [{"t": h["t"], "arm": "동시 투여", "cum_toxicity": h["cum_toxicity"]}
                      for h in hist_simul]
            tdom, tcolors = [treated_label, "동시 투여"], ["#27AE60", "#C0392B"]
        st.altair_chart(
            _line_chart(pd.DataFrame(trows), "cum_toxicity", tdom, tcolors,
                        "누적 독성"),
            use_container_width=True)
    st.caption("누적 독성 = 활성 물질 독성가중치 × 투여기간 적산 (예시값, 임상수치 아님). "
               "순차 사이클은 강한 약을 쉬어가며 돌려 **누적 부담을 낮춥니다**.")

    st.subheader("② 기전 판독 — 면역 침투 & 장벽 질량")
    m1, m2 = st.columns(2)
    with m1:
        st.markdown("**CD8 침투도** (종양 25µm 내 비율)")
        st.altair_chart(
            _line_chart(_long("cd8_infiltration"), "cd8_infiltration", dom,
                        ["#7F8C8D", "#27AE60"], "침투도"),
            use_container_width=True)
    with m2:
        st.markdown("**myCAF 장벽 세포 수**")
        st.altair_chart(
            _line_chart(_long("n_myCAF"), "n_myCAF", dom,
                        ["#7F8C8D", "#C0392B"], "myCAF 수"),
            use_container_width=True)
    st.caption("항섬유화 계열은 myCAF 장벽을 낮춰 CD8 침투를 높입니다(왼쪽↑, 오른쪽↓). "
               "세포독성 계열은 장벽과 무관하게 증식을 직접 억제합니다.")

    # ---------- snapshots (컬럼 분할로 한 화면에 반응형 배치) ----------
    st.subheader("③ 조직 진화 (처치군)")
    times = sorted(snaps_tr.keys())
    scols = st.columns(len(times))
    for col, tt in zip(scols, times):
        coords_s, labels_s = snaps_tr[tt]
        with col:
            st.pyplot(abm_plots.fig_snapshot_single(coords_s, labels_s, tt),
                      use_container_width=True)

    # ---------- downloadable trajectory ----------
    df = pd.DataFrame(hist_tr).add_prefix("treated_").join(
        pd.DataFrame(hist_ctrl).add_prefix("control_"))
    st.download_button("⬇ 궤적 CSV 다운로드", df.to_csv(index=False),
                       file_name="sim_trajectory.csv", mime="text/csv")
