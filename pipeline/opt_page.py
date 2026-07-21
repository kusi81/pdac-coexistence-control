"""
🎯 통제 최적화 페이지 — 용량 × band 동시 최적화를 대시보드에서.

저독성 자연물 배합의 (용량, 적응형 휴약 band) 격자를 돌려, '통제 유지 최소 부담'
조합을 찾는다. 무거운 계산이라 작은 기본 격자 + st.cache_data로 반응성 유지.
"""
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

from synthetic import make_tissue
from abm import simulate, control_metrics, context_params, SUBSTANCES

NATURAL = ["curcumin", "garlic", "mugwort", "wild_ginseng", "ginsenoside_rg3",
           "platycodon", "sea_cucumber", "danshen", "milk_thistle", "astragaloside"]
REGIMES = {
    "보통 (약물로 통제 가능)": dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5,
                              k_caf_activate=0.10),
    "공격적 (내성 위협·escape)": dict(k_prolif=0.20, cd8_recruit=8, k_kill=0.45,
                                k_caf_activate=0.10),
}
CONTROL_MAX = 1.3


@st.cache_data(show_spinner=False)
def _grid(combo, regime_name, dose_grid, off_grid, adapt_on, days, n_seeds,
          context, res_cost):
    seeds = [42, 7, 123][:n_seeds]
    base = dict(context_params(context), **REGIMES[regime_name],
                resistance_cost=res_cost)
    shape = (len(off_grid), len(dose_grid))
    tox = np.zeros(shape); fin = np.zeros(shape); ctrl = np.zeros(shape)
    for sd in seeds:
        c, l, _ = make_tissue("contained", seed=sd)
        P = dict(base, seed=sd)
        for jd, dose in enumerate(dose_grid):
            regimen = [(s, dose) for s in combo]
            for io, off in enumerate(off_grid):
                h, _ = simulate(c, l, days=days, params=P, regimen_subs=regimen,
                                adaptive=True, adapt_on=adapt_on, adapt_off=off,
                                synergy=0.3, snapshots=(0.0, 1.0))
                m = control_metrics(h, n0=h[0]["n_tumor"])
                tox[io, jd] += m["cum_toxicity"]; fin[io, jd] += m["final_frac"]
                if m["progression_censored"] and m["final_frac"] < CONTROL_MAX:
                    ctrl[io, jd] += 1
    ns = len(seeds)
    return tox / ns, fin / ns, ctrl >= (ns / 2.0)


def _fig(tox, fin, controlled, best, dose_grid, off_grid):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    def heat(ax, M, title, cmap, fmt):
        im = ax.imshow(M, origin="lower", aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(dose_grid)))
        ax.set_xticklabels([f"{d:.0%}" for d in dose_grid])
        ax.set_yticks(range(len(off_grid))); ax.set_yticklabels(off_grid)
        ax.set_xlabel("용량 (전량 대비)"); ax.set_ylabel("adapt_off (휴약 배수)")
        ax.set_title(title, fontsize=10.5, fontweight="bold")
        for io in range(len(off_grid)):
            for jd in range(len(dose_grid)):
                ax.text(jd, io, fmt.format(M[io, jd]), ha="center", va="center",
                        fontsize=8.5, color="black")
                if controlled[io, jd]:
                    ax.add_patch(plt.Rectangle((jd - .5, io - .5), 1, 1, fill=False,
                                               edgecolor="#27AE60", lw=2.4))
        fig.colorbar(im, ax=ax, fraction=0.046)

    heat(a1, tox, "누적 독성(부담) — 초록=통제됨\n낮을수록 좋음", "YlOrRd", "{:.0f}")
    heat(a2, fin, "최종 종양(초기 대비)\n통제 여부", "RdYlGn_r", "{:.1f}")
    if best:
        jb = dose_grid.index(best[1]); ib = off_grid.index(best[2])
        for ax in (a1, a2):
            ax.scatter([jb], [ib], marker="*", s=320, color="#2E86AB",
                       edgecolor="white", linewidth=1.2, zorder=5)
    fig.tight_layout()
    return fig


def render():
    st.sidebar.markdown("---")
    st.sidebar.markdown("**최적화 대상**")
    context = st.sidebar.selectbox("장기 맥락", ["pdac", "hcc"],
                                   format_func=lambda c: {"pdac": "췌장암",
                                                          "hcc": "간세포암"}[c])
    combo = st.sidebar.multiselect(
        "자연물 배합", NATURAL, default=["garlic", "mugwort"],
        format_func=lambda k: SUBSTANCES[k]["label"])
    regime = st.sidebar.selectbox("종양 체제", list(REGIMES.keys()), index=1)
    st.sidebar.markdown("**격자·설정**")
    adapt_on = st.sidebar.slider("투여 시작 band (on)", 1.05, 1.4, 1.15, 0.05)
    days = st.sidebar.slider("기간 (days)", 60, 200, 120, 10)
    n_seeds = st.sidebar.select_slider("seed 수 (평균)", [1, 2, 3], value=1,
                                       help="많을수록 매끄럽지만 느림")
    res_cost = st.sidebar.slider("내성 적합도 비용", 0.0, 0.5, 0.24, 0.02)
    fine = st.sidebar.checkbox("촘촘한 격자(느림)", value=False)
    go = st.sidebar.button("▶ 최적화 실행", type="primary", use_container_width=True)

    st.title("🎯 통제 최적화 — 용량 × band")
    st.markdown(
        "저독성 자연물 배합의 **(용량, 적응형 휴약 band)** 격자를 돌려 "
        "*통제를 유지하면서 부담(독성)이 최소*인 조합을 찾습니다. "
        "'가장 세게'가 아니라 '충분하되 최소로'를 정량으로.")
    st.warning("⚠️ 가설 샌드박스. 파라미터는 문헌 접지(내성)+기전 방향(물질)이며 사람 효과 아님. "
               "격자×seed 계산이라 수십 초 걸릴 수 있습니다(결과는 캐시됨).")

    if not go:
        st.info("왼쪽에서 배합·체제를 정하고 **▶ 최적화 실행**을 누르세요.")
        return
    if not combo:
        st.error("자연물을 하나 이상 선택하세요.")
        return

    dose_grid = ([0.4, 0.6, 0.8, 1.0] if fine else [0.5, 0.75, 1.0])
    off_grid = ([0.4, 0.5, 0.6, 0.7, 0.85] if fine else [0.4, 0.6, 0.8])
    n_cells = len(dose_grid) * len(off_grid) * n_seeds
    with st.spinner(f"{n_cells} 시뮬레이션 실행 중 (용량 {len(dose_grid)} × band "
                    f"{len(off_grid)} × seed {n_seeds})…"):
        tox, fin, controlled = _grid(tuple(combo), regime, tuple(dose_grid),
                                     tuple(off_grid), adapt_on, days, n_seeds,
                                     context, res_cost)

    # 최적: 통제됨 중 독성 최소
    best = None
    for jd, dose in enumerate(dose_grid):
        for io, off in enumerate(off_grid):
            if controlled[io, jd] and (best is None or tox[io, jd] < best[0]):
                best = (tox[io, jd], dose, off, fin[io, jd])

    if best:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("최적 용량", f"{best[1]:.0%}")
        c2.metric("최적 휴약 band", f"{best[2]:.2f}")
        c3.metric("누적 독성(부담)", f"{best[0]:.0f}")
        c4.metric("최종 종양", f"{best[3]:.2f}x")
        # 전량 참고
        if 1.0 in dose_grid:
            jd = dose_grid.index(1.0)
            io_full = min(range(len(off_grid)),
                          key=lambda i: tox[i, jd] if controlled[i, jd] else 1e9)
            if controlled[io_full, jd]:
                full_tox = tox[io_full, jd]
                if best[0] < full_tox:
                    cut = 100 * (1 - best[0] / full_tox)
                    st.success(f"✅ 최적(용량 {best[1]:.0%})이 전량 대비 부담을 "
                               f"**{cut:.0f}% 낮췄습니다** (독성 {best[0]:.0f} vs 전량 "
                               f"{full_tox:.0f}). 통제엔 전량이 불필요.")
    else:
        st.error("이 격자에선 통제되는 조합이 없습니다 — 용량을 높이거나 배합을 강화하세요 "
                 "(약물 바닥 미달).")

    st.pyplot(_fig(tox, fin, controlled, best, dose_grid, off_grid))
    st.caption("★=최적(통제 유지 최소 부담). 초록 테두리=통제됨. 통제 안 되는 낮은 용량 "
               "열은 '최소 유효 용량 바닥' 미달입니다. band와 용량은 부분적으로 서로 "
               "대체되지만 용량 바닥 아래로는 통제 불가.")
