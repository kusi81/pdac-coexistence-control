"""
용량 × band 동시 최적화 — 두 통제 레버를 함께 조정.

용량(dose)과 휴약 band(adapt_off) 둘 다 '얼마나 세게 미느냐'를 결정한다:
  - 용량↑ → 회당 효과↑·독성↑
  - off↓ (더 깊이 눌러야 휴약) → 투여 빈도↑·독성↑
질문: 저용량 + 적절한 band로도 통제하면서 부담을 더 낮출 수 있는가?

목적(제약 최적화): 통제됨(진행없음 & 최종<1.3x) 중 누적독성 최소.
공격적 종양 + 문헌 접지 저항성(cost 0.24 등), 3-seed 평균. ⚠️ 가설 샌드박스.
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

from synthetic import make_tissue
from abm import simulate, control_metrics

PARAMS = dict(k_prolif=0.20, cd8_recruit=8, k_kill=0.45, k_caf_activate=0.10,
              init_resistant_frac=0.01, mutation_rate=0.001,
              resistant_immune_evasion=0.45, resistance_cost=0.24)
SUBS = ["curcumin", "garlic", "ginsenoside_rg3"]   # 커큐민+마늘+Rg3
DAYS = 200
ADAPT_ON = 1.15
DOSE_GRID = [0.4, 0.6, 0.8, 1.0]
OFF_GRID = [0.4, 0.5, 0.6, 0.7, 0.85]
CONTROL_MAX = 1.3
SEEDS = [42, 7, 123]


def main():
    tox = np.zeros((len(OFF_GRID), len(DOSE_GRID)))
    resf = np.zeros_like(tox)
    fin = np.zeros_like(tox)
    ctrl_count = np.zeros_like(tox)

    print(f"용량{len(DOSE_GRID)} × off{len(OFF_GRID)} × seed{len(SEEDS)} "
          f"(커큐민+마늘+Rg3, on={ADAPT_ON}, {DAYS}일)...")
    for sd in SEEDS:
        c, l, _ = make_tissue("contained", seed=sd)
        P = dict(PARAMS, seed=sd)
        for jd, dose in enumerate(DOSE_GRID):
            regimen = [(s, dose) for s in SUBS]
            for io, off in enumerate(OFF_GRID):
                h, _ = simulate(c, l, days=DAYS, params=P, regimen_subs=regimen,
                                adaptive=True, adapt_on=ADAPT_ON, adapt_off=off,
                                synergy=0.3, snapshots=(0.0, 1.0))
                m = control_metrics(h, n0=h[0]["n_tumor"])
                tox[io, jd] += m["cum_toxicity"]
                resf[io, jd] += m["final_resistant_frac"]
                fin[io, jd] += m["final_frac"]
                if m["progression_censored"] and m["final_frac"] < CONTROL_MAX:
                    ctrl_count[io, jd] += 1
    ns = len(SEEDS)
    tox /= ns; resf /= ns; fin /= ns
    controlled = ctrl_count >= (ns / 2.0)

    best = None
    for jd, dose in enumerate(DOSE_GRID):
        for io, off in enumerate(OFF_GRID):
            if controlled[io, jd] and (best is None or tox[io, jd] < best[0]):
                best = (tox[io, jd], dose, off, resf[io, jd], fin[io, jd])
    if best:
        print(f"최적: 용량={best[1]} off={best[2]} → 독성 {best[0]:.0f}, "
              f"내성 {best[3]:.2f}, 최종 {best[4]:.1f}x")
    # 참고: 전량(1.0)+기본 band(off 0.5)
    jd = DOSE_GRID.index(1.0); io = OFF_GRID.index(0.5)
    print(f"전량(1.0)+off0.5: 독성 {tox[io, jd]:.0f}, 최종 {fin[io, jd]:.1f}x, "
          f"통제 {'O' if controlled[io, jd] else 'X'}")

    _figure(tox, fin, controlled, best)
    return dict(tox=tox, fin=fin, controlled=controlled, best=best)


def _figure(tox, fin, controlled, best):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2))

    def heat(ax, M, title, cmap, fmt):
        im = ax.imshow(M, origin="lower", aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(DOSE_GRID)))
        ax.set_xticklabels([f"{d:.0%}" for d in DOSE_GRID])
        ax.set_yticks(range(len(OFF_GRID))); ax.set_yticklabels(OFF_GRID)
        ax.set_xlabel("용량 (전량 대비)")
        ax.set_ylabel("adapt_off (휴약 배수)")
        ax.set_title(title, fontsize=11, fontweight="bold")
        for io in range(len(OFF_GRID)):
            for jd in range(len(DOSE_GRID)):
                ax.text(jd, io, fmt.format(M[io, jd]), ha="center", va="center",
                        fontsize=8.5, color="black")
                if controlled[io, jd]:
                    ax.add_patch(plt.Rectangle((jd - .5, io - .5), 1, 1, fill=False,
                                               edgecolor="#27AE60", lw=2.5))
        fig.colorbar(im, ax=ax, fraction=0.046)

    heat(a1, tox, "누적 독성(부담) — 초록테두리=통제됨\n낮을수록 좋음", "YlOrRd", "{:.0f}")
    heat(a2, fin, "최종 종양(초기 대비 배수)\n통제 여부", "RdYlGn_r", "{:.1f}")
    if best:
        jb = DOSE_GRID.index(best[1]); ib = OFF_GRID.index(best[2])
        for ax in (a1, a2):
            ax.scatter([jb], [ib], marker="*", s=380, color="#2E86AB",
                       edgecolor="white", linewidth=1.3, zorder=5)
    fig.suptitle(f"용량 × band 동시 최적화 (커큐민+마늘+Rg3) — ★ 통제 유지 최소 부담: "
                 f"용량 {best[1]:.0%}/off {best[2]} (독성 {best[0]:.0f})"
                 if best else "용량 × band 최적화", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "dose_band_optimization.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
