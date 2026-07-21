"""
적응형 band/타이밍 세부 최적화 — adapt_on × adapt_off 격자 탐색.

상위 자연물 조합(마늘+쑥)에 대해, 투여 band를 바꿔가며:
  - adapt_on: 종양이 초기의 몇 배를 넘으면 투여 시작 (높을수록 늦게 개입 → 부담↓·escape위험↑)
  - adapt_off: 몇 배 밑으로 떨어지면 휴약 (낮을수록 강하게 밀어붙임 → 부담↑)
를 조합해 '통제 유지(진행 없음) 하 최소 부담' band를 찾는다.

목적(제약 최적화): 통제됨(관측기간 내 진행 없음 & 최종<1.2x) 중 누적독성 최소.
⚠️ 가설 샌드박스.
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

PARAMS = dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10,
              init_resistant_frac=0.03, mutation_rate=0.003,
              resistant_immune_evasion=0.35)
COMBO = [("garlic", 1.0), ("mugwort", 1.0)]   # 마늘+쑥 (항CAF+세포독성)
DAYS = 150
ON_GRID = [1.05, 1.15, 1.25, 1.4, 1.6]
OFF_GRID = [0.4, 0.5, 0.6, 0.7, 0.85]
CONTROL_MAX = 1.2   # 최종 종양이 이 배수 미만 & 진행없음 = '통제됨'


def main():
    c, l, _ = make_tissue("contained", seed=42)
    n0 = simulate(c, l, days=2, params=PARAMS, snapshots=(0.0,))[0][0]["n_tumor"]

    tox = np.full((len(OFF_GRID), len(ON_GRID)), np.nan)
    fin = np.full_like(tox, np.nan)
    res = np.full_like(tox, np.nan)
    controlled = np.zeros_like(tox, dtype=bool)

    print(f"격자 {len(ON_GRID)}×{len(OFF_GRID)} 탐색 (마늘+쑥, {DAYS}일)...")
    for j, on in enumerate(ON_GRID):
        for i, off in enumerate(OFF_GRID):
            if off >= on:
                continue
            h, _ = simulate(c, l, days=DAYS, params=PARAMS, regimen_subs=COMBO,
                            adaptive=True, adapt_on=on, adapt_off=off,
                            synergy=0.3, snapshots=(0.0, 1.0))
            m = control_metrics(h, n0=n0)
            tox[i, j] = m["cum_toxicity"]
            fin[i, j] = m["final_frac"]
            res[i, j] = m["final_resistant_frac"]
            ok = m["progression_censored"] and m["final_frac"] < CONTROL_MAX
            controlled[i, j] = ok

    # 최적: 통제됨 중 독성 최소
    best = None
    for j, on in enumerate(ON_GRID):
        for i, off in enumerate(OFF_GRID):
            if controlled[i, j]:
                if best is None or tox[i, j] < best[0]:
                    best = (tox[i, j], on, off, fin[i, j], res[i, j])
    if best:
        print(f"\n최적 band: adapt_on={best[1]} adapt_off={best[2]} "
              f"→ 독성 {best[0]:.0f}, 최종종양 {best[3]:.2f}x, 내성 {best[4]:.2f}")
    # 기본값(1.1/0.7) 참고
    hb, _ = simulate(c, l, days=DAYS, params=PARAMS, regimen_subs=COMBO,
                     adaptive=True, adapt_on=1.1, adapt_off=0.7, synergy=0.3,
                     snapshots=(0.0, 1.0))
    mb = control_metrics(hb, n0=n0)
    print(f"기본 band(1.1/0.7): 독성 {mb['cum_toxicity']:.0f}, 최종 {mb['final_frac']:.2f}x")

    _figure(tox, fin, controlled, best, n0)
    return dict(tox=tox, fin=fin, controlled=controlled, best=best)


def _figure(tox, fin, controlled, best, n0):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2))
    ext = [0, len(ON_GRID), 0, len(OFF_GRID)]

    def heat(ax, M, title, cmap, fmt):
        im = ax.imshow(M, origin="lower", aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(ON_GRID))); ax.set_xticklabels(ON_GRID)
        ax.set_yticks(range(len(OFF_GRID))); ax.set_yticklabels(OFF_GRID)
        ax.set_xlabel("adapt_on (투여 시작 배수)")
        ax.set_ylabel("adapt_off (휴약 배수)")
        ax.set_title(title, fontsize=11, fontweight="bold")
        for i in range(len(OFF_GRID)):
            for j in range(len(ON_GRID)):
                if not np.isnan(M[i, j]):
                    ax.text(j, i, fmt.format(M[i, j]), ha="center", va="center",
                            fontsize=8,
                            color="white" if not controlled[i, j] else "black")
                # 통제됨 셀에 테두리
                if controlled[i, j]:
                    ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                               edgecolor="#27AE60", lw=2.5))
        fig.colorbar(im, ax=ax, fraction=0.046)

    heat(a1, tox, "누적 독성(부담) — 초록테두리=통제됨\n낮을수록 좋음", "YlOrRd", "{:.0f}")
    heat(a2, fin, "최종 종양(초기 대비 배수)\n통제 여부", "RdYlGn_r", "{:.1f}")
    if best:
        jb = ON_GRID.index(best[1]); ib = OFF_GRID.index(best[2])
        for ax in (a1, a2):
            ax.scatter([jb], [ib], marker="*", s=350, color="#2E86AB",
                       edgecolor="white", linewidth=1.2, zorder=5)
    fig.suptitle(f"적응형 band 최적화 (마늘+쑥) — ★=최적(통제 유지 최소 부담): "
                 f"on={best[1]}/off={best[2]}, 독성{best[0]:.0f}" if best else
                 "적응형 band 최적화 (마늘+쑥)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "band_optimization.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
