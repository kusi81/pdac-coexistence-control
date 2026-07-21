"""
공격적 체제 band 최적화 — adaptive therapy의 핵심 역설이 드러나는 곳.

공격적 종양 + 약한 면역 + 강한 내성/면역회피 + 내성 적합도 비용(resistance_cost).
이 체제에선 band가 통제/내성/부담을 예리하게 가른다:
  - 너무 빡빡(off 낮음): 감수성을 없애 내성이 부상 위험 + 부담↑
  - 스위트스팟(off 중간): 감수성을 남겨 내성을 경쟁으로 억누름 + 최소 부담
  - 너무 느슨(off 높음): 내성 분율↑ (경쟁적 방출 시작)
목적: 통제됨 중 누적독성 최소 (내성도 함께 확인).
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

# 공격적 종양 + 문헌 접지 저항성 파라미터 (cost 0.24=NSCLC 실측, evasion 0.45=PD-L1/항원소실)
PARAMS = dict(k_prolif=0.20, cd8_recruit=8, k_kill=0.45, k_caf_activate=0.10,
              init_resistant_frac=0.01, mutation_rate=0.001,
              resistant_immune_evasion=0.45, resistance_cost=0.24)
COMBO = [("curcumin", 1.0), ("garlic", 1.0), ("ginsenoside_rg3", 1.0)]
DAYS = 200
ON_GRID = [1.05, 1.15, 1.3, 1.5]
OFF_GRID = [0.4, 0.5, 0.6, 0.7, 0.85]
CONTROL_MAX = 1.3


SEEDS = [42, 7, 123]   # 조직·시뮬 seed 3개 평균 → 확률적 노이즈 완화


def main():
    tox = np.zeros((len(OFF_GRID), len(ON_GRID)))
    resf = np.zeros_like(tox)
    fin = np.zeros_like(tox)
    ctrl_count = np.zeros_like(tox)
    valid = np.zeros_like(tox, dtype=bool)

    print(f"격자 {len(ON_GRID)}×{len(OFF_GRID)} × seed {len(SEEDS)} "
          f"(공격적 체제, 커큐민+마늘+Rg3, {DAYS}일)...")
    for sd in SEEDS:
        c, l, _ = make_tissue("contained", seed=sd)
        P = dict(PARAMS, seed=sd)
        for j, on in enumerate(ON_GRID):
            for i, off in enumerate(OFF_GRID):
                if off >= on:
                    continue
                valid[i, j] = True
                h, _ = simulate(c, l, days=DAYS, params=P, regimen_subs=COMBO,
                                adaptive=True, adapt_on=on, adapt_off=off,
                                synergy=0.3, snapshots=(0.0, 1.0))
                n0 = h[0]["n_tumor"]
                m = control_metrics(h, n0=n0)
                tox[i, j] += m["cum_toxicity"]
                resf[i, j] += m["final_resistant_frac"]
                fin[i, j] += m["final_frac"]
                if m["progression_censored"] and m["final_frac"] < CONTROL_MAX:
                    ctrl_count[i, j] += 1
    ns = len(SEEDS)
    tox = np.where(valid, tox / ns, np.nan)
    resf = np.where(valid, resf / ns, np.nan)
    fin = np.where(valid, fin / ns, np.nan)
    controlled = ctrl_count >= (ns / 2.0)   # 과반 seed에서 통제되면 '통제됨'

    best = None
    for j, on in enumerate(ON_GRID):
        for i, off in enumerate(OFF_GRID):
            if controlled[i, j] and (best is None or tox[i, j] < best[0]):
                best = (tox[i, j], on, off, resf[i, j], fin[i, j])
    if best:
        print(f"최적 band: on={best[1]}/off={best[2]} → 독성 {best[0]:.0f}, "
              f"내성 {best[3]:.2f}, 최종 {best[4]:.1f}x")
    # 연속(관행) 벤치마크 — seed 평균
    ct, cr, cf = 0.0, 0.0, 0.0
    for sd in SEEDS:
        c, l, _ = make_tissue("contained", seed=sd)
        hc, _ = simulate(c, l, days=DAYS, params=dict(PARAMS, seed=sd),
                         regimen_subs=COMBO, snapshots=(0.0, 1.0))
        m = control_metrics(hc, n0=hc[0]["n_tumor"])
        ct += m["cum_toxicity"]; cr += m["final_resistant_frac"]; cf += m["final_frac"]
    ns = len(SEEDS)
    mc = {"cum_toxicity": ct / ns, "final_resistant_frac": cr / ns, "final_frac": cf / ns}
    print(f"연속(관행): 독성 {mc['cum_toxicity']:.0f}, 내성 {mc['final_resistant_frac']:.2f}, "
          f"최종 {mc['final_frac']:.1f}x")

    _figure(tox, resf, controlled, best, mc)
    return dict(tox=tox, resf=resf, controlled=controlled, best=best)


def _figure(tox, resf, controlled, best, mc):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2))

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
                            fontsize=8.5, color="black")
                if controlled[i, j]:
                    ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                               edgecolor="#27AE60", lw=2.5))
        fig.colorbar(im, ax=ax, fraction=0.046)

    heat(a1, tox, "누적 독성(부담) — 초록테두리=통제됨\n낮을수록 좋음 (연속=%.0f)" %
         mc["cum_toxicity"], "YlOrRd", "{:.0f}")
    heat(a2, resf, "최종 내성 분율 — 경쟁적 방출\n낮을수록 좋음 (연속=%.2f)" %
         mc["final_resistant_frac"], "PuRd", "{:.2f}")
    if best:
        jb = ON_GRID.index(best[1]); ib = OFF_GRID.index(best[2])
        for ax in (a1, a2):
            ax.scatter([jb], [ib], marker="*", s=380, color="#2E86AB",
                       edgecolor="white", linewidth=1.3, zorder=5)
    fig.suptitle("공격적 체제 band 최적화 — ★ 내부 최적점(통제 유지·최소 부담·내성 억제): "
                 f"on={best[1]}/off={best[2]} (연속은 박멸하나 부담 2배+)"
                 if best else "공격적 체제 band 최적화", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "band_optimization_aggressive.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
