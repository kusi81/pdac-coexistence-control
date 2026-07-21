"""
자연물 배합 × 적응형 최적화 — 사용자 궁극 목표.

질문: 저독성 자연물(전통식품 성분) 조합을 적응형 on/off로 투여하면, 관행 화학요법보다
      훨씬 적은 부담으로 종양을 오래 통제(공존)할 수 있는가? 어떤 배합이 최선인가?

평가: control_score = TTP / (누적독성+1)  (적은 부담으로 오래 통제할수록 큼).
      + TTP(진행까지 시간), 최종종양배수, 내성분율, 독성.

⚠️ 가설 샌드박스. 파라미터는 문헌에서 방향만 접지. 사람 효과 아님.
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

# escape 가능 체제 (약물이 '통제는 하되 박멸은 어려운' 창)
PARAMS = dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10,
              init_resistant_frac=0.03, mutation_rate=0.003,
              resistant_immune_evasion=0.35)
DAYS = 150
ADAPT_ON, ADAPT_OFF = 1.1, 0.7

# Food-medicine-homology candidate regimens (mechanism-based); all low-toxicity.
NATURAL = {
    "Curcumin": [("curcumin", 1.0)],
    "Garlic": [("garlic", 1.0)],
    "Mugwort (anti-CAF)": [("mugwort", 1.0)],
    "Wild ginseng": [("wild_ginseng", 1.0)],
    "Garlic + Mugwort": [("garlic", 1.0), ("mugwort", 1.0)],
    "Ginseng + Garlic + Mugwort": [("wild_ginseng", 1.0), ("garlic", 1.0), ("mugwort", 1.0)],
    "Curcumin + Garlic + Rg3": [("curcumin", 1.0), ("garlic", 1.0), ("ginsenoside_rg3", 1.0)],
    "Danshen + Astragaloside (anti-fibrotic)": [("danshen", 1.0), ("astragaloside", 1.0)],
    "5-compound mix": [("curcumin", 1.0), ("garlic", 1.0), ("mugwort", 1.0),
                       ("ginsenoside_rg3", 1.0), ("sea_cucumber", 1.0)],
}
# Standard-of-care benchmark
BENCH = {
    "Gemcitabine (continuous)": [("gemcitabine", 1.0)],
}
UNTREATED = "Untreated"


def main():
    c, l, _ = make_tissue("contained", seed=42)
    n0 = simulate(c, l, days=1, params=PARAMS, snapshots=(0.0,))[0][0]["n_tumor"]
    # 정확한 n0
    h_probe, _ = simulate(c, l, days=2, params=PARAMS, snapshots=(0.0,))
    n0 = h_probe[0]["n_tumor"]

    results = []
    # 무처치
    h0, _ = simulate(c, l, days=DAYS, params=PARAMS, snapshots=(0.0, 1.0))
    results.append((UNTREATED, "none", h0, control_metrics(h0, n0=n0), "#7F8C8D"))
    # 자연물: 적응형
    for name, reg in NATURAL.items():
        h, _ = simulate(c, l, days=DAYS, params=PARAMS, regimen_subs=reg,
                        adaptive=True, adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF,
                        synergy=0.3, snapshots=(0.0, 1.0))
        results.append((name, "natural-adaptive", h, control_metrics(h, n0=n0), "#27AE60"))
    # 벤치마크: 연속 관행
    for name, reg in BENCH.items():
        h, _ = simulate(c, l, days=DAYS, params=PARAMS, regimen_subs=reg,
                        snapshots=(0.0, 1.0))
        results.append((name, "conventional", h, control_metrics(h, n0=n0), "#C0392B"))

    # 랭킹: 통제 성공(censored) 먼저, 그 중 저독성 순 → '적은 부담으로 통제'가 최상위.
    print(f"{'배합':<20}{'유형':<20}{'TTP':>7}{'최종':>7}{'내성':>7}{'독성':>7}{'통제점수':>9}")
    for name, typ, h, m, _c in sorted(
            results, key=lambda r: (not r[3]["progression_censored"],
                                    r[3]["cum_toxicity"])):
        star = "*" if m["progression_censored"] else " "
        print(f"{name:<19}{typ:<20}{m['ttp_days']:>5.0f}d{star}{m['final_frac']:>6.1f}x"
              f"{m['final_resistant_frac']:>7.2f}{m['cum_toxicity']:>7.0f}"
              f"{m['control_score']:>9.1f}")

    _figure(results, n0)
    return results


def _figure(results, n0):
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 5.4))
    # A: 대표 궤적 (무처치, 젬시타빈, 상위 자연물 2개)
    nat = [r for r in results if r[1] == "natural-adaptive"]
    nat_sorted = sorted(nat, key=lambda r: -r[3]["control_score"])
    show = [r for r in results if r[0] == UNTREATED] + nat_sorted[:2] + \
           [r for r in results if r[1] == "conventional"]
    for name, typ, h, m, col in show:
        t = [x["t"] for x in h]; y = [x["n_tumor"] / n0 for x in h]
        ls = "-" if typ != "conventional" else "--"
        axA.plot(t, y, ls, lw=2, color=col, label=f"{name}")
    axA.axhline(1.5, color="#E74C3C", ls=":", lw=1)
    axA.text(2, 1.53, "progression threshold 1.5x", fontsize=8, color="#E74C3C")
    axA.set_xlabel("Time (days)"); axA.set_ylabel("Tumor / initial")
    axA.set_title("Tumor trajectory — natural adaptive regimens (green) control at low burden",
                  fontsize=10.5, fontweight="bold")
    axA.legend(frameon=False, fontsize=8)

    # B: 트레이드오프 프론티어 — TTP(통제 지속) vs 독성(부담).
    for name, typ, h, m, col in results:
        mk = "o" if typ == "natural-adaptive" else ("s" if typ == "conventional" else "^")
        axB.scatter(m["cum_toxicity"], m["ttp_days"], s=95, color=col,
                    marker=mk, edgecolor="black", linewidth=0.8, zorder=3)
        axB.annotate(name, (m["cum_toxicity"], m["ttp_days"]), fontsize=7.5,
                     xytext=(5, 3), textcoords="offset points")
    axB.axhline(DAYS, color="#27AE60", ls=":", lw=1)
    axB.text(60, DAYS - 6, "controlled over full horizon (no progression)",
             fontsize=8, color="#27AE60")
    axB.set_xlabel("Cumulative toxicity (patient burden)  →  lower is better")
    axB.set_ylabel("Time-to-progression (days)  →  longer is better")
    axB.set_title("Burden vs control-duration frontier\n"
                  "(upper-left = long control at low burden)",
                  fontsize=11, fontweight="bold")
    from matplotlib.lines import Line2D
    axB.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#27AE60",
               markeredgecolor="k", markersize=9, label="Natural, adaptive"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#C0392B",
               markeredgecolor="k", markersize=9, label="Standard chemo (continuous)"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor="#7F8C8D",
               markeredgecolor="k", markersize=9, label="Untreated")],
        frameon=False, fontsize=8, loc="upper right")

    fig.suptitle("Food-medicine-homology regimens × adaptive optimization — "
                 "low-toxicity tumor coexistence (in-silico)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "natural_adaptive_optim.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
