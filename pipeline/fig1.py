"""Figure 1 — 신규성 landscape (본문 §3.1).

  1a  차원별 PubMed 히트: 개별 축은 포화, 통합 쿼리는 극소(off-target)
  1b  3진영 벤: ①정적 네트워크약리 ②한약-CAF wet-lab ③공간/ABM 모델 — 교집합 공백=본 연구

데이터: docs/literature_search/results.json (lit_search.py, 체계검색).
"""
import sys, os, json, warnings
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
matplotlib.rcParams["axes.unicode_minus"] = False
from matplotlib.patches import Circle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Q = json.load(open(os.path.join(ROOT, "docs", "literature_search",
                                "results.json"), encoding="utf-8"))["queries"]

LABELS = {
    "D8_coexistence_natural_control": "Coexistence/control × natural cpd.",
    "D5_foodmed_cancer_model": "Food-medicine × cancer × model",
    "D7_spatial_abm_pdac": "Spatial omics × ABM × PDAC",
    "D4_mycaf_natural_pdac": "myCAF × natural cpd. × PDAC",
    "D3_netpharm_pdac_caf": "Network pharm. × PDAC × CAF",
    "D6_adaptive_pdac_resistance": "Adaptive therapy × PDAC × resist.",
    "D2_abm_natural_tme": "ABM/compute × natural cpd. × TME",
    "D1_adaptive_natural": "Adaptive therapy × natural cpd.",
    "INTEGRATED_A_model_natural_caf": "Integrated A: model × natural × CAF",
    "INTEGRATED_B_adaptive_natural_caf": "Integrated B: adaptive × natural × CAF",
}
INTEG = {"INTEGRATED_A_model_natural_caf", "INTEGRATED_B_adaptive_natural_caf"}


def panel_a(ax):
    items = [(LABELS[k], Q[k]["count"], k in INTEG) for k in LABELS]
    items.sort(key=lambda t: t[1], reverse=True)
    y = np.arange(len(items))
    colors = ["#C0392B" if ig else "#34618A" for _, _, ig in items]
    ax.barh(y, [c for _, c, _ in items], color=colors, edgecolor="black",
            linewidth=0.6)
    for i, (_, c, ig) in enumerate(items):
        ax.text(c + 1.2, i, f"{c}" + ("  (off-target)" if ig else ""),
                va="center", fontsize=8.5,
                color="#C0392B" if ig else "#333")
    ax.set_yticks(y); ax.set_yticklabels([lab for lab, _, _ in items], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("PubMed hit count (systematic survey)")
    ax.set_xlim(0, 66)
    ax.set_title("a  Literature by dimension — single axes saturated, integrated queries near-empty",
                 fontsize=10, fontweight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#34618A", label="Single axis (established field)"),
                       Patch(color="#C0392B", label="Integrated (3 axes) — off-target")],
              fontsize=8, frameon=False, loc="lower right")


def panel_b(ax):
    R = 1.35
    ctr = {"①": (-0.7, 0.45), "②": (0.7, 0.45), "③": (0.0, -0.75)}
    cols = {"①": "#E67E22", "②": "#16A085", "③": "#8E44AD"}
    names = {"①": "① Static network\npharmacology\n(food-medicine)\n~33 papers",
             "②": "② Herbal-CAF\nwet-lab (PDAC)\n~7 papers",
             "③": "③ Spatial/ABM\nmodels (PDAC TME)\n~31 papers"}
    labpos = {"①": (-1.95, 1.55), "②": (1.95, 1.55), "③": (0.0, -2.4)}
    for k, (cx, cy) in ctr.items():
        ax.add_patch(Circle((cx, cy), R, color=cols[k], alpha=0.22,
                            ec=cols[k], lw=1.5))
        ax.text(*labpos[k], names[k], ha="center", va="center", fontsize=8.5,
                fontweight="bold", color=cols[k])
    # central intersection = this work (empty)
    ax.scatter([0], [0.05], marker="*", s=650, color="#C0392B",
               edgecolor="black", linewidth=0.8, zorder=5)
    ax.text(0, -0.42, "This work\n(0 prior)", ha="center", va="center",
            fontsize=9.2, fontweight="bold", color="#C0392B")
    ax.text(0, 2.6, "The dynamic, spatial, adaptive-control framework bridging all three is empty",
            ha="center", fontsize=8.2, color="#555")
    ax.set_xlim(-3.1, 3.1); ax.set_ylim(-3.0, 3.0); ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("b  Three camps, no directly matching integration (integrated queries: 2-3, all off-target)",
                 fontsize=9.6, fontweight="bold")


def main():
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.6),
                            gridspec_kw=dict(width_ratios=[1.05, 1]))
    panel_a(axs[0]); panel_b(axs[1])
    fig.suptitle("Figure 1 — Novelty: components saturated, no directly matching integration "
                 "(targeted PubMed survey, 12 queries: 10 dimensional + 2 integrated)",
                 fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(ROOT, "assets", "fig1_novelty.png")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
