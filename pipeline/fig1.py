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
    "D8_coexistence_natural_control": "공존/통제 × 천연물",
    "D5_foodmed_cancer_model": "약식동원 × 암 × 모델",
    "D7_spatial_abm_pdac": "공간전사체 × ABM × PDAC",
    "D4_mycaf_natural_pdac": "myCAF × 천연물 × PDAC",
    "D3_netpharm_pdac_caf": "네트워크약리 × PDAC × CAF",
    "D6_adaptive_pdac_resistance": "적응요법 × PDAC × 내성",
    "D2_abm_natural_tme": "ABM/계산 × 천연물 × TME",
    "D1_adaptive_natural": "적응요법 × 천연물",
    "INTEGRATED_A_model_natural_caf": "통합A: 모델 × 천연물 × CAF",
    "INTEGRATED_B_adaptive_natural_caf": "통합B: 적응요법 × 천연물 × CAF",
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
    ax.set_xlabel("PubMed 히트 수 (체계검색)")
    ax.set_xlim(0, 66)
    ax.set_title("a  차원별 문헌량 — 개별 축은 포화, 통합 쿼리는 극소",
                 fontsize=10.5, fontweight="bold")
    # 범례
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#34618A", label="개별 축(기존 분야)"),
                       Patch(color="#C0392B", label="통합(3축 동시) — off-target")],
              fontsize=8, frameon=False, loc="lower right")


def panel_b(ax):
    R = 1.35
    ctr = {"①": (-0.7, 0.45), "②": (0.7, 0.45), "③": (0.0, -0.75)}
    cols = {"①": "#E67E22", "②": "#16A085", "③": "#8E44AD"}
    names = {"①": "① 정적 네트워크약리학\n(약식동원 화합물)\n약 33편",
             "②": "② 한약-CAF wet-lab\n(PDAC)\n약 7편",
             "③": "③ 공간/ABM 모델\n(PDAC TME)\n약 31편"}
    labpos = {"①": (-1.9, 1.5), "②": (1.9, 1.5), "③": (0.0, -2.35)}
    for k, (cx, cy) in ctr.items():
        ax.add_patch(Circle((cx, cy), R, color=cols[k], alpha=0.22,
                            ec=cols[k], lw=1.5))
        ax.text(*labpos[k], names[k], ha="center", va="center", fontsize=8.8,
                fontweight="bold", color=cols[k])
    # 중앙 교집합 = 본 연구 (공백)
    ax.scatter([0], [0.05], marker="*", s=650, color="#C0392B",
               edgecolor="black", linewidth=0.8, zorder=5)
    ax.text(0, -0.42, "본 연구\n(선행 0편)", ha="center", va="center",
            fontsize=9.2, fontweight="bold", color="#C0392B")
    ax.text(0, 2.55, "세 진영을 잇는 동적·공간·적응통제 프레임워크는 공백",
            ha="center", fontsize=9, color="#555")
    ax.set_xlim(-3.1, 3.1); ax.set_ylim(-3.0, 3.0); ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("b  3진영과 빈 교집합 (통합 쿼리 2·3편은 모두 off-target)",
                 fontsize=10.5, fontweight="bold")


def main():
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.6),
                            gridspec_kw=dict(width_ratios=[1.05, 1]))
    panel_a(axs[0]); panel_b(axs[1])
    fig.suptitle("Figure 1 — 신규성: 구성요소는 포화, 통합은 미개척 "
                 "(체계적 PubMed 검색, 10 쿼리)", fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(ROOT, "assets", "fig1_novelty.png")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
