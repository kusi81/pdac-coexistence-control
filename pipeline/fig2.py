"""Figure 2 조립 — 지표 검증 (본문 §3.2 참조).

  2a  합성 조직 특이도: contained(myCAF 링) vs diffuse(무작위) + barrier z
  2b  SCOTIA 저자주석 양성대조 통과: myCAF-iCAF 종양거리(µm), <0=myCAF 더 가까움(정답)
  2c  Xenium module-score 주석 실패: 같은 지표가 실패(주석 문제)

2b: data/scotia/scotia_posctrl.csv. 2c: refine_annotate.py 결과(임베드, 재현가능).
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt
matplotlib.rcParams["axes.unicode_minus"] = False   # − 글리프 문제 해결
from matplotlib import gridspec
from synthetic import make_tissue
from spatial_core import barrier_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CBG = {"Tumor": "#34495E", "myCAF": "#C0392B", "iCAF": "#E67E22",
       "CD8_T": "#2980B9"}
TREAT_C = {"Untreated": "#2980B9", "CRT": "#C0392B", "CRTL": "#8E44AD"}

# 2c: Xenium refined-annotation proximity (myCAF - iCAF, µm) — refine_annotate.py
XENIUM = [("31076", "naive", 11.8), ("28429", "naive", 7.0),
          ("38245", "CRT", 4.4), ("35406", "CRT", 73.3), ("39928", "naive", -7.2)]


def _scatter(ax, mode, title):
    coords, labels, _ = make_tissue(mode=mode, seed=42)
    bs = barrier_score(coords, labels, tumor="Tumor", barrier="myCAF",
                       immune="CD8_T", corridor_um=30, seed=0)
    z = bs.get("z_score") if bs else None
    for ct in ["Tumor", "myCAF", "CD8_T"]:
        m = labels == ct
        ax.scatter(coords[m, 0], coords[m, 1], s=2.5, c=CBG[ct],
                   label=ct, alpha=0.7, linewidths=0)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    zt = f"barrier z = {z:+.1f}" if z is not None else "barrier z = n/a"
    verdict = "barrier interposed (containment)" if (z or 0) > 2 else "random (no barrier)"
    ax.set_title(f"{title}\n{zt}  ·  {verdict}", fontsize=10, fontweight="bold")
    return ax


def _strip(ax, rows, title, xlim, note):
    """rows: list of (label, group, diff). diff<0 = myCAF closer = 정답."""
    ax.axvspan(xlim[0], 0, color="#2ECC71", alpha=0.08)   # 정답 영역
    ax.axvspan(0, xlim[1], color="#E74C3C", alpha=0.08)   # 오답 영역
    ax.axvline(0, color="black", lw=1.2)
    rng = np.random.default_rng(1)
    seen = {}
    for lab, grp, d in rows:
        y = rng.uniform(-0.35, 0.35)
        ax.scatter(d, y, s=48, c=TREAT_C.get(grp, "#666"),
                   edgecolor="black", linewidth=0.5, zorder=3,
                   label=grp if grp not in seen else None)
        seen[grp] = 1
    npass = sum(1 for _, _, d in rows if d < 0)
    ax.set_xlim(*xlim); ax.set_ylim(-0.6, 0.6); ax.set_yticks([])
    ax.set_title(f"{title}  ({npass}/{len(rows)} pass positive control)",
                 fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=7.5, loc="upper left", frameon=False)


def main():
    df = pd.read_csv(os.path.join(ROOT, "data", "scotia", "scotia_posctrl.csv"))
    scotia_rows = [(r["sample"], r["treat"], r["diff_um"]) for _, r in df.iterrows()]

    fig = plt.figure(figsize=(12.5, 6.4))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1.25], height_ratios=[1, 1],
                           hspace=0.42, wspace=0.22)
    # 2a
    _scatter(fig.add_subplot(gs[0, 0]), "contained",
             "a  Synthetic: contained (myCAF ring)")
    axd = _scatter(fig.add_subplot(gs[1, 0]), "diffuse",
                   "Synthetic: diffuse (random)")
    axd.legend(fontsize=7.5, loc="upper right", frameon=False,
               markerscale=2, ncol=3, bbox_to_anchor=(1.0, -0.02))
    # 2b SCOTIA
    XR = (-78, 80)
    _strip(fig.add_subplot(gs[0, 1]), scotia_rows,
           "b  SCOTIA author annotation (CosMx)", XR, "pass")
    # 2c Xenium
    ax2c = fig.add_subplot(gs[1, 1])
    _strip(ax2c, [(s, ("CRT" if g == "CRT" else "Untreated"), d)
                  for s, g, d in XENIUM],
           "c  Xenium module-score annotation", XR, "fail")
    ax2c.set_xlabel("myCAF - iCAF distance to tumor (µm)  ·  "
                    "<0 = myCAF closer (correct)", fontsize=9.5)

    fig.suptitle("Figure 2 — the containment metric is sound; annotation is the limiting factor "
                 "(synthetic specificity · SCOTIA positive control · Xenium annotation failure)",
                 fontsize=11.5, fontweight="bold", y=0.99)
    out = os.path.join(ROOT, "assets", "fig2_validation.png")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
