"""공간 ABM vs well-mixed(평균장) 모델 — 공간모델의 필요성 입증 (리뷰어 #2).

리뷰어 Major-1: "공간 geometry가 (단순 종양수가 아니라) 선택되는 전략을 바꾸는가?
공간모델이 아니면 못 얻는 결과인가?" 이를 well-mixed 대조로 입증한다.

well-mixed 모델(ODE): 같은 기전(증식·면역살상·내성·약물)을 갖되 **공간 없음** —
myCAF는 오직 aggregate 분율 M을 통해 면역살상을 감쇠(exp(-alpha·M))할 뿐, 국소
confinement(공간 가둠)나 corridor geometry가 없다. 두 가지를 보인다:

  A. 기하 맹목(geometry-blindness): contained와 diffuse는 세포수(=aggregate M)가
     동일 → well-mixed는 **동일 결과**를 낸다. 반면 공간 ABM은 면역요법 효능이
     ~12배 다르다(Fig S17/predict_biomarker). 평균장은 결정적 정보를 버린다.
  B. 전략 역전: well-mixed에서 stromal level M을 스윕하면 최적 M은 **항상 0**
     (M은 면역배제 비용만 있고 confinement 이득이 없으므로) → 평균장은 늘 '기질 제거'만
     권한다. 공간 ABM은 조건부 비영 최적(keep-stroma regime; Fig 4)을 낸다.

결과: assets/wellmixed_compare.png (+ 콘솔 요약). ODE만 쓰므로 가볍다.
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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def wellmixed(M, alpha=2.2, days=90, dt=0.25, immune=True, drug_kill=0.0,
              k_prolif=0.15, k_kill=0.6, K=1.0, cost=0.24, evasion=0.4,
              init_res=0.03):
    """평균장 ODE. M=aggregate myCAF(0-1). 면역살상은 exp(-alpha·M)로 감쇠(공간 없음).
    confinement 항 없음. 반환: 최종 종양/초기."""
    S = 1.0 - init_res; R = init_res
    n0 = S + R
    kill0 = k_kill if immune else 0.0
    kill_eff = kill0 * np.exp(-alpha * M)          # 평균장 면역배제(기하 무관)
    for _ in range(int(days / dt)):
        tot = S + R
        grow_s = k_prolif * S * (1 - tot / K)
        grow_r = k_prolif * (1 - cost) * R * (1 - tot / K)
        kill_s = kill_eff * S
        kill_r = kill_eff * (1 - evasion) * R
        drug_s = drug_kill * S
        S += (grow_s - kill_s - drug_s) * dt
        R += (grow_r - kill_r) * dt
        S = max(S, 0.0); R = max(R, 0.0)
    return (S + R) / n0


def main():
    # ── A. geometry-blindness ──
    # contained/diffuse는 make_tissue에서 세포수 동일(myCAF abundance≈0.29) → 같은 M.
    M_shared = 0.29
    wm_contained = wellmixed(M_shared, immune=True)
    wm_diffuse = wellmixed(M_shared, immune=True)     # 동일 M → 동일
    # 공간 ABM 실측(predict_biomarker.csv): contained-full vs diffuse-full 면역요법 최종
    abm_c, abm_d = np.nan, np.nan
    csvp = os.path.join(ROOT, "data", "predict_biomarker.csv")
    if os.path.exists(csvp):
        df = pd.read_csv(csvp)
        cc = df[df.arch.str.contains("Contained .full", regex=True)]
        dd = df[df.arch.str.contains("Diffuse")]
        abm_c = float(cc["immune_final"].mean()) if len(cc) else np.nan
        abm_d = float(dd["immune_final"].mean()) if len(dd) else np.nan
    print("[A] geometry-blindness (immune-therapy final tumor):")
    print(f"  well-mixed: contained={wm_contained:.2f}  diffuse={wm_diffuse:.2f}  "
          f"(identical — same aggregate M)")
    print(f"  spatial ABM: contained={abm_c:.2f}  diffuse={abm_d:.2f}  "
          f"({abm_c/abm_d:.1f}× difference)" if abm_d else "")

    # ── B. optimal stromal level vs immune-exclusion ──
    Ms = np.linspace(0, 0.5, 11)
    curves = {}
    for alpha in (1.0, 2.2, 3.5):
        curves[alpha] = [wellmixed(M, alpha=alpha, immune=True) for M in Ms]
        opt = Ms[int(np.argmin(curves[alpha]))]
        print(f"[B] well-mixed alpha={alpha}: optimal stromal M = {opt:.2f} "
              f"({'reduce (M=0)' if opt == 0 else 'interior'})")

    figure(wm_contained, wm_diffuse, abm_c, abm_d, Ms, curves)


def figure(wm_c, wm_d, abm_c, abm_d, Ms, curves):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.4))

    # a: geometry-blindness
    labels = ["Well-mixed\ncontained", "Well-mixed\ndiffuse",
              "Spatial ABM\ncontained", "Spatial ABM\ndiffuse"]
    vals = [wm_c, wm_d, abm_c, abm_d]
    cols = ["#95A5A6", "#95A5A6", "#C0392B", "#2980B9"]
    x = np.arange(4)
    a1.bar(x, vals, color=cols, edgecolor="black", linewidth=0.7, alpha=0.85)
    for i, v in enumerate(vals):
        a1.text(i, v + 0.02, f"{v:.2f}×", ha="center", fontsize=9)
    a1.axhline(1.5, color="#E74C3C", ls=":", lw=1)
    a1.set_xticks(x); a1.set_xticklabels(labels, fontsize=8.5)
    a1.set_ylabel("Tumor burden under immune therapy (fold)")
    a1.set_title("a  Geometry-blindness of the mean-field model\n"
                 "well-mixed can't tell contained from diffuse (same abundance);\n"
                 "the spatial ABM predicts a large difference", fontsize=9.8,
                 fontweight="bold")
    a1.annotate("", xy=(3, abm_d), xytext=(2, abm_c),
                arrowprops=dict(arrowstyle="<->", color="#555", lw=1.2))
    if abm_d:
        a1.text(2.5, (abm_c + abm_d) / 2 + 0.06, f"{abm_c/abm_d:.0f}×",
                ha="center", fontsize=9, color="#555", fontweight="bold")

    # b: optimal stromal level — well-mixed monotonic (no interior optimum)
    for alpha, col in zip((1.0, 2.2, 3.5), ("#27AE60", "#E67E22", "#8E44AD")):
        a2.plot(Ms, curves[alpha], "o-", color=col, lw=2,
                label=f"immune exclusion α={alpha}")
        oi = int(np.argmin(curves[alpha]))
        a2.scatter([Ms[oi]], [curves[alpha][oi]], s=90, facecolor="none",
                   edgecolor=col, linewidths=2, zorder=5)
    a2.set_xlabel("Aggregate stromal level M (myCAF)")
    a2.set_ylabel("Final tumor burden (fold)")
    a2.set_title("b  Mean-field always prefers stromal reduction (optimum at M=0)\n"
                 "no confinement benefit exists without geometry → no keep-stroma regime\n"
                 "(the spatial ABM does produce a non-zero optimum; Fig. 4)",
                 fontsize=9.6, fontweight="bold")
    a2.legend(frameon=False, fontsize=8.5)

    fig.suptitle("Why the spatial model is necessary: a well-mixed (mean-field) model with "
                 "the same mechanisms discards the geometry that drives the strategy "
                 "(in silico)", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "assets", "wellmixed_compare.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
