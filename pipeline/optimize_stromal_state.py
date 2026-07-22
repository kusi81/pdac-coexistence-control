"""L1 핵심 시연 — '최적 기질 상태(optimal stromal state)'.

리뷰어 핵심지적: myCAF containment가 모델에 없었음. 이제 abm.py에 구현:
  myCAF 밀도가 ① 종양 확장(딸세포 배치)을 물리적으로 가두고(confinement)
  ② 약물 침투를 저해하며 ③ CD8 침투를 배제한다(기존).
→ 트레이드오프: myCAF↑ = 가둠(좋음) vs 면역배제+약물차단(나쁨).

이 스크립트: 고정 치료 하에서 myCAF 장벽 수준(k_caf_activate)을 스윕 → 종양 부담이
U자형(양끝보다 중간이 최소)이면, '제거가 아닌 최적 조절'이라는 중심가설이 모델에서
직접 성립함을 보인다. (항섬유화 용량 = k_caf_activate를 낮추는 것과 동치)
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
matplotlib.rcParams["axes.unicode_minus"] = False

from synthetic import make_tissue
from abm import simulate, control_metrics

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# treatment strong enough to matter but not to eradicate alone → confinement의 기여가 드러남
PARAMS = dict(k_prolif=0.17, cd8_recruit=6, k_kill=0.40)
DRUG = [("generic_cytotoxic", 0.30)]
DAYS = 110
SEEDS = [1, 2, 3, 4, 5]
ACT_GRID = [0.0, 0.02, 0.04, 0.06, 0.08, 0.11, 0.15, 0.20, 0.28]


def main():
    c, l, _ = make_tissue("contained", seed=42)
    n0 = (l == "Tumor").sum()
    rows = []
    print("myCAF 장벽수준(k_caf_activate) → 최종 종양/초기 (5-seed 평균±sd)")
    for act in ACT_GRID:
        vals = []
        for s in SEEDS:
            p = dict(PARAMS); p["k_caf_activate"] = act; p["seed"] = s
            h, _ = simulate(c, l, days=DAYS, params=p, regimen_subs=DRUG,
                            snapshots=(0.0, 1.0))
            vals.append(control_metrics(h, n0=n0)["final_frac"])
        m, sd = float(np.mean(vals)), float(np.std(vals))
        rows.append((act, m, sd))
        print(f"  {act:.2f} → {m:.2f} ± {sd:.2f}x")

    acts = np.array([r[0] for r in rows])
    means = np.array([r[1] for r in rows])
    sds = np.array([r[2] for r in rows])
    imin = int(np.argmin(means))
    np.save(os.path.join(ROOT, "data", "stromal_sweep.npy"),
            np.column_stack([acts, means, sds]))

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    ax.errorbar(acts, means, yerr=sds, fmt="o-", color="#2E7D5B", lw=2,
                capsize=3, markersize=7, zorder=3)
    ax.scatter([acts[imin]], [means[imin]], marker="*", s=520, color="#C0392B",
               edgecolor="black", linewidth=0.8, zorder=5)
    ax.annotate(f"optimal stromal state\n(k_caf_activate={acts[imin]:.2f})",
                (acts[imin], means[imin]), xytext=(acts[imin] + 0.03, means[imin] - 0.18),
                fontsize=9, color="#C0392B", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#C0392B"))
    # arm annotations
    ax.annotate("stroma depletion →\ntumor spatial escape\n(no confinement)",
                (acts[0], means[0]), xytext=(acts[0] + 0.005, means[0] + 0.12),
                fontsize=8.5, color="#555")
    ax.annotate("excess stroma →\nimmune exclusion +\nimpaired drug delivery",
                (acts[-1], means[-1]), xytext=(acts[-1] - 0.085, means[-1] + 0.10),
                fontsize=8.5, color="#555", ha="left")
    ax.set_xlabel("myCAF barrier level  (k_caf_activate; ← stronger anti-fibrotic)")
    ax.set_ylabel("Final tumor burden (fold vs initial)  →  lower is better")
    ax.set_title("The myCAF barrier is a controllable resource with an optimum\n"
                 "under fixed treatment, intermediate stroma controls best (in-silico)",
                 fontsize=11.5, fontweight="bold")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out = os.path.join(ROOT, "assets", "stromal_state_optimum.png")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"\n최적: k_caf_activate={acts[imin]:.2f} → {means[imin]:.2f}x "
          f"(양끝: {means[0]:.2f}x / {means[-1]:.2f}x)")
    print("wrote", out)


if __name__ == "__main__":
    main()
