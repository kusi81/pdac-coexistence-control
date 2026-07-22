"""위상도 작도 — data/phase_map.npz에서 읽어 그림만 생성(재계산 없음).

  (a) 각 조건의 최적 myCAF 수준: >0 = 기질유지(자원) 성립, 0 = 기질제거 최선
  (b) 기질 유지의 이득: tumor(myCAF=0) - min tumor. 양수 = 기질이 통제에 도움.
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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = np.load(os.path.join(ROOT, "data", "phase_map.npz"))
T, PRES, ALPHA, MY = d["T"], d["pressure"], d["alpha"], d["mycaf"]

opt = np.array([[MY[int(np.argmin(T[xi, yi]))] for xi in range(len(PRES))]
                for yi in range(len(ALPHA))])          # [alpha, pressure]
benefit = np.array([[T[xi, yi, 0] - T[xi, yi].min() for xi in range(len(PRES))]
                    for yi in range(len(ALPHA))])

fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.6))


def heat(ax, M, title, cmap, fmt, cbar_label):
    im = ax.imshow(M, origin="lower", aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(PRES))); ax.set_xticklabels([f"{p:.1f}" for p in PRES])
    ax.set_yticks(range(len(ALPHA))); ax.set_yticklabels([f"{a:.1f}" for a in ALPHA])
    ax.set_xlabel("Confinement strength (caf_pressure)  →")
    ax.set_ylabel("Immune-exclusion strength (cd8_barrier_alpha)  →")
    for yi in range(len(ALPHA)):
        for xi in range(len(PRES)):
            ax.text(xi, yi, fmt.format(M[yi, xi]), ha="center", va="center",
                    fontsize=9, color="black")
    ax.set_title(title, fontsize=11, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.046)
    cb.set_label(cbar_label, fontsize=8.5)


heat(axs[0], opt,
     "a  Optimal myCAF level per regime\n>0 = preserve/modulate stroma (resource);  0 = remove stroma",
     "YlGn", "{:.2f}", "optimal k_caf_activate")
heat(axs[1], benefit,
     "b  Benefit of keeping stroma\ntumor(no stroma) - tumor(optimal stroma);  >0 = stroma aids control",
     "RdYlGn", "{:+.2f}", "tumor-fold reduction")

fig.suptitle("Phase map — the myCAF-as-controllable-resource hypothesis holds where "
             "confinement outweighs immunosuppression\n(lower-right region); elsewhere stromal "
             "reduction is favored (in-silico)", fontsize=11.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(ROOT, "assets", "phase_map.png")
fig.savefig(out, dpi=120, bbox_inches="tight")
print("wrote", out)
# 요약
res = (opt > 0).sum()
print(f"'기질=자원' 성립 셀: {res}/{opt.size}  (최적 myCAF>0)")
