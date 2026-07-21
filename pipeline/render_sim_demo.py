"""시뮬레이션 데모 그림: 여러 처치군의 종양 궤적 + 억제율 비교 한 장."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

from synthetic import make_tissue
from abm import simulate, compose_regimen, SUBSTANCES

coords, labels, _ = make_tissue("contained", seed=42)
DAYS = 20

arms = [
    ("대조(무처치)", [], "#7F8C8D"),
    ("항섬유화 단독", [("generic_antifibrotic", 1.0)], "#2E86AB"),
    ("커큐민", [("curcumin", 1.0)], "#E67E22"),
    ("커큐민+Rg3 조합", [("curcumin", 1.0), ("ginsenoside_rg3", 1.0)], "#8E44AD"),
    ("세포독성(젬시타빈)", [("gemcitabine", 1.0)], "#C0392B"),
]

results = []
for name, sel, col in arms:
    hist, _ = simulate(coords, labels, days=DAYS,
                       regimen=compose_regimen(sel), snapshots=(0.0, 1.0))
    results.append((name, sel, col, hist))

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.2))

n0 = results[0][3][0]["n_tumor"]
nc = results[0][3][-1]["n_tumor"]   # 대조 최종
for name, sel, col, hist in results:
    t = [h["t"] for h in hist]
    y = [h["n_tumor"] for h in hist]
    axL.plot(t, y, "-o", ms=3, lw=1.9, color=col, label=name)
axL.axhline(n0, color="#CCC", ls=":", lw=1)
axL.set_xlabel("시간 (days)"); axL.set_ylabel("종양 세포 수")
axL.set_title("종양 성장 궤적 — 처치군 비교", fontsize=12, fontweight="bold")
axL.legend(frameon=False, fontsize=9)

names = [r[0] for r in results]
inh = [100.0 * (nc - r[3][-1]["n_tumor"]) / nc for r in results]
cols = [r[2] for r in results]
axR.barh(names, inh, color=cols, height=0.6)
axR.axvline(0, color="black", lw=1)
for i, v in enumerate(inh):
    axR.text(v + (0.6 if v >= 0 else -0.6), i, f"{v:+.0f}%", va="center",
             ha="left" if v >= 0 else "right", fontsize=10, fontweight="bold")
axR.set_xlabel("대조 대비 종양 억제율 % (양수 = 억제)")
axR.set_title("성장 억제율 — 대조 기준", fontsize=12, fontweight="bold")

fig.suptitle("In-silico 교란 시뮬레이션 (contained 조직, 20일) — 가설 탐색 샌드박스",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "assets", "sim_demo.png")
fig.savefig(out, dpi=110, bbox_inches="tight")
print("wrote", out)
for name, v in zip(names, inh):
    print(f"  {name}: 억제율 {v:+.0f}%")
