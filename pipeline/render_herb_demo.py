"""전통 약재 시뮬레이션 데모: 종양 궤적 + 억제율(근거 수준 색) 비교."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

from synthetic import make_tissue
from abm import simulate, compose_regimen, SUBSTANCES

coords, labels, _ = make_tissue("contained", seed=42)
DAYS = 20

# (표시명, 조합, 근거수준색)
EVCOL = {"strong": "#27AE60", "moderate": "#F1C40F", "weak": "#E67E22",
         "reference": "#95A5A6"}
arms = [
    ("대조(무처치)", [], "#7F8C8D"),
    ("마늘 [강]", [("garlic", 1.0)], EVCOL["strong"]),
    ("쑥·항CAF [강]", [("mugwort", 1.0)], EVCOL["strong"]),
    ("산삼 [강]", [("wild_ginseng", 1.0)], EVCOL["strong"]),
    ("도라지 [보통]", [("platycodon", 1.0)], EVCOL["moderate"]),
    ("해삼 [보통]", [("sea_cucumber", 1.0)], EVCOL["moderate"]),
    ("녹용 [약함]", [("deer_antler", 1.0)], EVCOL["weak"]),
    ("산삼+마늘+쑥 조합", [("wild_ginseng", 1.0), ("garlic", 1.0),
                     ("mugwort", 1.0)], "#8E44AD"),
]

results = []
for name, sel, col in arms:
    hist, _ = simulate(coords, labels, days=DAYS,
                       regimen=compose_regimen(sel), snapshots=(0.0, 1.0))
    results.append((name, col, hist))

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 6))
n0 = results[0][2][0]["n_tumor"]
nc = results[0][2][-1]["n_tumor"]
for name, col, hist in results:
    t = [h["t"] for h in hist]; y = [h["n_tumor"] for h in hist]
    axL.plot(t, y, "-", lw=2 if "조합" in name else 1.6, color=col, label=name)
axL.axhline(n0, color="#CCC", ls=":", lw=1)
axL.set_xlabel("시간 (days)"); axL.set_ylabel("종양 세포 수")
axL.set_title("종양 성장 궤적 — 전통 약재", fontsize=12, fontweight="bold")
axL.legend(frameon=False, fontsize=8.5)

names = [r[0] for r in results]
inh = [100.0 * (nc - r[2][-1]["n_tumor"]) / nc for r in results]
cols = [r[1] for r in results]
axR.barh(names, inh, color=cols, height=0.62)
axR.axvline(0, color="black", lw=1)
for i, v in enumerate(inh):
    axR.text(v + (0.5 if v >= 0 else -0.5), i, f"{v:+.0f}%", va="center",
             ha="left" if v >= 0 else "right", fontsize=10, fontweight="bold")
axR.set_xlabel("대조 대비 종양 억제율 % (양수=억제)")
axR.set_title("억제율 · 색=근거수준(초록=강 노랑=보통 주황=약)",
              fontsize=12, fontweight="bold")
axR.invert_yaxis()

fig.suptitle("전통 약재 In-silico 시뮬레이션 (가설 샌드박스 — 사람 효과 아님)",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "assets", "herb_demo.png")
fig.savefig(out, dpi=110, bbox_inches="tight")
print("wrote", out)
