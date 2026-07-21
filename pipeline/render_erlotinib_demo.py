"""에를로티닙 병용 데모: 표적 항암제 단독 vs 젬시타빈 병용 시너지."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

from synthetic import make_tissue
from abm import simulate, compose_regimen, context_params

c, l, _ = make_tissue("contained", seed=42)
p = context_params("pdac")
DAYS = 20


def run(sel, syn=0.0):
    h, _ = simulate(c, l, days=DAYS, params=p,
                    regimen=compose_regimen(sel, synergy=syn), snapshots=(0.0, 1.0))
    return h


arms = [
    ("대조(무처치)", [], 0.0, "#7F8C8D"),
    ("젬시타빈 (표준)", [("gemcitabine", 1.0)], 0.0, "#95A5A6"),
    ("에를로티닙 (EGFR 표적)", [("erlotinib", 1.0)], 0.0, "#2E86AB"),
    ("젬시타빈 + 에를로티닙 병용", [("gemcitabine", 1.0), ("erlotinib", 1.0)], 0.5, "#8E44AD"),
    ("엔테카비르 (가설·미검증)", [("entecavir", 1.0)], 0.0, "#C0392B"),
]

results = [(nm, col, run(sel, syn)) for nm, sel, syn, col in arms]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 5.6))
n0 = results[0][2][0]["n_tumor"]
nc = results[0][2][-1]["n_tumor"]
for nm, col, hist in results:
    t = [h["t"] for h in hist]; y = [h["n_tumor"] for h in hist]
    ls = "--" if "가설" in nm else "-"
    axL.plot(t, y, ls, lw=2.2 if "병용" in nm else 1.7, color=col, label=nm)
axL.axhline(n0, color="#CCC", ls=":", lw=1)
axL.set_xlabel("시간 (days)"); axL.set_ylabel("종양 세포 수")
axL.set_title("종양 궤적 — 표적 항암제 & 병용", fontsize=12, fontweight="bold")
axL.legend(frameon=False, fontsize=9)

names = [r[0] for r in results]
inh = [100.0 * (nc - r[2][-1]["n_tumor"]) / nc for r in results]
cols = [r[1] for r in results]
axR.barh(names, inh, color=cols, height=0.6)
axR.axvline(0, color="black", lw=1)
for i, v in enumerate(inh):
    axR.text(v + 0.5, i, f"{v:+.0f}%", va="center", ha="left",
             fontsize=10, fontweight="bold")
axR.set_xlabel("대조 대비 종양 억제율 %")
axR.set_title("병용이 단독을 능가 (실제 임상서도 검증)", fontsize=12, fontweight="bold")
axR.invert_yaxis()

fig.suptitle("EGFR 표적 병용 In-silico (PDAC, 20일) — 빨간 점선=미검증 가설",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "assets", "erlotinib_demo.png")
fig.savefig(out, dpi=115, bbox_inches="tight")
print("wrote", out)
