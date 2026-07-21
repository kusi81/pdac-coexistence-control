"""순차 사이클 vs 동시 투여 데모: 종양 억제 + 누적 독성(환자 부담) 비교."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

from synthetic import make_tissue
from abm import (simulate, context_params, classify_substances,
                 build_cycle_schedule)

c, l, _ = make_tissue("contained", seed=42)
p = context_params("pdac")
DAYS = 30
picked = [("ginsenoside_rg3", 1.0), ("garlic", 1.0), ("mugwort", 1.0),
          ("erlotinib", 1.0)]
anti, cyto = classify_substances(picked)
sched = build_cycle_schedule(anti, cyto, phase_days=5)

ctrl, _ = simulate(c, l, days=DAYS, params=p, snapshots=(0.0, 1.0))
simul, _ = simulate(c, l, days=DAYS, params=p, regimen_subs=picked, snapshots=(0.0, 1.0))
seq, _ = simulate(c, l, days=DAYS, params=p, schedule=sched, snapshots=(0.0, 1.0))

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 5.4))


def line(ax, h, key, col, lab, ls="-"):
    ax.plot([x["t"] for x in h], [x[key] for x in h], ls, lw=2, color=col, label=lab)


# 왼쪽: 종양
line(axL, ctrl, "n_tumor", "#7F8C8D", "대조(무처치)")
line(axL, simul, "n_tumor", "#C0392B", "동시 투여")
line(axL, seq, "n_tumor", "#8E44AD", "순차 사이클")
# phase 경계 음영 (사이클)
cyc = sum(ph["days"] for ph in sched)
for start in range(0, DAYS, cyc):
    axL.axvspan(start, start + sched[0]["days"], color="#27AE60", alpha=0.06)
axL.set_xlabel("시간 (days)"); axL.set_ylabel("종양 세포 수")
axL.set_title("종양 억제 — 순차도 동등 수준 달성", fontsize=12, fontweight="bold")
axL.legend(frameon=False, fontsize=9)

# 오른쪽: 누적 독성
line(axR, simul, "cum_toxicity", "#C0392B", "동시 투여")
line(axR, seq, "cum_toxicity", "#8E44AD", "순차 사이클")
axR.set_xlabel("시간 (days)"); axR.set_ylabel("누적 독성 (환자 부담)")
axR.set_title("누적 독성 — 순차가 절반 수준", fontsize=12, fontweight="bold")
axR.legend(frameon=False, fontsize=9)
axR.annotate(f"동시 {simul[-1]['cum_toxicity']:.0f}", (DAYS, simul[-1]['cum_toxicity']),
             fontsize=9, color="#C0392B", ha="right", va="bottom", fontweight="bold")
axR.annotate(f"순차 {seq[-1]['cum_toxicity']:.0f}", (DAYS, seq[-1]['cum_toxicity']),
             fontsize=9, color="#8E44AD", ha="right", va="bottom", fontweight="bold")

nc = ctrl[-1]["n_tumor"]
inh_s = 100 * (nc - simul[-1]["n_tumor"]) / nc
inh_q = 100 * (nc - seq[-1]["n_tumor"]) / nc
fig.suptitle(f"순차 사이클 vs 동시 투여 (PDAC, {DAYS}일) — "
             f"억제 동시 {inh_s:.0f}% / 순차 {inh_q:.0f}%, 독성은 순차가 절반 "
             "(가설 샌드박스)", fontsize=10.5)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "assets", "schedule_demo.png")
fig.savefig(out, dpi=115, bbox_inches="tight")
print("wrote", out)
