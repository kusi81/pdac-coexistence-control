"""L4 공정 2×2 — 스케줄 효과와 약물 효과를 분리.

리뷰어 지적: 자연물=적응형 vs 젬=연속 비교는 불공정. 이 스크립트는
  {gemcitabine, natural combo} × {continuous, adaptive}
4팔을 같은 체제·기간에 돌려 주효과(스케줄/약물)와 상호작용을 분리한다.
자연물 조합의 사전주입 이점을 없애기 위해 synergy=0.
공정 비교는 (a) 동일 스케줄 내 약물 비교, (b) 동일 약물 내 스케줄 비교로 읽는다.
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
PARAMS = dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10,
              init_resistant_frac=0.03, mutation_rate=0.003,
              resistant_immune_evasion=0.35)
DAYS = 150
ADAPT_ON, ADAPT_OFF = 1.1, 0.7
SEEDS = [1, 2, 3, 4, 5]
GEM = [("gemcitabine", 1.0)]
NAT = [("curcumin", 1.0), ("garlic", 1.0), ("ginsenoside_rg3", 1.0)]
AGENTS = {"Gemcitabine": GEM, "Natural combo": NAT}
SCHEDULES = ["continuous", "adaptive"]


def run(subs, sched, seed):
    p = dict(PARAMS); p["seed"] = seed
    kw = dict(days=DAYS, params=p, regimen_subs=subs, synergy=0.0,
              snapshots=(0.0, 1.0))
    if sched == "adaptive":
        kw.update(adaptive=True, adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF)
    h, _ = simulate(c, l, **kw)
    return control_metrics(h, n0=n0)


c, l, _ = make_tissue("contained", seed=42)
n0 = (l == "Tumor").sum()


def main():
    # cell[(agent,sched)] = dict of mean metrics
    cell = {}
    print(f"{'arm':<32}{'final':>7}{'expo':>7}{'resist':>8}{'ctrl':>7}")
    for ag, subs in AGENTS.items():
        for sc in SCHEDULES:
            ms = [run(subs, sc, s) for s in SEEDS]
            agg = {k: float(np.mean([m[k] for m in ms]))
                   for k in ("final_frac", "cum_toxicity",
                             "final_resistant_frac", "control_score")}
            cell[(ag, sc)] = agg
            print(f"{ag+' / '+sc:<32}{agg['final_frac']:>6.2f}x"
                  f"{agg['cum_toxicity']:>7.0f}{agg['final_resistant_frac']:>8.2f}"
                  f"{agg['control_score']:>7.1f}")

    # 주효과 (final tumor 기준)
    def ff(a, s):
        return cell[(a, s)]["final_frac"]
    sched_eff = 0.5 * ((ff("Gemcitabine", "adaptive") - ff("Gemcitabine", "continuous"))
                       + (ff("Natural combo", "adaptive") - ff("Natural combo", "continuous")))
    agent_eff = 0.5 * ((ff("Natural combo", "continuous") - ff("Gemcitabine", "continuous"))
                       + (ff("Natural combo", "adaptive") - ff("Gemcitabine", "adaptive")))
    inter = (ff("Natural combo", "adaptive") - ff("Natural combo", "continuous")) \
        - (ff("Gemcitabine", "adaptive") - ff("Gemcitabine", "continuous"))
    print(f"\n[주효과] 스케줄(적응−연속): {sched_eff:+.2f}x  "
          f"약물(자연−젬): {agent_eff:+.2f}x  상호작용: {inter:+.2f}x")

    _figure(cell)
    return cell


def _figure(cell):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    arms = [("Gemcitabine", "continuous"), ("Gemcitabine", "adaptive"),
            ("Natural combo", "continuous"), ("Natural combo", "adaptive")]
    labels = ["Gem\ncontinuous", "Gem\nadaptive", "Natural\ncontinuous",
              "Natural\nadaptive"]
    cols = ["#C0392B", "#E59866", "#1E8449", "#7DCEA0"]
    x = np.arange(4)
    fin = [cell[a]["final_frac"] for a in arms]
    expo = [cell[a]["cum_toxicity"] for a in arms]
    a1.bar(x, fin, color=cols, edgecolor="black", linewidth=0.7)
    for i, v in enumerate(fin):
        a1.text(i, v + 0.01, f"{v:.2f}x", ha="center", fontsize=9)
    a1.set_xticks(x); a1.set_xticklabels(labels, fontsize=9)
    a1.set_ylabel("Final tumor burden (fold vs initial)")
    a1.set_title("a  Tumor control (lower = better)", fontsize=11, fontweight="bold")
    a2.bar(x, expo, color=cols, edgecolor="black", linewidth=0.7)
    for i, v in enumerate(expo):
        a2.text(i, v + 1, f"{v:.0f}", ha="center", fontsize=9)
    a2.set_xticks(x); a2.set_xticklabels(labels, fontsize=9)
    a2.set_ylabel("Cumulative exposure (patient burden)")
    a2.set_title("b  Treatment exposure (lower = better)", fontsize=11, fontweight="bold")
    fig.suptitle("Fair 2×2 — agent (gemcitabine vs natural) × schedule (continuous vs "
                 "adaptive), synergy off\nadaptive lowers exposure for both agents; "
                 "compare within-schedule to isolate the agent effect", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(ROOT, "assets", "fair_2x2.png")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
