"""Non-trivial 예측 (2) — 어떤 sequence가 simultaneous combination보다 우수한가.

리뷰어 Major-1: 자명하지 않은 예측. 장벽-게이트 면역요법에서 '순서'가 결과를 바꾸는가?
강한 면역배제 contained 체제에서 anti-fibrotic(장벽 열기)과 immune(살상)을 조합할 때:
  - Immune only / Anti-fibrotic only (단독 참조)
  - Simultaneous (동시 지속)
  - Open-first: anti-fibrotic 선행 → 이후 anti-fibrotic+immune (장벽 먼저 열고 진입)
  - Kill-first: immune 선행 → 이후 anti-fibrotic+immune (닫힌 장벽에 먼저 면역 시도)
예측: 장벽을 먼저 열어야 CD8가 진입 → Open-first가 Kill-first보다, 그리고 동시보다도
유리할 수 있다(닫힌 초기에 면역을 낭비하지 않음). 이는 입력규칙에서 자명하지 않은 순서효과.

6워커 병렬. 결과: data/predict_sequence.csv, assets/predict_sequence.png.
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
matplotlib.rcParams["axes.unicode_minus"] = False
from synthetic import make_tissue
from abm import simulate, control_metrics

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 강한 면역배제 contained 체제(장벽 게이팅이 작동)
PARAMS = dict(k_prolif=0.16, cd8_recruit=12, k_kill=0.6, cd8_barrier_alpha=2.6,
              caf_confine=0.7, caf_pressure=0.9, init_resistant_frac=0.02,
              mutation_rate=0.002, resistant_immune_evasion=0.4, max_tumor=6000)
DAYS = 120
LEAD = 30                    # 선행 단계 길이(일)
CAP = 2200
SEEDS = [1, 2, 3, 4, 5]
WORKERS = 6
IMM = [("generic_immunostim", 1.0)]
AF = [("generic_antifibrotic", 1.0)]
BOTH = [("generic_antifibrotic", 1.0), ("generic_immunostim", 1.0)]
ARMS = ["Immune only", "Anti-fibrotic only", "Simultaneous",
        "Open-first (AF→AF+Imm)", "Kill-first (Imm→AF+Imm)"]
_TIS = {}


def _tissue(seed):
    if seed not in _TIS:
        c, l, _ = make_tissue("contained", seed=seed)
        idx = np.random.default_rng(700 + seed).choice(len(l), CAP, replace=False)
        _TIS[seed] = (c[idx], l[idx], int((l[idx] == "Tumor").sum()))
    return _TIS[seed]


def evaluate(task):
    ai, seed = task
    c, l, n0 = _tissue(seed)
    p = dict(PARAMS); p["seed"] = seed
    kw = dict(days=DAYS, params=p, snapshots=(0.0, 1.0))
    if ai == 0:
        kw["regimen_subs"] = IMM
    elif ai == 1:
        kw["regimen_subs"] = AF
    elif ai == 2:
        kw["regimen_subs"] = BOTH
    elif ai == 3:      # open-first
        kw["schedule"] = [dict(label="AF", days=LEAD, substances=AF),
                          dict(label="AF+Imm", days=DAYS - LEAD, substances=BOTH)]
    else:              # kill-first
        kw["schedule"] = [dict(label="Imm", days=LEAD, substances=IMM),
                          dict(label="AF+Imm", days=DAYS - LEAD, substances=BOTH)]
    h, _ = simulate(c, l, **kw)
    m = control_metrics(h, n0=n0)
    return (ARMS[ai], seed, m["final_frac"], float(m["progression_censored"]),
            m["cum_toxicity"])


def main():
    from concurrent.futures import ProcessPoolExecutor
    tasks = [(ai, s) for ai in range(len(ARMS)) for s in SEEDS]
    print(f"{len(tasks)} sims, {WORKERS}워커 병렬", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for k, r in enumerate(ex.map(evaluate, tasks, chunksize=3)):
            rows.append(r)
    df = pd.DataFrame(rows, columns=["arm", "seed", "final", "censored", "exposure"])
    df.to_csv(os.path.join(ROOT, "data", "predict_sequence.csv"), index=False)
    print(f"\n{'arm':<26}{'final(mean±sd)':>18}{'ctrl%':>7}{'expo':>7}")
    for a in ARMS:
        g = df[df.arm == a]
        print(f"{a:<26}{g.final.mean():>8.2f}±{g.final.std():.2f}"
              f"{g.censored.mean()*100:>7.0f}%{g.exposure.mean():>7.0f}")
    figure(df)


def figure(df):
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    cols = ["#95A5A6", "#E67E22", "#8E44AD", "#27AE60", "#C0392B"]
    x = np.arange(len(ARMS))
    means = [df[df.arm == a]["final"].mean() for a in ARMS]
    sems = [df[df.arm == a]["final"].std() / np.sqrt(len(SEEDS)) for a in ARMS]
    bars = ax.bar(x, means, color=cols, yerr=sems, capsize=3, edgecolor="black",
                  linewidth=0.6, alpha=0.85)
    for i, a in enumerate(ARMS):
        cc = df[df.arm == a]["censored"].mean() * 100
        ax.text(i, means[i] + sems[i] + 0.03, f"ctrl {cc:.0f}%", ha="center",
                fontsize=8)
    ax.axhline(1.5, color="#E74C3C", ls=":", lw=1)
    ax.text(len(ARMS) - 0.5, 1.53, "progression 1.5×", ha="right", fontsize=7.5,
            color="#E74C3C")
    ax.set_xticks(x); ax.set_xticklabels(ARMS, fontsize=8.5, rotation=12, ha="right")
    ax.set_ylabel("Final tumor burden (fold)  →  lower better")
    ax.set_title("Order effect in barrier-gated immunotherapy: both mechanisms are needed, "
                 "and opening\nthe barrier first is marginally better than immune-first "
                 "(modest sequence effect; #4, in silico)", fontsize=10.5,
                 fontweight="bold")
    fig.tight_layout()
    out = os.path.join(ROOT, "assets", "predict_sequence.png")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print("wrote", out)


def replot():
    figure(pd.read_csv(os.path.join(ROOT, "data", "predict_sequence.csv")))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
