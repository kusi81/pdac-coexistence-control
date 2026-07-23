"""단일 encoded agent(gemcitabine)의 스케줄 비교 — 리뷰어/사용자 지적.

지적: continuous maximum-dose gemcitabine은 임상 표준이 아니며, 자연물 우월성으로
오해될 수 있다. 핵심 메시지는 '어떤 agent가 낫다'가 아니라:
  동일 encoded agent·동일 효능계수에서도, adaptive scheduling이 continuous full-
  intensity보다 적은 modeled exposure로 control을 유지할 수 있다.

같은 gemcitabine을 스케줄만 바꿔 비교(효능계수 고정):
  1) Continuous full-intensity        — 매일 투여(모델 참조 arm; 임상 표준 아님)
  2) Intermittent (clinical approx.)   — 3주 투여/1주 휴약 q28(임상 근사 간헐)
  3) Ideal adaptive                    — 부담 band 촉발(이상적 관측)
  4) Observed CA19-9 adaptive          — 지연·잡음·이산측정 관측(현실적)

20 seed, 6워커 병렬. 결과: data/gem_schedule.csv, assets/gem_schedule.png.
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
PARAMS = dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10,
              init_resistant_frac=0.03, mutation_rate=0.003,
              resistant_immune_evasion=0.35, max_tumor=6000)
GEM = [("gemcitabine", 1.0)]            # 동일 agent·동일 효능계수(스케줄만 변경)
DAYS = 150
ADAPT_ON, ADAPT_OFF = 1.1, 0.7
CAP = 1800
SEEDS = list(range(1, 21))
WORKERS = 6
OBS = dict(obs_model=True, obs_interval=28.0, obs_noise_cv=0.25,
           obs_lag_days=14.0, min_on_days=14.0, min_off_days=21.0,
           obs_confirm=2, obs_safety_mult=1.4)
# 임상 근사 간헐: 3주 투여 / 1주 휴약 (q28)
INTERMITTENT = [dict(label="on", days=21, substances=GEM),
                dict(label="off", days=7, substances=[])]
ARMS = ["Continuous", "Intermittent q28", "Ideal adaptive", "Observed CA19-9"]


def _tissue(seed):
    c, l, _ = make_tissue("contained", seed=seed)
    idx = np.random.default_rng(1000 + seed).choice(len(l), CAP, replace=False)
    return c[idx], l[idx], int((l[idx] == "Tumor").sum())


def _run(seed, ai):
    c, l, n0 = _tissue(seed)
    p = dict(PARAMS); p["seed"] = seed
    kw = dict(days=DAYS, params=p, snapshots=(0.0, 1.0))
    if ai == 0:                                   # continuous full-intensity
        kw["regimen_subs"] = GEM
    elif ai == 1:                                 # intermittent clinical q28
        kw["schedule"] = INTERMITTENT
    elif ai == 2:                                 # ideal adaptive
        kw.update(regimen_subs=GEM, adaptive=True, adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF)
    else:                                         # observed CA19-9 adaptive
        p.update(OBS)
        kw.update(regimen_subs=GEM, adaptive=True, adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF)
    h, _ = simulate(c, l, **kw)
    return h, control_metrics(h, n0=n0)


def evaluate(task):
    ai, seed = task
    _, m = _run(seed, ai)
    return (ai, seed, m["cum_toxicity"], float(m["progression_censored"]),
            m["final_frac"])


def main():
    from concurrent.futures import ProcessPoolExecutor
    tasks = [(ai, s) for ai in range(len(ARMS)) for s in SEEDS]
    print(f"{len(tasks)} runs, {WORKERS}워커 병렬", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for k, r in enumerate(ex.map(evaluate, tasks, chunksize=4)):
            rows.append(r)
            if (k + 1) % 40 == 0 or k == len(tasks) - 1:
                print(f"  {k+1}/{len(tasks)}", flush=True)
    df = pd.DataFrame(rows, columns=["ai", "seed", "exposure", "censored", "final"])
    df.to_csv(os.path.join(ROOT, "data", "gem_schedule.csv"), index=False)
    print(f"\n{'schedule (same agent)':<24}{'expo med[IQR]':>18}{'control%':>10}{'final med':>11}")
    for ai, label in enumerate(ARMS):
        g = df[df.ai == ai]
        lo, md, hi = np.percentile(g.exposure, [25, 50, 75])
        print(f"{label:<24}{md:>7.0f}[{lo:.0f}-{hi:.0f}]{g.censored.mean()*100:>9.0f}%"
              f"{np.median(g['final']):>10.2f}x")
    figure(df)


def figure(df):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2))
    cols = ["#C0392B", "#E67E22", "#27AE60", "#2980B9"]
    # a: exposure boxplots
    data = [df[df.ai == ai]["exposure"].values for ai in range(len(ARMS))]
    bp = a1.boxplot(data, patch_artist=True, widths=0.6,
                    medianprops=dict(color="black", lw=1.4))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    for ai in range(len(ARMS)):
        g = df[df.ai == ai]
        a1.text(ai + 1, np.median(g.exposure) + 4,
                f"ctrl {g.censored.mean()*100:.0f}%", ha="center", fontsize=8)
    a1.set_xticklabels(["Continuous\nfull-intensity", "Intermittent\nq28 (clinical)",
                        "Ideal\nadaptive", "Observed\nCA19-9"], fontsize=8)
    a1.set_ylabel("Cumulative modeled exposure  →  lower better")
    a1.set_title("a  Same encoded agent (gemcitabine), schedule only\n"
                 "adaptive holds control at lower exposure than fixed schedules",
                 fontsize=10, fontweight="bold")
    # b: final burden boxplots (control quality)
    dataf = [df[df.ai == ai]["final"].values for ai in range(len(ARMS))]
    bp2 = a2.boxplot(dataf, patch_artist=True, widths=0.6,
                     medianprops=dict(color="black", lw=1.4))
    for patch, c in zip(bp2["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    a2.axhline(1.5, color="#E74C3C", ls=":", lw=1)
    a2.text(len(ARMS), 1.53, "progression 1.5×", ha="right", fontsize=7.5,
            color="#E74C3C")
    a2.set_xticklabels(["Continuous", "Intermittent", "Ideal\nadaptive",
                        "Observed\nCA19-9"], fontsize=8)
    a2.set_ylabel("Final tumor burden (fold)  →  lower better")
    a2.set_title("b  Tumor control is maintained across schedules\n"
                 "(exposure differs; control does not)", fontsize=10,
                 fontweight="bold")
    fig.suptitle("Schedule comparison for one encoded agent (gemcitabine, fixed efficacy "
                 "coefficient): adaptive scheduling maintains control at lower modeled "
                 "exposure than continuous full-intensity dosing — in silico",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "assets", "gem_schedule.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    figure(pd.read_csv(os.path.join(ROOT, "data", "gem_schedule.csv")))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
