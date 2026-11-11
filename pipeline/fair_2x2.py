"""L4 공정 2×2 — 스케줄 효과와 약물 효과를 분리 (Fig. S5).

리뷰어 지적: 자연물=적응형 vs 젬=연속 비교는 불공정. 이 스크립트는
  {gemcitabine, natural combo} × {continuous, adaptive}
4팔을 같은 체제·기간에 돌려 주효과(스케줄/약물)와 상호작용을 분리한다.
자연물 조합의 사전주입 이점을 없애기 위해 synergy=0.
공정 비교는 (a) 동일 스케줄 내 약물 비교, (b) 동일 약물 내 스케줄 비교로 읽는다.

300 seeds (원래 5 seeds — 5시드 평균은 적응형 잔여 종양을 과소추정했다).
조직은 고정(make_tissue seed=42), 시뮬레이션 seed 만 1..300.
출력: data/fair_2x2_seeds.csv (팔 × seed 한 줄), assets/fair_2x2.png
  python pipeline/fair_2x2.py           # 1,200 runs (병렬) → CSV + 그림
  python pipeline/fair_2x2.py replot    # 저장된 CSV 로 그림만
"""
import csv
import math
import os
import sys
import warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "data", "fair_2x2_seeds.csv")
PARAMS = dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10,
              init_resistant_frac=0.03, mutation_rate=0.003,
              resistant_immune_evasion=0.35)
DAYS = 150
ADAPT_ON, ADAPT_OFF = 1.1, 0.7
SEEDS = list(range(1, 301))
WORKERS = 6
GEM = [("gemcitabine", 1.0)]
NAT = [("curcumin", 1.0), ("garlic", 1.0), ("ginsenoside_rg3", 1.0)]
AGENTS = {"Gemcitabine": GEM, "Natural combo": NAT}
SCHEDULES = ["continuous", "adaptive"]
ARMS = [(a, s) for a in AGENTS for s in SCHEDULES]
FIELDS = ["agent", "schedule", "seed", "final_frac", "cum_toxicity",
          "final_resistant_frac", "peak_frac", "auc_burden", "ttp_days",
          "progression_censored"]
_T = {}


def _tissue():
    if not _T:
        from synthetic import make_tissue
        cells, labels, _ = make_tissue("contained", seed=42)
        _T.update(c=cells, l=labels, n0=int((labels == "Tumor").sum()))
    return _T


def run(task):
    agent, sched, seed = task
    from abm import simulate, control_metrics
    t = _tissue()
    p = dict(PARAMS)
    p["seed"] = seed
    kw = dict(days=DAYS, params=p, regimen_subs=AGENTS[agent], synergy=0.0,
              snapshots=(0.0, 1.0))
    if sched == "adaptive":
        kw.update(adaptive=True, adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF)
    h, _ = simulate(t["c"], t["l"], **kw)
    m = control_metrics(h, n0=t["n0"])
    row = {"agent": agent, "schedule": sched, "seed": seed}
    row.update({k: float(m[k]) for k in FIELDS[3:]})
    return row


def load():
    with open(CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r["seed"] = int(r["seed"])
        for k in FIELDS[3:]:
            r[k] = float(r[k])
    return rows


def _ci(xs):
    xs = np.asarray(xs, float)
    return xs.mean(), 1.96 * xs.std(ddof=1) / math.sqrt(len(xs))


def summarize(rows):
    by = {arm: [r for r in rows if (r["agent"], r["schedule"]) == arm] for arm in ARMS}
    print(f"{'arm':<28}{'n':>4}{'final (mean±95%CI)':>22}{'median':>8}{'≥0.5x':>7}"
          f"{'expo':>8}{'resist':>8}{'progressed':>11}")
    for arm, rs in by.items():
        ff = [r["final_frac"] for r in rs]
        m, h = _ci(ff)
        prog = sum(r["progression_censored"] < 0.5 for r in rs)
        print(f"{arm[0]+' / '+arm[1]:<28}{len(rs):>4}{m:>14.3f} ±{h:>6.3f}"
              f"{np.median(ff):>8.2f}{sum(x >= 0.5 for x in ff):>7}"
              f"{np.mean([r['cum_toxicity'] for r in rs]):>8.1f}"
              f"{np.mean([r['final_resistant_frac'] for r in rs]):>8.3f}"
              f"{prog:>8}/{len(rs)}")
    # 같은 seed 끼리 짝지은 주효과 (final tumor 기준)
    val = {(r["agent"], r["schedule"], r["seed"]): r["final_frac"] for r in rows}
    seeds = sorted({r["seed"] for r in rows})

    def d(a, b):
        return [val[a + (s,)] - val[b + (s,)] for s in seeds]
    sched = [0.5 * (x + y) for x, y in zip(d(("Gemcitabine", "adaptive"), ("Gemcitabine", "continuous")),
                                           d(("Natural combo", "adaptive"), ("Natural combo", "continuous")))]
    agent = [0.5 * (x + y) for x, y in zip(d(("Natural combo", "continuous"), ("Gemcitabine", "continuous")),
                                           d(("Natural combo", "adaptive"), ("Gemcitabine", "adaptive")))]
    inter = [x - y for x, y in zip(d(("Natural combo", "adaptive"), ("Natural combo", "continuous")),
                                   d(("Gemcitabine", "adaptive"), ("Gemcitabine", "continuous")))]
    for name, xs in (("스케줄(적응−연속)", sched), ("약물(자연−젬)", agent), ("상호작용", inter)):
        m, h = _ci(xs)
        print(f"[주효과·짝] {name}: {m:+.3f}x ±{h:.3f}")
    return by


def main():
    from concurrent.futures import ProcessPoolExecutor
    tasks = [(a, s, seed) for seed in SEEDS for a, s in ARMS]
    print(f"{len(ARMS)} arms × {len(SEEDS)} seeds = {len(tasks)} runs, {WORKERS}워커 병렬",
          flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for i, r in enumerate(ex.map(run, tasks, chunksize=4)):
            rows.append(r)
            if (i + 1) % 100 == 0 or i == len(tasks) - 1:
                print(f"  {i+1}/{len(tasks)}", flush=True)
    rows.sort(key=lambda r: (r["seed"], ARMS.index((r["agent"], r["schedule"]))))
    with open(CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print("wrote", CSV)
    _figure(summarize(rows))


def replot():
    _figure(summarize(load()))


def _figure(by):
    import matplotlib
    matplotlib.use("Agg")
    import fonts  # noqa
    import matplotlib.pyplot as plt
    matplotlib.rcParams["axes.unicode_minus"] = False
    n = len(next(iter(by.values())))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    labels = ["Gem\ncontinuous", "Gem\nadaptive", "Natural\ncontinuous",
              "Natural\nadaptive"]
    cols = ["#C0392B", "#E59866", "#1E8449", "#7DCEA0"]
    x = np.arange(4)
    rng = np.random.default_rng(0)
    for i, arm in enumerate(ARMS):
        ff = np.array([r["final_frac"] for r in by[arm]])
        m, h = _ci(ff)
        a1.bar(i, m, color=cols[i], edgecolor="black", linewidth=0.7, alpha=0.85, zorder=1)
        a1.scatter(i + rng.uniform(-0.28, 0.28, len(ff)), ff, s=5, color="black",
                   alpha=0.25, linewidths=0, zorder=2)
        a1.errorbar(i, m, yerr=h, color="black", capsize=5, lw=1.4, zorder=3)
        a1.text(i, max(ff.max(), m + h) + 0.03,
                f"mean {m:.2f}x\nmedian {np.median(ff):.2f}x", ha="center", fontsize=8.5)
        ex = np.array([r["cum_toxicity"] for r in by[arm]])
        em, eh = _ci(ex)
        a2.bar(i, em, yerr=eh, capsize=5, color=cols[i], edgecolor="black", linewidth=0.7)
        a2.text(i, em + eh + 1.5, f"{em:.0f}", ha="center", fontsize=9)
    a1.axhline(1.5, color="grey", ls="--", lw=0.9)
    a1.text(3.45, 1.5, "progression\n(1.5x)", va="center", ha="right", fontsize=8, color="grey")
    a1.set_ylim(0, 1.75)
    a1.set_xticks(x)
    a1.set_xticklabels(labels, fontsize=9)
    a1.set_ylabel("Final tumor burden (fold vs initial)")
    a1.set_title("a  Tumor control (lower = better; dots = seeds)", fontsize=11,
                 fontweight="bold")
    a2.set_xticks(x)
    a2.set_xticklabels(labels, fontsize=9)
    a2.set_ylabel("Cumulative exposure (patient burden)")
    a2.set_title("b  Treatment exposure (lower = better)", fontsize=11, fontweight="bold")
    fig.suptitle(f"Fair 2×2 — agent (gemcitabine vs natural) × schedule (continuous vs "
                 f"adaptive), synergy off, {n} seeds (mean ± 95% CI)\nadaptive lowers "
                 "exposure for both agents at the cost of residual tumor; no run progressed",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "assets", "fair_2x2.png")
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
