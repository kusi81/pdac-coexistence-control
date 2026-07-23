"""다중 seed 재계산 + Pareto frontier — 리뷰어 통계 지적 방어.

리뷰어 지적:
  (1) 3-seed 평균으로 regimen ranking을 주장하기엔 확률성 평가가 부족 →
      20~50 seed, median·95% 구간, rank stability.
  (2) 단일 control_score는 무처치(독성0)를 과대평가 → composite 대신
      progression-vs-exposure Pareto frontier로 표시.

설계: 각 regimen을 SEEDS개 seed로 (조직 seed + 시뮬 seed 동시 변동) 재계산.
      조합의 synergy=0(공정). 폭주 방지 max_tumor 상한. 6워커 병렬.
출력: data/pareto_seeds.csv (regimen별 median·95%CI·censored_frac·rank 분포),
      assets/pareto_seeds.png (a: Pareto frontier w/ 95% CI, b: rank stability).
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
DAYS = 150
ADAPT_ON, ADAPT_OFF = 1.1, 0.7
CAP = 2500                 # 조직 서브샘플(성능; 구조 보존)
SEEDS = list(range(1, 31))  # 30 seeds
WORKERS = 6

# regimen 목록: (표시명, 유형, 서브스턴스, adaptive?)
REGIMENS = [
    ("Untreated", "none", None, False),
    ("Curcumin", "natural", [("curcumin", 1.0)], True),
    ("Garlic", "natural", [("garlic", 1.0)], True),
    ("Mugwort", "natural", [("mugwort", 1.0)], True),
    ("Wild ginseng", "natural", [("wild_ginseng", 1.0)], True),
    ("Garlic+Mugwort", "natural", [("garlic", 1.0), ("mugwort", 1.0)], True),
    ("Gin+Gar+Mug", "natural",
     [("wild_ginseng", 1.0), ("garlic", 1.0), ("mugwort", 1.0)], True),
    ("Cur+Gar+Rg3", "natural",
     [("curcumin", 1.0), ("garlic", 1.0), ("ginsenoside_rg3", 1.0)], True),
    ("Danshen+Astra", "natural",
     [("danshen", 1.0), ("astragaloside", 1.0)], True),
    ("Gemcitabine (cont.)", "conventional", [("gemcitabine", 1.0)], False),
]
COL = {"none": "#7F8C8D", "natural": "#27AE60", "conventional": "#C0392B"}
MARK = {"none": "^", "natural": "o", "conventional": "s"}


def _tissue(seed):
    c, l, _ = make_tissue("contained", seed=seed)
    idx = np.random.default_rng(1000 + seed).choice(len(l), CAP, replace=False)
    return c[idx], l[idx]


def evaluate(task):
    """task=(reg_index, seed) → metrics dict. 모듈 최상위(병렬 picklable)."""
    ri, seed = task
    name, typ, subs, adaptive = REGIMENS[ri]
    c, l = _tissue(seed)
    n0 = int((l == "Tumor").sum())
    p = dict(PARAMS); p["seed"] = seed
    kw = dict(days=DAYS, params=p, snapshots=(0.0, 1.0))
    if subs is not None:
        kw["regimen_subs"] = subs
        kw["synergy"] = 0.0                 # 공정: 조합 시너지 off
    if adaptive:
        kw.update(adaptive=True, adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF)
    h, _ = simulate(c, l, **kw)
    m = control_metrics(h, n0=n0, crit_mult=1.5)
    return dict(ri=ri, seed=seed,
                censored=bool(m["progression_censored"]),
                ttp=float(m["ttp_days"]),
                exposure=float(m["cum_toxicity"]),
                final=float(m["final_frac"]),
                resistant=float(m["final_resistant_frac"]))


def main():
    from concurrent.futures import ProcessPoolExecutor
    tasks = [(ri, s) for ri in range(len(REGIMENS)) for s in SEEDS]
    print(f"{len(REGIMENS)} regimens × {len(SEEDS)} seeds = {len(tasks)} runs, "
          f"{WORKERS}워커 병렬", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for i, r in enumerate(ex.map(evaluate, tasks, chunksize=4)):
            rows.append(r)
            if (i + 1) % 60 == 0 or i == len(tasks) - 1:
                print(f"  {i+1}/{len(tasks)}", flush=True)
    df = pd.DataFrame(rows)

    # ── rank stability: 각 seed 내에서 (censored 우선, 그다음 exposure 낮은 순) 랭크 ──
    df["rank"] = np.nan
    for s in SEEDS:
        sub = df[df.seed == s].copy()
        sub = sub.sort_values(["censored", "exposure"],
                              ascending=[False, True])   # censored=True 먼저
        sub["rk"] = np.arange(1, len(sub) + 1)
        df.loc[sub.index, "rank"] = sub["rk"].values

    # ── regimen별 집계 ──
    agg = []
    for ri, (name, typ, subs, adaptive) in enumerate(REGIMENS):
        g = df[df.ri == ri]
        ex_lo, ex_md, ex_hi = np.percentile(g.exposure, [2.5, 50, 97.5])
        fn_lo, fn_md, fn_hi = np.percentile(g.final, [2.5, 50, 97.5])
        agg.append(dict(
            regimen=name, type=typ,
            censored_frac=float(g.censored.mean()),
            exposure_med=ex_md, exposure_lo=ex_lo, exposure_hi=ex_hi,
            final_med=fn_md, final_lo=fn_lo, final_hi=fn_hi,
            resistant_med=float(np.median(g.resistant)),
            rank_med=float(np.median(g["rank"])),
            rank_lo=float(np.percentile(g["rank"], 2.5)),
            rank_hi=float(np.percentile(g["rank"], 97.5)),
            top3_frac=float((g["rank"] <= 3).mean())))
    A = pd.DataFrame(agg)
    os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
    A.to_csv(os.path.join(ROOT, "data", "pareto_seeds.csv"), index=False)
    df.to_csv(os.path.join(ROOT, "data", "pareto_seeds_raw.csv"), index=False)

    print(f"\n{'regimen':<20}{'cens%':>6}{'expo(med[95%])':>20}"
          f"{'final(med)':>11}{'rank med[95%]':>15}{'top3%':>7}")
    for _, r in A.sort_values(["censored_frac", "exposure_med"],
                              ascending=[False, True]).iterrows():
        print(f"{r['regimen']:<20}{r['censored_frac']*100:>5.0f}%"
              f"{r['exposure_med']:>8.0f}[{r['exposure_lo']:.0f}-{r['exposure_hi']:.0f}]"
              f"{r['final_med']:>10.2f}x"
              f"{r['rank_med']:>7.0f}[{r['rank_lo']:.0f}-{r['rank_hi']:.0f}]"
              f"{r['top3_frac']*100:>6.0f}%")

    figure(A, df)


def _pareto_front(pts):
    """pts=[(x,y,idx)] lower-left 우세. 비지배 점 idx 집합 반환."""
    front = []
    for x, y, idx in pts:
        dominated = any((x2 <= x and y2 <= y and (x2 < x or y2 < y))
                        for x2, y2, _ in pts)
        if not dominated:
            front.append(idx)
    return set(front)


def _plot_points(ax, A, fs=7.2, dx=5, dy=4, label=True, stagger=False):
    Asrt = A.sort_values("exposure_med") if stagger else A
    tiers = [26, 40, 54, 12]           # staggered label heights (offset pts)
    for k, (i, r) in enumerate(Asrt.iterrows()):
        c = COL[r["type"]]; m = MARK[r["type"]]
        ax.errorbar(r["exposure_med"], r["final_med"],
                    xerr=[[max(0, r["exposure_med"] - r["exposure_lo"])],
                          [max(0, r["exposure_hi"] - r["exposure_med"])]],
                    yerr=[[max(0, r["final_med"] - r["final_lo"])],
                          [max(0, r["final_hi"] - r["final_med"])]],
                    fmt=m, color=c, ms=9, mec="black", mew=0.8,
                    ecolor=c, elinewidth=0.9, capsize=2, alpha=0.9, zorder=3)
        if not label:
            continue
        if stagger:
            yo = tiers[k % len(tiers)]
            ax.annotate(r["regimen"], (r["exposure_med"], r["final_med"]),
                        xytext=(0, yo), textcoords="offset points", fontsize=fs,
                        ha="center", color="#1E6B3A", rotation=0,
                        arrowprops=dict(arrowstyle="-", color="#BBB", lw=0.6))
        else:
            ax.annotate(r["regimen"], (r["exposure_med"], r["final_med"]),
                        xytext=(dx, dy), textcoords="offset points", fontsize=fs,
                        color="#333")


def figure(A, df):
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14.5, 6.0))

    # ── Panel a: Pareto frontier (exposure vs final burden), 95% CI ──
    pts = [(r["exposure_med"], r["final_med"], i) for i, r in A.iterrows()]
    front = _pareto_front(pts)
    _plot_points(a1, A, label=False)
    # gemcitabine + untreated + danshen labeled in main axes; cluster labeled in inset
    for i, r in A.iterrows():
        if r["exposure_med"] > 12 or r["final_med"] > 0.9:
            a1.annotate(r["regimen"], (r["exposure_med"], r["final_med"]),
                        xytext=(6, 4), textcoords="offset points", fontsize=7.8,
                        color="#333")
    fx = sorted([(A.loc[i, "exposure_med"], A.loc[i, "final_med"]) for i in front])
    a1.plot([p[0] for p in fx], [p[1] for p in fx], "--", color="#555",
            lw=1.2, zorder=1, label="Pareto frontier")
    a1.axhline(1.5, color="#E74C3C", ls=":", lw=1)
    a1.text(a1.get_xlim()[1] * 0.98, 1.57, "progression 1.5×", ha="right",
            fontsize=7.5, color="#E74C3C")
    a1.set_xlabel("Cumulative modeled exposure  →  lower better")
    a1.set_ylabel("Final tumor burden (fold)  →  lower better")
    a1.set_title("a  Progression-vs-exposure Pareto frontier\n"
                 "median ± 95% interval over 30 seeds (synergy off)",
                 fontsize=10.5, fontweight="bold")
    a1.legend(frameon=False, fontsize=8, loc="upper center")

    # inset: zoom the crowded low-exposure controlling cluster (exposure 0-16)
    axin = inset_axes(a1, width="55%", height="52%", loc="center right",
                      borderpad=1.2)
    _plot_points(axin, A[A.exposure_med <= 16], fs=7.0, stagger=True)
    axin.set_xlim(-1, 16); axin.set_ylim(-0.05, 0.62)
    axin.tick_params(labelsize=7)
    axin.set_title("low-exposure controllers (zoom)", fontsize=7.5)
    for sp in axin.spines.values():
        sp.set_edgecolor("#999")
    mark_inset(a1, axin, loc1=2, loc2=3, fc="none", ec="#BBB", lw=0.7)

    # ── Panel b: rank stability (per-seed rank distribution) ──
    order = A.sort_values(["rank_med", "exposure_med"]).index.tolist()
    data = [df[df.ri == i]["rank"].values for i in order]
    labels = [A.loc[i, "regimen"] for i in order]
    cols = [COL[A.loc[i, "type"]] for i in order]
    bp = a2.boxplot(data, vert=False, patch_artist=True, widths=0.6,
                    medianprops=dict(color="black", lw=1.4),
                    flierprops=dict(marker=".", ms=3, alpha=0.4))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.55)
    a2.set_yticklabels(labels, fontsize=8)
    a2.invert_yaxis()
    a2.set_xlabel("Per-seed rank (1 = best: controlled at lowest exposure)")
    a2.set_title("b  Rank stability across 30 seeds\n"
                 "(rank by progression control, then exposure)",
                 fontsize=10.5, fontweight="bold")

    fig.suptitle("Multi-seed re-analysis (n=30): progression-vs-exposure Pareto frontier "
                 "and ranking stability — in silico, under assigned exposure weights",
                 fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(ROOT, "assets", "pareto_seeds.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    A = pd.read_csv(os.path.join(ROOT, "data", "pareto_seeds.csv"))
    df = pd.read_csv(os.path.join(ROOT, "data", "pareto_seeds_raw.csv"))
    figure(A, df)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
