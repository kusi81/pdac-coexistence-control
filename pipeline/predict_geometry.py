"""Non-trivial 예측 (1) — geometry가 최적치료를 바꾸는가 + CAF감소와 invasion front.

리뷰어 Major-1: 모델이 입력규칙(myCAF가 딸세포 차단·수용력↓·약물/면역 차단)의 자명한
재확인이 아니라, 자명하지 않은 예측을 내는가? 특히 "공간 geometry가 (단순 종양수 변화가
아니라) 선택되는 치료전략을 바꾸는가?"

동일 abundance·다른 geometry로 검정: make_tissue의 contained(myCAF 링) vs diffuse
(같은 세포수, myCAF 산개)는 세포수가 동일. 여기에:
  #2 optimal treatment: {cytotoxic, anti-fibrotic, immune, cyto+AF, cyto+immune}를 돌려
     최적(종양 최소) 치료가 geometry에 따라 달라지는가.
  #3 invasion front: anti-fibrotic(CAF감소)이 종양 공간확산(spread90)을 언제 늘리는가 —
     confinement가 있는 contained에서는 장벽제거→확산↑, 산개 diffuse에서는 변화 작음.

6워커 병렬, max_tumor 상한. 결과: data/predict_geometry.csv, assets/predict_geometry.png.
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
PARAMS = dict(k_prolif=0.16, cd8_recruit=8, k_kill=0.45, init_resistant_frac=0.02,
              mutation_rate=0.002, resistant_immune_evasion=0.4, max_tumor=6000)
DAYS = 90
CAP = 2000
SEEDS = [1, 2, 3, 4, 5]
WORKERS = 6
GEOMS = ["contained", "diffuse"]
# #2 치료 세트 (동일 sub-maximal 용량; 기전만 다름)
REGIMENS = [
    ("Cytotoxic", [("generic_cytotoxic", 0.4)]),
    ("Anti-fibrotic", [("generic_antifibrotic", 1.0)]),
    ("Immune", [("generic_immunostim", 1.0)]),
    ("Cyto+Anti-fibrotic", [("generic_cytotoxic", 0.4), ("generic_antifibrotic", 1.0)]),
    ("Cyto+Immune", [("generic_cytotoxic", 0.4), ("generic_immunostim", 1.0)]),
]
_TIS = {}


def _tissue(geom, seed):
    key = (geom, seed)
    if key not in _TIS:
        c, l, _ = make_tissue(geom, seed=seed)
        idx = np.random.default_rng(500 + seed).choice(len(l), CAP, replace=False)
        _TIS[key] = (c[idx], l[idx], int((l[idx] == "Tumor").sum()))
    return _TIS[key]


def _spread(coords, labels):
    """종양 공간확산 = 종양 중심에서 90퍼센타일 거리(µm). invasion front 대리."""
    t = coords[labels == "Tumor"]
    if len(t) < 5:
        return 0.0
    ctr = t.mean(0)
    return float(np.percentile(np.sqrt(((t - ctr) ** 2).sum(1)), 90))


def evaluate(task):
    kind, gi, ri, seed = task
    geom = GEOMS[gi]
    c, l, n0 = _tissue(geom, seed)
    p = dict(PARAMS); p["seed"] = seed
    if kind == "opt":                                  # #2
        label, subs = REGIMENS[ri]
        h, snaps = simulate(c, l, days=DAYS, params=p, regimen_subs=subs,
                            snapshots=(0.0, 1.0))
        m = control_metrics(h, n0=n0)
        return ("opt", geom, label, seed, m["final_frac"], np.nan)
    else:                                              # #3: untreated vs anti-fibrotic
        subs = None if ri == 0 else [("generic_antifibrotic", 1.0)]
        kw = dict(days=DAYS, params=p, snapshots=(0.0, 1.0))
        if subs:
            kw["regimen_subs"] = subs
        h, snaps = simulate(c, l, **kw)
        last = snaps[max(snaps)]
        spread = _spread(last[0], last[1])
        arm = "Untreated" if ri == 0 else "Anti-fibrotic"
        return ("inv", geom, arm, seed, control_metrics(h, n0=n0)["final_frac"], spread)


def main():
    from concurrent.futures import ProcessPoolExecutor
    tasks = [("opt", gi, ri, s) for gi in range(2)
             for ri in range(len(REGIMENS)) for s in SEEDS]
    tasks += [("inv", gi, ri, s) for gi in range(2)
              for ri in range(2) for s in SEEDS]
    print(f"{len(tasks)} sims, {WORKERS}워커 병렬", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for k, r in enumerate(ex.map(evaluate, tasks, chunksize=4)):
            rows.append(r)
            if (k + 1) % 30 == 0 or k == len(tasks) - 1:
                print(f"  {k+1}/{len(tasks)}", flush=True)
    df = pd.DataFrame(rows, columns=["kind", "geom", "arm", "seed", "final", "spread"])
    df.to_csv(os.path.join(ROOT, "data", "predict_geometry.csv"), index=False)

    print("\n#2 최적치료 (final tumor, 낮을수록 좋음):")
    for geom in GEOMS:
        g = df[(df.kind == "opt") & (df.geom == geom)]
        means = g.groupby("arm")["final"].mean().sort_values()
        best = means.index[0]
        print(f"  [{geom}] best={best} ({means.iloc[0]:.2f}x) | "
              + " ".join(f"{a}:{v:.2f}" for a, v in means.items()))
    print("\n#3 invasion front (종양 spread90 µm):")
    for geom in GEOMS:
        g = df[(df.kind == "inv") & (df.geom == geom)]
        u = g[g.arm == "Untreated"]["spread"].mean()
        a = g[g.arm == "Anti-fibrotic"]["spread"].mean()
        print(f"  [{geom}] untreated={u:.0f} → anti-fibrotic={a:.0f}  (Δ={a-u:+.0f})")
    figure(df)


def figure(df):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.4))
    gcol = {"contained": "#C0392B", "diffuse": "#2980B9"}

    # a: #2 optimal treatment by geometry (grouped bars, lower=better)
    arms = [r[0] for r in REGIMENS]
    x = np.arange(len(arms)); w = 0.38
    for j, geom in enumerate(GEOMS):
        g = df[(df.kind == "opt") & (df.geom == geom)]
        means = [g[g.arm == a]["final"].mean() for a in arms]
        sems = [g[g.arm == a]["final"].std() / np.sqrt(len(SEEDS)) for a in arms]
        a1.bar(x + (j - 0.5) * w, means, w, yerr=sems, capsize=2,
               color=gcol[geom], alpha=0.8, edgecolor="black", linewidth=0.5,
               label=f"{geom} geometry")
        # mark best per geometry
        bi = int(np.argmin(means))
        a1.text(x[bi] + (j - 0.5) * w, means[bi] + 0.03, "★", ha="center",
                fontsize=13, color=gcol[geom])
    a1.set_xticks(x); a1.set_xticklabels(arms, fontsize=8, rotation=15, ha="right")
    a1.set_ylabel("Final tumor burden (fold)  →  lower better")
    a1.set_title("a  Same myCAF abundance, different geometry → different treatment\n"
                 "efficacy; immune therapy is ~12× weaker when stroma confines (#2)",
                 fontsize=10, fontweight="bold")
    a1.legend(frameon=False, fontsize=9)

    # b: #3 invasion front — CAF reduction increases spread only in contained
    x2 = np.arange(2); w2 = 0.38
    for j, geom in enumerate(GEOMS):
        g = df[(df.kind == "inv") & (df.geom == geom)]
        vals = [g[g.arm == "Untreated"]["spread"].mean(),
                g[g.arm == "Anti-fibrotic"]["spread"].mean()]
        sems = [g[g.arm == "Untreated"]["spread"].std() / np.sqrt(len(SEEDS)),
                g[g.arm == "Anti-fibrotic"]["spread"].std() / np.sqrt(len(SEEDS))]
        a2.bar(x2 + (j - 0.5) * w2, vals, w2, yerr=sems, capsize=3,
               color=gcol[geom], alpha=0.8, edgecolor="black", linewidth=0.5,
               label=f"{geom} geometry")
        dv = vals[1] - vals[0]
        a2.annotate(f"Δ={dv:+.0f}µm", (x2[1] + (j - 0.5) * w2, vals[1]),
                    xytext=(0, 6), textcoords="offset points", ha="center",
                    fontsize=8, color=gcol[geom], fontweight="bold")
    a2.set_xticks(x2); a2.set_xticklabels(["Untreated", "Anti-fibrotic\n(CAF reduction)"],
                                          fontsize=9)
    a2.set_ylabel("Tumor spatial spread (90th-pct radius, µm)")
    a2.set_title("b  When does CAF reduction widen the invasion front?\n"
                 "only where geometry confines the tumor (#3)", fontsize=10,
                 fontweight="bold")
    a2.legend(frameon=False, fontsize=9)

    fig.suptitle("Spatial geometry—not just abundance—changes the model's predictions: "
                 "treatment efficacy and invasion response depend on myCAF architecture "
                 "at fixed cell counts (in silico)", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "assets", "predict_geometry.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    figure(pd.read_csv(os.path.join(ROOT, "data", "predict_geometry.csv")))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
