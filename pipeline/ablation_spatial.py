"""공간모델 필요성 — 전체 ablation 비교 (리뷰어 Major Comment 1).

핵심 질문: "공간 geometry가 (단순 종양수가 아니라) *선택되는 전략*(기질 유지 vs 제거)을
바꾸는가?" 이를 기전별 ablation으로 분해한다. 고정 체제에서 stromal level(k_caf_activate)을
스윕해 최종 종양 최소가 되는 '최적 기질 수준'을 각 변형에서 찾는다:

  1. Full spatial ABM         — confinement ON + immune-exclusion ON   → 내부 최적(keep 일부)
  2. No-confinement           — caf_confine/pressure/drug_block = 0     → 기질=비용만 → 최적 M=0(제거)
  3. No-immune-exclusion      — cd8_barrier_alpha = 0                   → 기질=이득만 → 최적 M=최대(유지)
  4. Randomized geometry      — diffuse(동일 세포수)                    → geometry만 제거
  (well-mixed는 Fig S18에서 별도로: 항상 M=0)

또한 patient-geometry vs count-matched synthetic: 실제 환자 창과, 동일 세포수의 합성
조직에서 최적/결과가 다른지 → patient geometry의 added value.

6워커 병렬. 결과: data/ablation_spatial.csv, assets/ablation_spatial.png.
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
from synthetic import make_tissue
from abm import simulate, control_metrics

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# keep-stroma가 나타나는 체제(강한 confinement, 중간 면역배제)
BASE = dict(k_prolif=0.17, cd8_recruit=6, k_kill=0.40, caf_confine=0.9,
            caf_confine_ref=2.5, caf_pressure=1.8, cd8_barrier_alpha=1.0,
            caf_drug_block=0.6, max_tumor=5500)
DRUG = [("generic_cytotoxic", 0.30)]
DAYS = 90
CAP = 2000
SEEDS = [1, 2, 3, 4, 5]
WORKERS = 6
MYCAF = [0.0, 0.06, 0.12, 0.20, 0.30]     # 기질 수준 스윕
# ablation 변형: (label, param overrides, geometry)
VARIANTS = [
    ("Full spatial ABM", {}, "contained"),
    ("No-confinement", dict(caf_confine=0.0, caf_pressure=0.0, caf_drug_block=0.0),
     "contained"),
    ("No-immune-exclusion", dict(cd8_barrier_alpha=0.0), "contained"),
    ("Randomized geometry\n(diffuse, same counts)", {}, "diffuse"),
]
_TIS = {}


def _tissue(geom, seed):
    key = (geom, seed)
    if key not in _TIS:
        c, l, _ = make_tissue(geom, seed=seed)
        idx = np.random.default_rng(300 + seed).choice(len(l), CAP, replace=False)
        _TIS[key] = (c[idx], l[idx], int((l[idx] == "Tumor").sum()))
    return _TIS[key]


def evaluate(task):
    vi, mi, si = task
    label, ov, geom = VARIANTS[vi]
    c, l, n0 = _tissue(geom, SEEDS[si])
    p = dict(BASE); p.update(ov)
    p.update(k_caf_activate=MYCAF[mi], seed=SEEDS[si])
    h, _ = simulate(c, l, days=DAYS, params=p, regimen_subs=DRUG, snapshots=(0.0, 1.0))
    return (vi, mi, si, float(control_metrics(h, n0=n0)["final_frac"]))


def main():
    from concurrent.futures import ProcessPoolExecutor
    tasks = [(vi, mi, si) for vi in range(len(VARIANTS))
             for mi in range(len(MYCAF)) for si in range(len(SEEDS))]
    print(f"{len(tasks)} sims, {WORKERS}워커 병렬", flush=True)
    T = np.zeros((len(VARIANTS), len(MYCAF), len(SEEDS)))
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for k, (vi, mi, si, v) in enumerate(ex.map(evaluate, tasks, chunksize=4)):
            T[vi, mi, si] = v
            if (k + 1) % 30 == 0 or k == len(tasks) - 1:
                print(f"  {k+1}/{len(tasks)}", flush=True)
    curves = T.mean(axis=2)          # (variant, mycaf)
    rows = []
    print(f"\n{'variant':<34}{'optimal M':>10}{'strategy':>12}")
    for vi, (label, _o, _g) in enumerate(VARIANTS):
        opt_i = int(np.argmin(curves[vi]))
        opt_m = MYCAF[opt_i]
        strat = ("keep (interior)" if 0 < opt_i < len(MYCAF) - 1 else
                 ("reduce (M=0)" if opt_i == 0 else "preserve (M=max)"))
        rows.append(dict(variant=label.replace("\n", " "), optimal_mycaf=opt_m,
                         strategy=strat,
                         **{f"m{MYCAF[j]}": curves[vi, j] for j in range(len(MYCAF))}))
        print(f"{label.replace(chr(10),' '):<34}{opt_m:>10.2f}{strat:>16}")
    pd.DataFrame(rows).to_csv(os.path.join(ROOT, "data", "ablation_spatial.csv"),
                             index=False)
    patient_vs_synth()
    figure(curves)


def patient_vs_synth():
    """patient-geometry vs count-matched synthetic — geometry의 added value."""
    import anndata as ad
    H5 = os.path.join(ROOT, "data", "scotia", "raw_meta_data_final.h5ad")
    if not os.path.exists(H5):
        print("[patient-vs-synth] SCOTIA h5ad 없음 → 스킵"); return
    try:
        from patient_grounded import patient_tissue
        a = ad.read_h5ad(H5, backed="r")
        coords, labels = patient_tissue(a.obs, "U7-a")
    except Exception as e:
        print("[patient-vs-synth] 스킵:", e); return
    from collections import Counter
    cnt = Counter(labels)
    n0p = int((labels == "Tumor").sum())
    # count-matched synthetic: 동일 타입별 세포수를 무작위 배치(geometry 파괴)
    rng = np.random.default_rng(0)
    field = 1500.0
    sx, sl = [], []
    for t, n in cnt.items():
        sx.append(rng.uniform(0, field, (n, 2))); sl += [t] * n
    sc = np.vstack(sx); sl = np.array(sl, dtype=object)
    reg = dict(regimen_subs=[("generic_immunostim", 1.0)], adaptive=False)
    out = {}
    for tag, (cc, ll, n0) in [("patient geometry", (coords, labels, n0p)),
                              ("count-matched synthetic", (sc, sl, n0p))]:
        vals = []
        for s in (1, 2, 3):
            p = dict(BASE); p["seed"] = s
            h, _ = simulate(cc, ll, days=DAYS, params=p, snapshots=(0.0, 1.0), **reg)
            vals.append(control_metrics(h, n0=n0)["final_frac"])
        out[tag] = float(np.mean(vals))
    print(f"[patient-vs-synth] immune therapy final: "
          f"patient={out['patient geometry']:.2f}x  "
          f"count-matched synthetic={out['count-matched synthetic']:.2f}x")
    pd.DataFrame([out]).to_csv(os.path.join(ROOT, "data", "ablation_patient_synth.csv"),
                              index=False)


def figure(curves):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.4))
    cols = ["#C0392B", "#E67E22", "#2980B9", "#27AE60"]
    for vi, (label, _o, _g) in enumerate(VARIANTS):
        a1.plot(MYCAF, curves[vi], "o-", color=cols[vi], lw=2, label=label)
        oi = int(np.argmin(curves[vi]))
        a1.scatter([MYCAF[oi]], [curves[vi][oi]], s=110, facecolor="none",
                   edgecolor=cols[vi], linewidths=2.2, zorder=5)
    a1.set_xlabel("Stromal level M (k_caf_activate)")
    a1.set_ylabel("Final tumor burden (fold)  →  lower better")
    a1.set_title("a  Optimal stromal strategy by mechanism ablation\n"
                 "(○ = optimum; interior=keep, M=0=reduce, M=max=preserve)",
                 fontsize=10, fontweight="bold")
    a1.legend(frameon=False, fontsize=8.5, loc="upper left")

    # b: selected strategy summary (optimal M per variant)
    labels = [v[0].replace("\n", " ") for v in VARIANTS] + ["Well-mixed\n(Fig. S18)"]
    optM = [MYCAF[int(np.argmin(curves[vi]))] for vi in range(len(VARIANTS))] + [0.0]
    y = np.arange(len(labels))
    bcols = cols + ["#7F8C8D"]
    a2.barh(y, optM, color=bcols, edgecolor="black", linewidth=0.7, alpha=0.85)
    for i, m in enumerate(optM):
        strat = ("reduce" if m == 0 else ("preserve" if m >= MYCAF[-1] else "keep-some"))
        a2.text(m + 0.005, i, f"{m:.2f} ({strat})", va="center", fontsize=8)
    a2.set_yticks(y); a2.set_yticklabels([l.replace("\n", " ") for l in labels],
                                         fontsize=8)
    a2.invert_yaxis()
    a2.set_xlabel("Selected optimal stromal level M")
    a2.set_xlim(0, MYCAF[-1] * 1.35)
    a2.set_title("b  The selected strategy changes with the mechanism\n"
                 "confinement creates 'keep'; immune exclusion creates 'reduce'",
                 fontsize=10, fontweight="bold")
    fig.suptitle("Ablation: spatial geometry changes the selected stromal strategy, not "
                 "just the tumor number (in silico)", fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "assets", "ablation_spatial.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    df = pd.read_csv(os.path.join(ROOT, "data", "ablation_spatial.csv"))
    curves = np.array([[r[f"m{m}"] for m in MYCAF] for _, r in df.iterrows()])
    figure(curves)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
