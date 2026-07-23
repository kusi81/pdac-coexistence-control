"""Figure 4 재계산 — 더 촘촘한 grid + 다중 seed + phase uncertainty (리뷰어 #1).

기존 phase map은 4×4 grid·cell당 3 seed로 해상도·불확실성이 부족했다. 여기서는
6×6 grid, cell당 5 seed로 재계산하고, 각 (confinement, immune-exclusion) cell에서
'기질 유지가 유리(자원 regime)'일 확률 P(keep)을 seed 간 변동으로 정량화한다.
단순 keep/reduce 이분이 아니라 확률·이득의 신뢰구간으로 phase 경계 불확실성을 표시.

6워커 병렬, max_tumor 상한. 결과: data/phase_map_dense.npz, assets/phase_dense.png.
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
BASE = dict(k_prolif=0.17, cd8_recruit=6, k_kill=0.40,
            caf_confine=0.9, caf_confine_ref=2.5, max_tumor=5500)
DRUG = [("generic_cytotoxic", 0.30)]
DAYS = 90
CAP = 2000
WORKERS = 6
PRESSURE = [0.0, 0.6, 1.2, 1.8, 2.4, 3.0]     # x: confinement (6)
ALPHA = [0.2, 0.6, 1.0, 1.4, 1.8, 2.4]        # y: immune exclusion (6)
MYCAF = [0.0, 0.06, 0.12, 0.20, 0.30]         # myCAF level sweep (5)
SEEDS = [1, 2, 3, 4, 5]                        # 5 seeds/cell
_TIS = {}


def _tissue(seed):
    if seed not in _TIS:
        c, l, _ = make_tissue("contained", seed=42)
        idx = np.random.default_rng(seed).choice(len(l), CAP, replace=False)
        _TIS[seed] = (c[idx], l[idx], int((l[idx] == "Tumor").sum()))
    return _TIS[seed]


def evaluate(task):
    xi, yi, mi, si = task
    c, l, n0 = _tissue(SEEDS[si])
    p = dict(BASE)
    p.update(caf_pressure=PRESSURE[xi], cd8_barrier_alpha=ALPHA[yi],
             k_caf_activate=MYCAF[mi], seed=SEEDS[si])
    h, _ = simulate(c, l, days=DAYS, params=p, regimen_subs=DRUG, snapshots=(0.0, 1.0))
    return (xi, yi, mi, si, float(control_metrics(h, n0=n0)["final_frac"]))


def main():
    from concurrent.futures import ProcessPoolExecutor
    tasks = [(xi, yi, mi, si) for xi in range(len(PRESSURE))
             for yi in range(len(ALPHA)) for mi in range(len(MYCAF))
             for si in range(len(SEEDS))]
    print(f"{len(tasks)} sims ({len(PRESSURE)}×{len(ALPHA)}×{len(MYCAF)}×{len(SEEDS)}), "
          f"{WORKERS}워커 병렬", flush=True)
    T = np.zeros((len(PRESSURE), len(ALPHA), len(MYCAF), len(SEEDS)))
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for k, (xi, yi, mi, si, v) in enumerate(ex.map(evaluate, tasks, chunksize=4)):
            T[xi, yi, mi, si] = v
            if (k + 1) % 90 == 0 or k == len(tasks) - 1:
                print(f"  {k+1}/{len(tasks)}", flush=True)
    np.savez(os.path.join(ROOT, "data", "phase_map_dense.npz"),
             T=T, pressure=PRESSURE, alpha=ALPHA, mycaf=MYCAF, seeds=SEEDS)
    # per-cell, per-seed: optimal mycaf>0? and benefit
    keep_prob = np.zeros((len(PRESSURE), len(ALPHA)))
    ben_mean = np.zeros_like(keep_prob); ben_lo = np.zeros_like(keep_prob)
    ben_hi = np.zeros_like(keep_prob)
    for xi in range(len(PRESSURE)):
        for yi in range(len(ALPHA)):
            keeps, bens = [], []
            for si in range(len(SEEDS)):
                curve = T[xi, yi, :, si]
                opt = int(np.argmin(curve))
                keeps.append(1.0 if opt > 0 and curve[0] - curve[opt] > 0.02 else 0.0)
                bens.append(curve[0] - curve.min())
            keep_prob[xi, yi] = np.mean(keeps)
            ben_mean[xi, yi] = np.mean(bens)
            ben_lo[xi, yi], ben_hi[xi, yi] = np.percentile(bens, [10, 90])
    print(f"\n자원 regime(P(keep)≥0.6) 셀수 = {(keep_prob>=0.6).sum()}/{keep_prob.size}")
    figure(keep_prob, ben_mean, ben_lo, ben_hi)


def figure(keep_prob, ben_mean, ben_lo, ben_hi):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.6))
    ext = [-.5, len(PRESSURE) - .5, -.5, len(ALPHA) - .5]
    # a: P(keep) heatmap — phase uncertainty
    im1 = a1.imshow(keep_prob.T, origin="lower", cmap="RdYlGn", vmin=0, vmax=1,
                    aspect="auto")
    for xi in range(len(PRESSURE)):
        for yi in range(len(ALPHA)):
            a1.text(xi, yi, f"{keep_prob[xi,yi]:.1f}", ha="center", va="center",
                    fontsize=7.5, color="black")
    a1.set_xticks(range(len(PRESSURE))); a1.set_xticklabels(PRESSURE)
    a1.set_yticks(range(len(ALPHA))); a1.set_yticklabels(ALPHA)
    a1.set_xlabel("Physical confinement strength (caf_pressure)")
    a1.set_ylabel("Immune-exclusion strength (cd8_barrier_alpha)")
    a1.set_title("a  P(keep-stroma is optimal) across 5 seeds\n"
                 "(phase-boundary uncertainty, not a hard line)", fontsize=10,
                 fontweight="bold")
    fig.colorbar(im1, ax=a1, fraction=0.045, pad=0.02, label="P(keep stroma)")
    # b: mean benefit with 10–90% band width as hatching alpha
    im2 = a2.imshow(ben_mean.T, origin="lower", cmap="viridis", aspect="auto",
                    vmin=0)
    for xi in range(len(PRESSURE)):
        for yi in range(len(ALPHA)):
            ci = ben_hi[xi, yi] - ben_lo[xi, yi]
            a2.text(xi, yi, f"±{ci/2:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white")
    a2.set_xticks(range(len(PRESSURE))); a2.set_xticklabels(PRESSURE)
    a2.set_yticks(range(len(ALPHA))); a2.set_yticklabels(ALPHA)
    a2.set_xlabel("Physical confinement strength (caf_pressure)")
    a2.set_ylabel("Immune-exclusion strength (cd8_barrier_alpha)")
    a2.set_title("b  Mean benefit of keeping stroma\n(±half the 10–90% seed interval)",
                 fontsize=10, fontweight="bold")
    fig.colorbar(im2, ax=a2, fraction=0.045, pad=0.02,
                 label="tumor(no stroma) − tumor(optimal)")
    fig.suptitle("Figure 4 (dense) — condition-dependent stromal optimum on a 6×6 grid "
                 "with 5-seed phase uncertainty (in silico)", fontsize=11.5,
                 fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(ROOT, "assets", "phase_dense.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    z = np.load(os.path.join(ROOT, "data", "phase_map_dense.npz"))
    T = z["T"]
    kp = np.zeros((T.shape[0], T.shape[1])); bm = np.zeros_like(kp)
    bl = np.zeros_like(kp); bh = np.zeros_like(kp)
    for xi in range(T.shape[0]):
        for yi in range(T.shape[1]):
            keeps, bens = [], []
            for si in range(T.shape[3]):
                curve = T[xi, yi, :, si]; opt = int(np.argmin(curve))
                keeps.append(1.0 if opt > 0 and curve[0] - curve[opt] > 0.02 else 0.0)
                bens.append(curve[0] - curve.min())
            kp[xi, yi] = np.mean(keeps); bm[xi, yi] = np.mean(bens)
            bl[xi, yi], bh[xi, yi] = np.percentile(bens, [10, 90])
    figure(kp, bm, bl, bh)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
