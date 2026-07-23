"""CAF pro-tumor 축 복원 — phase map 강건성 검정 (리뷰어 지적 (ii)).

리뷰어 지적: PDAC에서 caf_protumor=0으로 두면 CAF의 종양지원 생물학(파라크린 증식,
IL-6/JAK-STAT 치료유발 생존신호=약물내성)을 빠뜨려 결과가 'containment(기질=자원)'
쪽으로 구조적으로 편향된다.

이 스크립트는 동일 조직·격자에서 phase map을 두 조건으로 계산:
  OFF : caf_protumor=0,   caf_survival=0    (원래 가정)
  ON  : caf_protumor=0.25, caf_survival=0.4 (복원한 PDAC pro-tumor/내성 생물학)
각 (confinement=caf_pressure, immune-exclusion=cd8_barrier_alpha)에서 myCAF 수준을
스윕해 '기질 유지의 이득' = tumor(무기질) − tumor(최적기질) 를 구한다.
이득>0 = 자원 regime. ON에서 이 regime이 얼마나 줄어드는지(그러나 존속하는지)를 본다.

6워커 병렬 + max_tumor 상한. 결과: data/phase_map_protumor.npz, assets/phase_protumor.png.
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
SEEDS = [1, 2, 3]
CAP = 2500
WORKERS = 6
PRESSURE = [0.0, 0.9, 1.8, 2.8]          # x: confinement 강도
ALPHA = [0.2, 0.6, 1.2, 2.0]             # y: 면역배제 강도
MYCAF = [0.0, 0.08, 0.16, 0.26]          # myCAF 수준 스윕
# 두 조건: 복원 전/후
CONDS = {"off": dict(caf_protumor=0.0, caf_survival=0.0),
         "on":  dict(caf_protumor=0.25, caf_survival=0.4)}
CONDKEYS = ["off", "on"]

_T = {}


def _tissue():
    if "c" not in _T:
        c, l, _ = make_tissue("contained", seed=42)
        idx = np.random.default_rng(7).choice(len(l), CAP, replace=False)
        _T["c"], _T["l"], _T["n0"] = c[idx], l[idx], int((l[idx] == "Tumor").sum())
    return _T["c"], _T["l"], _T["n0"]


def evaluate(task):
    ci, xi, yi, mi, si = task
    c, l, n0 = _tissue()
    p = dict(BASE); p.update(CONDS[CONDKEYS[ci]])
    p.update(caf_pressure=PRESSURE[xi], cd8_barrier_alpha=ALPHA[yi],
             k_caf_activate=MYCAF[mi], seed=SEEDS[si])
    h, _ = simulate(c, l, days=DAYS, params=p, regimen_subs=DRUG, snapshots=(0.0, 1.0))
    return (ci, xi, yi, mi, si, float(control_metrics(h, n0=n0)["final_frac"]))


def main():
    from concurrent.futures import ProcessPoolExecutor
    tasks = [(ci, xi, yi, mi, si)
             for ci in range(2) for xi in range(len(PRESSURE))
             for yi in range(len(ALPHA)) for mi in range(len(MYCAF))
             for si in range(len(SEEDS))]
    print(f"{len(tasks)} sims (2 conds × {len(PRESSURE)}×{len(ALPHA)}×{len(MYCAF)} "
          f"× {len(SEEDS)} seeds), {WORKERS}워커 병렬", flush=True)
    T = np.zeros((2, len(PRESSURE), len(ALPHA), len(MYCAF), len(SEEDS)))
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for k, (ci, xi, yi, mi, si, v) in enumerate(
                ex.map(evaluate, tasks, chunksize=4)):
            T[ci, xi, yi, mi, si] = v
            if (k + 1) % 60 == 0 or k == len(tasks) - 1:
                print(f"  {k+1}/{len(tasks)}", flush=True)
    Tm = T.mean(axis=4)                       # seed 평균 → (2, x, y, mycaf)
    np.savez(os.path.join(ROOT, "data", "phase_map_protumor.npz"),
             T=Tm, pressure=PRESSURE, alpha=ALPHA, mycaf=MYCAF, conds=CONDKEYS)

    # 자원 regime 판정: 이득 = tumor(mycaf=0) − min_over_mycaf(tumor)
    for ci, ck in enumerate(CONDKEYS):
        opt = np.argmin(Tm[ci], axis=2)                       # (x,y) 최적 mycaf index
        benefit = Tm[ci, :, :, 0] - Tm[ci].min(axis=2)        # ≥0
        resource = (opt > 0) & (benefit > 0.02)
        print(f"[{ck}] 자원 regime 셀수 = {int(resource.sum())}/{resource.size}  "
              f"평균이득={benefit[resource].mean() if resource.any() else 0:.3f}")
    figure(Tm)


def figure(Tm):
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
    ben = {}
    vmax = 0
    for ci in range(2):
        ben[ci] = Tm[ci, :, :, 0] - Tm[ci].min(axis=2)   # benefit of keeping stroma
        vmax = max(vmax, ben[ci].max())
    for ci, ck in enumerate(("off", "on")):
        ax = axs[ci]
        B = ben[ci].T                                    # rows=alpha(y), cols=pressure(x)
        im = ax.imshow(B, origin="lower", cmap="RdYlGn", vmin=0, vmax=max(vmax, 0.05),
                       aspect="auto")
        opt = np.argmin(Tm[ci], axis=2).T
        for yi in range(len(ALPHA)):
            for xi in range(len(PRESSURE)):
                res = opt[yi, xi] > 0 and B[yi, xi] > 0.02
                ax.text(xi, yi, "keep" if res else "reduce",
                        ha="center", va="center", fontsize=8.5,
                        color="#0B3D1A" if res else "#4A4A4A",
                        fontweight="bold" if res else "normal")
        ax.set_xticks(range(len(PRESSURE))); ax.set_xticklabels(PRESSURE)
        ax.set_yticks(range(len(ALPHA))); ax.set_yticklabels(ALPHA)
        ax.set_xlabel("Physical confinement strength (caf_pressure)")
        if ci == 0:
            ax.set_ylabel("Immune-exclusion strength (cd8_barrier_alpha)")
        nres = int(((np.argmin(Tm[ci], axis=2) > 0) &
                    (ben[ci] > 0.02)).sum())
        title = ("a  CAF pro-tumor biology OFF\n(original assumption)"
                 if ci == 0 else
                 "b  CAF pro-tumor biology ON\n(paracrine + IL-6 drug tolerance)")
        ax.set_title(f"{title}   —   resource regime: "
                     f"{nres}/{ben[ci].size} cells", fontsize=9.8, fontweight="bold")
    fig.subplots_adjust(right=0.88)
    cax = fig.add_axes([0.90, 0.15, 0.018, 0.62])
    fig.colorbar(im, cax=cax,
                 label="Benefit of keeping stroma\n= tumor(no stroma) − tumor(optimal)")
    fig.suptitle("Restoring CAF pro-tumor/drug-tolerance biology shrinks—but does not "
                 "erase—the 'stroma-as-resource' regime\n(green/keep = preserving stroma "
                 "aids control; red/reduce = deplete; in silico)", fontsize=11.5,
                 fontweight="bold")
    fig.subplots_adjust(left=0.07, right=0.88, top=0.82, bottom=0.12, wspace=0.12)
    out = os.path.join(ROOT, "assets", "phase_protumor.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    z = np.load(os.path.join(ROOT, "data", "phase_map_protumor.npz"))
    figure(z["T"])


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
