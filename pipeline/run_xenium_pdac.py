"""
췌장암 Xenium 단일세포 실데이터 장벽 점수 — GSE274673 (10x Xenium PDAC).

Zhou Visium(spot·약함)와 달리 **단일세포 해상도 + µm 좌표**라 corridor_um이 실제
마이크론으로 해석되고, myCAF/iCAF 검정이 원리적으로 유효하다. 세포타입 주석이 없어
마커 기반 단일세포 타이핑 후:
  종양=Tumor, 면역=CD8_T, 장벽=myCAF, 대조 장벽=iCAF
로 containment를 검정한다. 양성대조: myCAF가 iCAF보다 종양에 가까운가(문헌).
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

import data_loader as dl
from spatial_core import barrier_score, proximity_test, median_nn_distance

COLORS = {"Tumor": "#3D3D3D", "myCAF": "#C0392B", "iCAF": "#2E86AB",
          "CD8_T": "#27AE60", "Macrophage": "#E67E22", "apCAF": "#8E44AD",
          "Endothelial": "#16A085"}
CORRIDOR = 30.0   # µm (Xenium 좌표가 µm이므로 실제 30µm 회랑)


def analyze(bundle):
    adata = dl.load_xenium_bundle(bundle)
    labels, cov = dl.annotate_pdac(adata)   # 저자 2단계 module-score 방법
    coords = np.asarray(adata.obsm["spatial"], float)
    res = {"coords": coords, "labels": labels, "n": len(labels), "cov": cov,
           "tag": os.path.basename(bundle)[:32]}
    res["bs_my"] = barrier_score(coords, labels, tumor="Tumor", barrier="myCAF",
                                 immune="CD8_T", corridor_um=CORRIDOR, seed=0)
    res["bs_ic"] = barrier_score(coords, labels, tumor="Tumor", barrier="iCAF",
                                 immune="CD8_T", corridor_um=CORRIDOR, seed=0)
    res["pt"] = proximity_test(coords, labels, "Tumor", "myCAF", "iCAF", n_perms=500)
    return res


def main():
    bundles = sorted(glob.glob(os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "data", "xenium", "output-*")))
    bundles = [b for b in bundles if os.path.isdir(b)]
    print(f"Xenium 시료 {len(bundles)}개 발견")
    results = []
    import pandas as pd
    for b in bundles:
        r = analyze(b)
        results.append(r)
        print(f"\n[{r['tag']}] n={r['n']} 세포, 마커 {r['cov']:.0%}")
        print("  세포타입:", dict(pd.Series(r["labels"]).value_counts()))
        for nm, bs in [("myCAF", r["bs_my"]), ("iCAF", r["bs_ic"])]:
            if bs and bs["z_score"] is not None:
                print(f"  {nm} 장벽 z={bs['z_score']:+.2f}  obs={bs['observed_corridor_density']:.2f} "
                      f"null={bs['null_mean']:.2f} p={bs['p_value']:.3g} → {bs['interpretation']}")
        pt = r["pt"]
        if pt:
            ok = "✓문헌일치" if pt["interpretation"].startswith("myCAF") else "✗문헌불일치"
            print(f"  근접성: myCAF {pt['median_dist_a']:.0f}µm vs iCAF {pt['median_dist_b']:.0f}µm "
                  f"(p={pt['p_permutation']:.3g}) → {pt['interpretation']} [{ok}]")
    _figure(results)
    return results


def _figure(results):
    nc = len(results)
    fig, axes = plt.subplots(2, nc, figsize=(6.0 * nc, 10),
                             gridspec_kw={"height_ratios": [3, 1.4]})
    axes = np.atleast_2d(axes)
    if nc == 1:
        axes = axes.reshape(2, 1)
    for j, r in enumerate(results):
        coords, labels = r["coords"], r["labels"]
        ax = axes[0, j]
        for ct in ["Macrophage", "iCAF", "myCAF", "CD8_T", "Tumor"]:
            m = labels == ct
            if m.any():
                ax.scatter(coords[m, 0], coords[m, 1], s=2, c=COLORS[ct],
                           alpha=0.6, linewidths=0,
                           label=f"{ct} ({int(m.sum())})" if j == nc - 1 else None)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_title(f"{r['tag']}\nn={r['n']} 세포 (단일세포)", fontsize=10)
        ax.set_xlabel("x (µm)")
        if j == 0:
            ax.set_ylabel("y (µm)")
        if j == nc - 1:
            ax.legend(frameon=False, fontsize=7.5, markerscale=3,
                      loc="center left", bbox_to_anchor=(1.01, 0.5))

        axb = axes[1, j]
        names, obs, nul, sd, zs = [], [], [], [], []
        for nm, bs in [("myCAF", r["bs_my"]), ("iCAF", r["bs_ic"])]:
            if bs and bs["z_score"] is not None:
                names.append(nm); obs.append(bs["observed_corridor_density"])
                nul.append(bs["null_mean"]); sd.append(bs["null_sd"])
                zs.append(bs["z_score"])
        x = np.arange(len(names))
        axb.bar(x - 0.19, obs, 0.38, color=[COLORS[n] for n in names], label="observed")
        axb.bar(x + 0.19, nul, 0.38, yerr=sd, color="#BBB", capsize=3, label="null")
        axb.set_xticks(x)
        axb.set_xticklabels([f"{n}\nz={z:+.1f}" for n, z in zip(names, zs)])
        axb.set_ylabel("corridor density")
        pt = r["pt"]
        if pt:
            axb.set_title(f"myCAF {pt['median_dist_a']:.0f} vs iCAF "
                          f"{pt['median_dist_b']:.0f}µm→종양", fontsize=9)
        if j == nc - 1:
            axb.legend(frameon=False, fontsize=7.5)

    fig.suptitle("췌장암 Xenium 단일세포 실데이터 (GSE274673) — myCAF containment 검정 "
                 "(µm 좌표, 마커 타이핑)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "realdata_pdac_xenium.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("\nwrote", out)


if __name__ == "__main__":
    main()
