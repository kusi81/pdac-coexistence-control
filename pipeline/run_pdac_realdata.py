"""
췌장암(PDAC) 실데이터 장벽 점수 — Figshare Zhou et al. Visium 데이터셋.

MIBI-TOF(대장암)와 달리 실제 PDAC 조직. 단, Visium(spot 단위)이라:
  - 세포타입 주석이 없으면 마커 기반 argmax로 타이핑(spot당 우세 타입 근사)
  - 해상도 낮아 containment 주장은 '약함'(README 플랫폼 상한)
  - 좌표는 픽셀 → corridor를 spot 간격에 맞춰 자동 스케일

종양=Tumor, 면역=CD8_T, 장벽=myCAF, 대조 장벽=iCAF 로 검정.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

import data_loader as dl
from spatial_core import barrier_score, proximity_test, median_nn_distance

COLORS = {"Tumor": "#3D3D3D", "myCAF": "#C0392B", "iCAF": "#2E86AB",
          "CD8_T": "#27AE60", "Macrophage": "#E67E22"}


def _auto_corridor(coords):
    """spot 간 중앙 최근접거리의 ~0.75배를 corridor로."""
    d, _ = cKDTree(coords).query(coords, k=2)
    med = float(np.median(d[:, 1]))
    return max(10.0, 0.75 * med), med


def main(index=0):
    print("Zhou PDAC 데이터 로드/압축해제 중...")
    adata, files = dl.load_zhou_pdac(index=index)
    print(f"h5ad {len(files)}개 중 [{index}] = {os.path.basename(files[index])}")
    print(f"shape: {adata.shape}")

    # 세포타입 주석 확인 → 없으면 마커 기반 타이핑
    ck = None
    try:
        ck = dl.find_celltype_key(adata)
    except Exception:
        pass
    has_needed = ck and set(["Tumor", "myCAF", "CD8_T"]).issubset(
        set(np.asarray(adata.obs[ck]).astype(str))) if ck else False
    if not has_needed:
        labels, cov = dl.annotate_by_markers(adata)
        print(f"마커 기반 타이핑 적용 (패널 유전자 커버리지 {cov:.0%})")
    else:
        labels = np.asarray(adata.obs[ck]).astype(str)
        print(f"기존 주석 사용: {ck}")

    import pandas as pd
    print("세포타입 분포:", dict(pd.Series(labels).value_counts()))
    coords = np.asarray(adata.obsm["spatial"], float)[:, :2]
    corridor, med = _auto_corridor(coords)
    print(f"spot 중앙간격 {med:.0f}px → corridor {corridor:.0f}px")

    bs_my = barrier_score(coords, labels, tumor="Tumor", barrier="myCAF",
                          immune="CD8_T", corridor_um=corridor, seed=0)
    bs_ic = barrier_score(coords, labels, tumor="Tumor", barrier="iCAF",
                          immune="CD8_T", corridor_um=corridor, seed=0)
    pt = proximity_test(coords, labels, "Tumor", "myCAF", "iCAF", n_perms=500)

    for name, bs in [("myCAF", bs_my), ("iCAF", bs_ic)]:
        if bs and bs["z_score"] is not None:
            print(f"  {name} 장벽 z={bs['z_score']:+.2f}  "
                  f"obs={bs['observed_corridor_density']:.2f} null={bs['null_mean']:.2f} "
                  f"p={bs['p_value']:.3g} → {bs['interpretation']}")
        else:
            print(f"  {name}: 세포 수 부족")
    if pt:
        print(f"  근접성: myCAF {pt['median_dist_a']:.0f}px vs iCAF "
              f"{pt['median_dist_b']:.0f}px (p={pt['p_permutation']:.3g}) "
              f"→ {pt['interpretation']}")

    _figure(coords, labels, bs_my, bs_ic, os.path.basename(files[index]))
    return dict(bs_my=bs_my, bs_ic=bs_ic, pt=pt)


def _figure(coords, labels, bs_my, bs_ic, tag):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.6),
                                   gridspec_kw={"width_ratios": [3, 2]})
    for ct in ["Macrophage", "iCAF", "myCAF", "CD8_T", "Tumor"]:
        m = labels == ct
        if m.any():
            ax1.scatter(coords[m, 0], coords[m, 1], s=14, c=COLORS.get(ct, "#999"),
                        alpha=0.8, linewidths=0, label=f"{ct} ({int(m.sum())})")
    ax1.set_aspect("equal"); ax1.invert_yaxis()
    ax1.set_title(f"PDAC Visium (마커 타이핑)\n{tag}", fontsize=10)
    ax1.set_xlabel("x (px)"); ax1.set_ylabel("y (px)")
    ax1.legend(frameon=False, fontsize=7.5, markerscale=1.6,
               loc="center left", bbox_to_anchor=(1.01, 0.5))

    names, zs, obs, nul, sd = [], [], [], [], []
    for nm, bs in [("myCAF", bs_my), ("iCAF", bs_ic)]:
        if bs and bs["z_score"] is not None:
            names.append(nm); zs.append(bs["z_score"])
            obs.append(bs["observed_corridor_density"]); nul.append(bs["null_mean"])
            sd.append(bs["null_sd"])
    x = np.arange(len(names))
    ax2.bar(x - 0.19, obs, 0.38, color=[COLORS[n] for n in names], label="observed")
    ax2.bar(x + 0.19, nul, 0.38, yerr=sd, color="#BBB", capsize=3, label="null")
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{n}\nz={z:+.1f}" for n, z in zip(names, zs)])
    ax2.set_ylabel("corridor density")
    ax2.set_title("장벽 점수 (observed vs matched-null)", fontsize=10)
    ax2.legend(frameon=False, fontsize=8)

    fig.suptitle("췌장암 실데이터 (Zhou Visium) — myCAF containment 검정 "
                 "(⚠️ Visium=약함, 마커 argmax 근사)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "realdata_pdac_zhou.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main(index=int(sys.argv[1]) if len(sys.argv) > 1 else 0)
