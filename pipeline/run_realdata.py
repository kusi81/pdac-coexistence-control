"""
실데이터 장벽 점수 — squidpy MIBI-TOF(대장암)에 barrier score 적용.

합성 조직이 아닌 실제 공간 단일세포 데이터. 시료(point)별로 분리해:
  종양 = Epithelial, 면역 = Tcell_CD8, 장벽 = Fibroblast
로 매핑하고, Fibroblast가 종양↔CD8 경로에 개재되는지(containment)를 검정한다.

주의: 좌표는 픽셀 단위(FOV 1024px ≈ 800µm, ~0.8µm/px). corridor는 픽셀 기준.
합성 검증과 달리 정답(ground truth)이 없으므로, 관측 vs matched-count 랜덤기질
귀무 비교로 '위치 우연 이상' 여부만 판정한다.
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
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

import data_loader as dl
from spatial_core import barrier_score, radial_profile, median_nn_distance

CORRIDOR = 25          # 픽셀 (~20µm)
COLORS = {
    "Epithelial": "#3D3D3D", "Fibroblast": "#C0392B", "Tcell_CD8": "#27AE60",
    "Tcell_CD4": "#2ECC71", "Myeloid_CD68": "#E67E22", "Myeloid_CD11c": "#F39C12",
    "Endothelial": "#2E86AB", "Imm_other": "#95A5A6",
}


def main():
    print("MIBI-TOF 로드 중...")
    adata = dl.load_mibitof()
    ck = dl.find_celltype_key(adata)
    roles = dl.guess_roles(np.asarray(adata.obs[ck]).astype(str))
    tumor, immune, barrier = roles["tumor"], roles["immune"], roles["barrier"]
    print(f"매핑 → 종양={tumor}, 면역={immune}, 장벽={barrier}")

    libs = dl.libraries(adata)
    results = {}
    for lib in libs:
        coords, labels = dl.to_coords_labels(adata, celltype_key=ck, library=lib)
        bs = barrier_score(coords, labels, tumor=tumor, barrier=barrier,
                           immune=immune, corridor_um=CORRIDOR, seed=0)
        med, n = median_nn_distance(coords, labels, tumor, immune)
        results[lib] = dict(bs=bs, coords=coords, labels=labels,
                            cd8_med=med, n=len(labels))
        if bs:
            z = bs["z_score"]
            print(f"[{lib}] n={len(labels)}  {barrier} 장벽 z="
                  f"{'N/A' if z is None else round(z,2)}  "
                  f"obs={bs['observed_corridor_density']:.2f} null={bs['null_mean']:.2f} "
                  f"p={bs['p_value']:.3g} → {bs['interpretation']}")
        else:
            print(f"[{lib}] n={len(labels)}  세포 수 부족으로 판정 불가")

    _figure(results, tumor, immune, barrier)
    return results


def _figure(results, tumor, immune, barrier):
    libs = list(results)
    fig, axes = plt.subplots(2, len(libs), figsize=(5.2 * len(libs), 9),
                             gridspec_kw={"height_ratios": [3, 1.3]})
    if len(libs) == 1:
        axes = axes.reshape(2, 1)
    order = ["Endothelial", "Imm_other", "Myeloid_CD11c", "Myeloid_CD68",
             "Tcell_CD4", "Fibroblast", "Tcell_CD8", "Epithelial"]
    for j, lib in enumerate(libs):
        r = results[lib]
        coords, labels = r["coords"], r["labels"]
        ax = axes[0, j]
        present = [c for c in order if c in set(labels)] + \
                  [c for c in set(labels) if c not in order]
        for ct in present:
            m = labels == ct
            ax.scatter(coords[m, 0], coords[m, 1], s=9,
                       c=COLORS.get(ct, "#999"), alpha=0.8, linewidths=0,
                       label=f"{ct} ({int(m.sum())})" if j == len(libs) - 1 else None)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_title(f"{lib}  (n={r['n']})", fontsize=11, fontweight="bold")
        ax.set_xlabel("x (px)")
        if j == 0:
            ax.set_ylabel("y (px)")
        if j == len(libs) - 1:
            ax.legend(frameon=False, fontsize=7, markerscale=1.8,
                      loc="center left", bbox_to_anchor=(1.01, 0.5))

        axb = axes[1, j]
        bs = r["bs"]
        if bs and bs["z_score"] is not None:
            z = bs["z_score"]
            axb.bar([0], [bs["observed_corridor_density"]], 0.4, color=COLORS[barrier],
                    label="observed")
            axb.bar([0.5], [bs["null_mean"]], 0.4, yerr=bs["null_sd"], color="#BBB",
                    capsize=3, label="null")
            axb.set_xticks([0.25]); axb.set_xticklabels([f"z={z:+.1f}"])
            axb.set_title("interposed ✓" if z > 2 else "no evidence", fontsize=9.5,
                          color=("#C0392B" if z > 2 else "#7F8C8D"))
            if j == 0:
                axb.set_ylabel("corridor density")
            if j == len(libs) - 1:
                axb.legend(frameon=False, fontsize=7)
        else:
            axb.axis("off")

    fig.suptitle(f"실데이터 (MIBI-TOF 대장암) — {barrier} 장벽 점수: "
                 f"{barrier}가 {tumor}↔{immune} 사이에 개재되는가", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "realdata_mibitof.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
