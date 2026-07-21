"""대시보드 rim 패널을 MIBI-TOF(대장암)로 검증 — run_rim_panel + fig_rim 그대로 사용.
종양=Epithelial. 각 시료(point)에서 '누가 Epithelial에 인접한가'를 계산·시각화."""
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

import data_loader as dl
from analysis import run_rim_panel      # 대시보드와 동일 계산
from plots import fig_rim, COLORS       # 대시보드와 동일 그림

SHELL = 30   # px (대시보드 회랑 폭 기본값과 동일)


def main():
    adata = dl.load_mibitof()
    ck = dl.find_celltype_key(adata)
    libs = dl.libraries(adata)
    fig, axes = plt.subplots(1, len(libs), figsize=(6.0 * len(libs), 4.6))
    axes = np.atleast_1d(axes)
    for ax, lib in zip(axes, libs):
        coords, labels = dl.to_coords_labels(adata, celltype_key=ck, library=lib)
        rows = run_rim_panel(coords, labels, tumor="Epithelial",
                             shell_um=SHELL, n_perms=500)
        rows = sorted(rows, key=lambda d: d["z"])
        names = [r["cell_type"] for r in rows]
        zs = [r["z"] for r in rows]
        cols = [COLORS.get(n, "#999") for n in names]
        ax.barh(range(len(names)), zs, color=cols, height=0.62)
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
        ax.axvline(0, color="black", lw=1)
        ax.axvline(2, color="#888", ls="--", lw=1); ax.axvline(-2, color="#888", ls="--", lw=1)
        for i, z in enumerate(zs):
            ax.text(z + (0.5 if z >= 0 else -0.5), i, f"{z:+.0f}", va="center",
                    ha="left" if z >= 0 else "right", fontsize=8.5, fontweight="bold")
        top = max(rows, key=lambda d: d["z"])
        ax.set_title(f"{lib}  (최인접: {top['cell_type']} z={top['z']:+.0f})", fontsize=10)
        ax.set_xlabel("Epithelial rim(30px) 농축 z")
        print(f"[{lib}] " + ", ".join(f"{r['cell_type']} {r['z']:+.0f}" for r in
                                      sorted(rows, key=lambda d: -d['z'])))
    fig.suptitle("대시보드 rim 패널 검증 — MIBI-TOF 대장암: 누가 Epithelial(종양)에 인접한가",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "rim_panel_mibitof.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
