"""개선 지표 시각화: 종양 주변(rim) 세포타입 농축 — 실제 PDAC에서 뭐가 juxtatumoral인가."""
import sys, os, glob, warnings
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
from spatial_core import rim_enrichment

COLORS = {"Tumor": "#3D3D3D", "myCAF": "#C0392B", "iCAF": "#2E86AB",
          "CD8_T": "#27AE60", "Macrophage": "#E67E22", "apCAF": "#8E44AD",
          "Endothelial": "#16A085"}
CTS = ["myCAF", "iCAF", "apCAF", "Macrophage", "CD8_T", "Endothelial"]


def main():
    bundles = sorted(glob.glob(os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "data", "xenium", "output-*")))
    bundles = [b for b in bundles if os.path.isdir(b)]
    fig, axes = plt.subplots(1, len(bundles), figsize=(6.2 * len(bundles), 5),
                             sharex=True)
    axes = np.atleast_1d(axes)
    for ax, b in zip(axes, bundles):
        ad = dl.load_xenium_bundle(b); labels, _ = dl.annotate_pdac(ad)
        coords = np.asarray(ad.obsm["spatial"], float)
        zs, names, cols = [], [], []
        for ct in CTS:
            if (labels == ct).sum() < 5:
                continue
            r = rim_enrichment(coords, labels, tumor="Tumor", barrier=ct,
                               shell_um=30, n_perms=500)
            if r and r["z_score"] is not None:
                zs.append(r["z_score"]); names.append(ct); cols.append(COLORS[ct])
        order = np.argsort(zs)
        zs = np.array(zs)[order]; names = np.array(names)[order]
        cols = np.array(cols)[order]
        ax.barh(range(len(names)), zs, color=cols, height=0.62)
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
        ax.axvline(0, color="black", lw=1)
        ax.axvline(2, color="#888", ls="--", lw=1)
        ax.axvline(-2, color="#888", ls="--", lw=1)
        for i, z in enumerate(zs):
            ax.text(z + (0.6 if z >= 0 else -0.6), i, f"{z:+.0f}", va="center",
                    ha="left" if z >= 0 else "right", fontsize=9, fontweight="bold")
        ax.set_title(os.path.basename(b)[14:26], fontsize=10, fontweight="bold")
        ax.set_xlabel("종양 rim(30µm) 농축 z  (>2 = 종양주변 농축, <−2 = 배제)")
    fig.suptitle("실제 PDAC Xenium — 종양 주변에 실제로 뭐가 있나 (개선된 rim 지표)\n"
                 "myCAF는 rim에서 배제, CD8·내피·apCAF가 종양 인접 = 면역 침윤형",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "rim_enrichment_pdac.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
