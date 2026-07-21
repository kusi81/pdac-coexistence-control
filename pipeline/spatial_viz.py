"""
Visualization for the spatial layer.

  Fig 9  tissue map        -> see the architecture; sanity-check cell typing
  Fig 10 radial profile    -> is CD8 excluded from the tumor bed
  Fig 11 barrier score     -> the containment test, observed vs null
  Fig 12 nhood enrichment  -> who neighbors whom
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from spatial_core import radial_profile

COLORS = {
    "Tumor": "#3D3D3D", "myCAF": "#C0392B", "iCAF": "#2E86AB",
    "CD8_T": "#27AE60", "Macrophage": "#E67E22", "apCAF": "#8E44AD",
}
plt.rcParams.update({"figure.dpi": 150, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})


def fig_tissue(coords, labels, title="Tissue map", outdir="results", tag=""):
    fig, ax = plt.subplots(figsize=(6.6, 6.4))
    order = ["Macrophage", "iCAF", "myCAF", "CD8_T", "Tumor"]
    present = [c for c in order if c in set(labels)] + \
              [c for c in set(labels) if c not in order]
    for ct in present:
        m = labels == ct
        ax.scatter(coords[m, 0], coords[m, 1], s=5,
                   c=COLORS.get(ct, "#999"), alpha=0.75,
                   linewidths=0, label=f"{ct} ({m.sum()})")
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, markerscale=2.4,
              loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    p = f"{outdir}/fig9_tissue_map{tag}.png"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


def fig_radial(coords, labels, tumor="Tumor", targets=("CD8_T", "myCAF", "iCAF"),
               max_dist=250, outdir="results", tag=""):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for ct in targets:
        prof = radial_profile(coords, labels, tumor, ct,
                              max_dist=max_dist, n_bins=18)
        if prof is None:
            continue
        ax.plot(prof["distance_um"], prof["density_norm"], "-o", ms=3.2,
                lw=1.6, color=COLORS.get(ct, "#999"), label=ct)
    ax.set_xlabel(f"Distance from nearest {tumor} cell (µm)")
    ax.set_ylabel("Normalized density")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Radial density profile\n"
                 "CD8 rising with distance = immune exclusion", fontsize=10)
    fig.tight_layout()
    p = f"{outdir}/fig10_radial_profile{tag}.png"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


def fig_barrier(barrier_results, outdir="results", tag=""):
    """Observed corridor density vs matched-count random-stroma null."""
    items = [(k, v) for k, v in barrier_results.items() if v]
    if not items:
        return None
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.4))

    names = [k for k, _ in items]
    obs = [v["observed_corridor_density"] for _, v in items]
    nul = [v["null_mean"] for _, v in items]
    sd = [v["null_sd"] for _, v in items]
    x = np.arange(len(names))

    ax1.bar(x - 0.19, obs, 0.38, label="observed",
            color=[COLORS.get(n, "#777") for n in names])
    ax1.bar(x + 0.19, nul, 0.38, yerr=sd, label="null (random stroma)",
            color="#BBB", capsize=3)
    ax1.set_xticks(x); ax1.set_xticklabels(names)
    ax1.set_ylabel("Corridor density (barrier cells / sample point)")
    ax1.legend(frameon=False, fontsize=8)
    ax1.set_title("Barrier mass on tumor↔T-cell paths", fontsize=10)

    z = [v["z_score"] if v["z_score"] is not None else 0 for _, v in items]
    cols = ["#C0392B" if zz > 2 else "#7F8C8D" for zz in z]
    ax2.barh(names, z, color=cols, height=0.5)
    ax2.axvline(2, color="black", ls="--", lw=1)
    ax2.text(2, -0.42, " z=2", fontsize=7.5, va="top")
    ax2.set_xlabel("z vs matched-count random stroma")
    ax2.set_title("z > 2 = interposed beyond chance\n"
                  "(null matches cell COUNT, so this is positioning)",
                  fontsize=9.5)
    fig.tight_layout()
    p = f"{outdir}/fig11_barrier_score{tag}.png"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


def fig_nhood(z, outdir="results", tag=""):
    if z is None:
        return None
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    v = np.nanmax(np.abs(z.values))
    im = ax.imshow(z.values, cmap="RdBu_r", vmin=-v, vmax=v)
    ax.set_xticks(range(len(z.columns)))
    ax.set_xticklabels(z.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(z.index)))
    ax.set_yticklabels(z.index, fontsize=8)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            ax.text(j, i, f"{z.values[i, j]:.0f}", ha="center", va="center",
                    fontsize=6.5,
                    color="white" if abs(z.values[i, j]) > v * 0.55 else "black")
    ax.set_title("Neighborhood enrichment (z)\nsymmetric — shows adjacency, "
                 "not direction", fontsize=9.5)
    fig.colorbar(im, ax=ax, fraction=0.046, label="z-score")
    fig.tight_layout()
    p = f"{outdir}/fig12_nhood_enrichment{tag}.png"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


def plot_all_spatial(adata, results, outdir="results", tag=""):
    coords = np.asarray(adata.obsm["spatial"])
    labels = np.asarray(adata.obs["cell_type"])
    made = [
        fig_tissue(coords, labels, f"Tissue map{tag}", outdir, tag),
        fig_radial(coords, labels, outdir=outdir, tag=tag),
        fig_barrier(results.get("barrier", {}), outdir, tag),
        fig_nhood(results.get("nhood_zscore"), outdir, tag),
    ]
    made = [p for p in made if p]
    for p in made:
        print(f"[INFO] wrote {p}", flush=True)
    return made
