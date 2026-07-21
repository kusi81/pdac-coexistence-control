"""
Plot builders that RETURN matplotlib Figures (for Streamlit st.pyplot),
rather than saving to disk like the original spatial_viz.py.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa: F401  (import 부작용으로 한글 폰트 등록)
import matplotlib.pyplot as plt

from spatial_core import radial_profile

COLORS = {
    "Tumor": "#3D3D3D", "myCAF": "#C0392B", "iCAF": "#2E86AB",
    "CD8_T": "#27AE60", "Macrophage": "#E67E22", "apCAF": "#8E44AD",
    "Endothelial": "#16A085", "Epithelial": "#3D3D3D", "Fibroblast": "#C0392B",
    "Tcell_CD8": "#27AE60", "Tcell_CD4": "#2ECC71", "Imm_other": "#95A5A6",
    "Myeloid_CD68": "#E67E22", "Myeloid_CD11c": "#F39C12", "B": "#9B59B6",
}
plt.rcParams.update({"figure.dpi": 130, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})


def fig_tissue(coords, labels, title="Tissue map"):
    fig, ax = plt.subplots(figsize=(6.4, 6.2))
    order = ["Macrophage", "iCAF", "myCAF", "CD8_T", "Tumor"]
    present = [c for c in order if c in set(labels)] + \
              [c for c in set(labels) if c not in order]
    for ct in present:
        m = labels == ct
        ax.scatter(coords[m, 0], coords[m, 1], s=5,
                   c=COLORS.get(ct, "#999"), alpha=0.75,
                   linewidths=0, label=f"{ct} ({int(m.sum())})")
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, markerscale=2.4,
              loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    return fig


def fig_radial(coords, labels, tumor="Tumor",
               targets=("CD8_T", "myCAF", "iCAF"), max_dist=250):
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    plotted = False
    for ct in targets:
        prof = radial_profile(coords, labels, tumor, ct,
                              max_dist=max_dist, n_bins=18)
        if prof is None:
            continue
        ax.plot(prof["distance_um"], prof["density_norm"], "-o", ms=3.2,
                lw=1.6, color=COLORS.get(ct, "#999"), label=ct)
        plotted = True
    ax.set_xlabel(f"Distance from nearest {tumor} cell (µm)")
    ax.set_ylabel("Normalized density")
    if plotted:
        ax.legend(frameon=False, fontsize=8)
    ax.set_title("Radial density profile — CD8 rising with distance = immune exclusion",
                 fontsize=9.5)
    fig.tight_layout()
    return fig


def fig_barrier(barrier_results):
    items = [(k, v) for k, v in barrier_results.items() if v]
    if not items:
        return None
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.2))

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
    ax1.set_ylabel("Corridor density (barrier cells / sample pt)")
    ax1.legend(frameon=False, fontsize=8)
    ax1.set_title("Barrier mass on tumor↔T-cell paths", fontsize=10)

    z = [v["z_score"] if v["z_score"] is not None else 0 for _, v in items]
    cols = ["#C0392B" if zz > 2 else "#7F8C8D" for zz in z]
    ax2.barh(names, z, color=cols, height=0.5)
    ax2.axvline(2, color="black", ls="--", lw=1)
    ax2.text(2, -0.42, " z=2", fontsize=7.5, va="top")
    ax2.set_xlabel("z vs matched-count random stroma")
    ax2.set_title("z > 2 = interposed beyond chance\n(null matches cell COUNT — this is positioning)",
                  fontsize=9)
    fig.tight_layout()
    return fig


def fig_nhood(z):
    if z is None:
        return None
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    v = np.nanmax(np.abs(z.values)) or 1.0
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
    ax.set_title("Neighborhood enrichment (z)\nsymmetric — adjacency, not direction",
                 fontsize=9.5)
    fig.colorbar(im, ax=ax, fraction=0.046, label="z-score")
    fig.tight_layout()
    return fig


def fig_rim(rim_rows, shell_um=30, tumor="Tumor"):
    """종양-rim 농축 z를 세포타입별 가로 막대로 (누가 juxtatumoral인가)."""
    if not rim_rows:
        return None
    rows = sorted(rim_rows, key=lambda d: d["z"])   # 아래가 배제, 위가 농축
    names = [r["cell_type"] for r in rows]
    zs = [r["z"] for r in rows]
    cols = [COLORS.get(n, "#999") for n in names]
    fig, ax = plt.subplots(figsize=(6.8, 0.5 * len(names) + 1.4))
    ax.barh(range(len(names)), zs, color=cols, height=0.62)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
    ax.axvline(0, color="black", lw=1)
    ax.axvline(2, color="#888", ls="--", lw=1)
    ax.axvline(-2, color="#888", ls="--", lw=1)
    for i, z in enumerate(zs):
        ax.text(z + (0.4 if z >= 0 else -0.4), i, f"{z:+.0f}", va="center",
                ha="left" if z >= 0 else "right", fontsize=9, fontweight="bold")
    ax.set_xlabel(f"{tumor} 주변 rim({shell_um:g}) 농축 z  "
                  "(>2 종양인접, <−2 배제)")
    ax.set_title(f"누가 {tumor}에 인접한가 — rim 농축", fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig


def fig_cooccurrence(cooc, focus=("Tumor", "CD8_T")):
    """occ has shape (n_cats, n_cats, n_intervals-1). Plot P(focus[1]|focus[0])
    ratio vs distance."""
    if not cooc:
        return None
    occ = np.asarray(cooc["occ"])
    interval = np.asarray(cooc["interval"])
    cats = list(cooc["categories"])
    if focus[0] not in cats or focus[1] not in cats:
        return None
    i, j = cats.index(focus[0]), cats.index(focus[1])
    centers = (interval[:-1] + interval[1:]) / 2
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(centers, occ[i, j], "-o", ms=3, lw=1.6,
            color=COLORS.get(focus[1], "#333"))
    ax.axhline(1.0, color="#999", ls="--", lw=1)
    ax.set_xlabel("Distance (µm)")
    ax.set_ylabel(f"P({focus[1]} | {focus[0]}) / P({focus[1]})")
    ax.set_title(f"Co-occurrence vs radius — {focus[0]} ↔ {focus[1]}", fontsize=9.5)
    fig.tight_layout()
    return fig
