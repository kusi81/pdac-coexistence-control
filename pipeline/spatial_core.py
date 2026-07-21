"""
Layer 4: spatial analysis of CAF architecture in PDAC.

WHAT THIS LAYER EXISTS TO ANSWER
--------------------------------
Bulk composition (Layer 3b) cannot distinguish myCAF *surrounding* tumor islets
from myCAF merely *abundant*. The containment hypothesis -- "CAFs wall the tumor
off from T cells" -- is a claim about geometry, not abundance. This module tests
it directly.

The literature anchor: myCAFs are more proximal to malignant cells than iCAFs
(Nat Genet 2024, SCOTIA study), and myCAFs surround tumor islets while iCAFs sit
at greater distance in the stroma. So a correct implementation should reproduce
myCAF-proximal / iCAF-distal as a positive control. If it does not, distrust the
cell typing before distrusting the biology.

FOUR METRICS, EACH ANSWERING A DIFFERENT QUESTION
-------------------------------------------------
1. Neighborhood enrichment  -> who sits next to whom (squidpy, permutation-based)
2. Radial distance profile  -> how does CD8 density vary with distance from tumor
3. Barrier score            -> is CAF mass positioned BETWEEN tumor and T cells
4. Co-occurrence vs radius  -> at what length scale does association appear

Metric 3 is the one that actually tests containment. Enrichment and co-occurrence
are symmetric: they say "these are near each other" without direction. A barrier
is inherently directional -- CAF must be interposed on the path from T cell to
tumor. Reporting only 1/2/4 and calling it "containment" would overclaim.

PLATFORM MATTERS
----------------
Visium spots are 55um and contain 1-10 cells -> spot-level deconvolution
proportions, NOT cell types. Distances are spot-to-spot; sub-spot geometry is
invisible. Containment claims from Visium are weak.
Xenium/CosMx are single-cell -> real cell coordinates, real distances. Only these
support a strong containment claim.
This module handles both but reports the platform's inferential ceiling.
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu

CELLTYPE_KEY = "cell_type"
SPATIAL_KEY = "spatial"


def log(msg, level="INFO"):
    import sys
    print(f"[{level}] {msg}", file=sys.stderr)


# ==========================================================================
# Platform handling
# ==========================================================================

PLATFORM_SPECS = {
    "visium": dict(resolution="spot", spot_um=55, spacing_um=100,
                   single_cell=False,
                   caveat="55um spots contain 1-10 cells. Cell 'types' are "
                          "deconvolution proportions. Sub-spot geometry invisible; "
                          "containment claims are WEAK at this resolution."),
    "visium_hd": dict(resolution="bin", spot_um=8, spacing_um=8,
                      single_cell=False,
                      caveat="8um bins approach single-cell but still binned. "
                             "Containment claims are MODERATE."),
    "xenium": dict(resolution="cell", spot_um=None, spacing_um=None,
                   single_cell=True,
                   caveat="Single-cell resolution. Containment claims supported, "
                          "but panel is targeted (~300-500 genes) so CAF subtyping "
                          "depends on panel coverage."),
    "cosmx": dict(resolution="cell", spot_um=None, spacing_um=None,
                  single_cell=True,
                  caveat="Single-cell resolution. Containment claims supported. "
                         "Check FOV edge effects -- cells near FOV borders have "
                         "truncated neighborhoods."),
}


def check_platform(platform):
    if platform not in PLATFORM_SPECS:
        raise ValueError(f"platform must be one of {list(PLATFORM_SPECS)}")
    spec = PLATFORM_SPECS[platform]
    log(f"platform={platform} resolution={spec['resolution']}")
    log(spec["caveat"], "WARN" if not spec["single_cell"] else "INFO")
    return spec


# ==========================================================================
# Metric 1: neighborhood enrichment (squidpy wrapper)
# ==========================================================================

def neighborhood_enrichment(adata, celltype_key=CELLTYPE_KEY, n_perms=1000,
                            radius=None, n_neighs=6, seed=0):
    """Permutation-based neighbor enrichment. Symmetric: no directionality.

    Uses radius-based graph when radius is given (preferred for single-cell
    platforms, since it is a real length scale), else kNN.
    """
    import squidpy as sq

    if radius is not None:
        sq.gr.spatial_neighbors(adata, radius=radius, coord_type="generic",
                                delaunay=False)
        log(f"neighbor graph: radius={radius}um")
    else:
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighs, coord_type="generic")
        log(f"neighbor graph: kNN k={n_neighs}")

    sq.gr.nhood_enrichment(adata, cluster_key=celltype_key, n_perms=n_perms,
                           seed=seed, show_progress_bar=False)

    res = adata.uns[f"{celltype_key}_nhood_enrichment"]
    cats = list(adata.obs[celltype_key].cat.categories)
    z = pd.DataFrame(res["zscore"], index=cats, columns=cats)
    return z


# ==========================================================================
# Metric 2: radial distance profile
# ==========================================================================

def radial_profile(coords, labels, source, target, max_dist=200, n_bins=20):
    """Density of `target` cells as a function of distance from `source` cells.

    This is the classic immune-exclusion readout: if CD8 density rises with
    distance from tumor, T cells are excluded from the tumor bed.
    """
    src = coords[labels == source]
    tgt = coords[labels == target]
    if len(src) < 5 or len(tgt) < 5:
        return None

    tree = cKDTree(src)
    d, _ = tree.query(tgt, k=1)
    d = d[d <= max_dist]
    if len(d) == 0:
        return None

    edges = np.linspace(0, max_dist, n_bins + 1)
    counts, _ = np.histogram(d, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2

    # Normalize by annulus area so this is a density, not a raw count.
    # Without this, counts trivially increase with radius and every tissue
    # looks "excluded".
    areas = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    density = counts / areas

    return pd.DataFrame({
        "distance_um": centers, "count": counts,
        "density": density,
        "density_norm": density / density.max() if density.max() > 0 else density,
    })


def median_nn_distance(coords, labels, source, target):
    """Median nearest-neighbor distance from each target to nearest source."""
    src = coords[labels == source]
    tgt = coords[labels == target]
    if len(src) < 3 or len(tgt) < 3:
        return np.nan, 0
    d, _ = cKDTree(src).query(tgt, k=1)
    return float(np.median(d)), len(d)


def proximity_test(coords, labels, source, target_a, target_b, n_perms=1000,
                   seed=0):
    """Is target_a closer to source than target_b? (e.g. myCAF vs iCAF -> tumor)

    Positive control for the whole module: literature says myCAF is closer to
    malignant cells than iCAF. Permutation test on the label assignment keeps
    the spatial arrangement fixed and only shuffles which cells carry which of
    the two labels, so the null preserves tissue geometry.
    """
    src = coords[labels == source]
    a = coords[labels == target_a]
    b = coords[labels == target_b]
    if len(src) < 5 or len(a) < 5 or len(b) < 5:
        return None

    tree = cKDTree(src)
    da, _ = tree.query(a, k=1)
    db, _ = tree.query(b, k=1)
    obs = float(np.median(da) - np.median(db))

    pooled = np.vstack([a, b])
    na = len(a)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perms)
    for i in range(n_perms):
        idx = rng.permutation(len(pooled))
        pa, pb = pooled[idx[:na]], pooled[idx[na:]]
        null[i] = np.median(tree.query(pa, k=1)[0]) - np.median(tree.query(pb, k=1)[0])

    p = float((np.sum(np.abs(null) >= abs(obs)) + 1) / (n_perms + 1))
    u, p_mw = mannwhitneyu(da, db, alternative="two-sided")
    return {
        "source": source, "target_a": target_a, "target_b": target_b,
        "median_dist_a": float(np.median(da)),
        "median_dist_b": float(np.median(db)),
        "difference": obs, "p_permutation": p, "p_mannwhitney": float(p_mw),
        "n_a": len(a), "n_b": len(b),
        "interpretation": (f"{target_a} closer to {source}" if obs < 0
                           else f"{target_b} closer to {source}"),
    }


# ==========================================================================
# Metric 3: BARRIER SCORE -- the actual containment test
# ==========================================================================

def barrier_score(coords, labels, tumor="Tumor", barrier="myCAF",
                  immune="CD8_T", n_samples=2000, corridor_um=30,
                  tumor_shell_um=None, seed=0):
    """Fraction of tumor-immune straight-line paths that pass through barrier cells.

    THE LOGIC
    ---------
    Containment means the barrier is INTERPOSED. For each (immune, nearest tumor)
    pair, walk the segment between them and ask whether barrier cells lie within
    `corridor_um` of that path. High fraction = physically interposed.

    RIM FOCUS (tumor_shell_um)
    --------------------------
    Real containment is a THIN barrier at the tumor rim, but averaging the corridor
    over the whole path dilutes that rim signal when immune cells are far away (the
    path is mostly intervening stroma). When `tumor_shell_um` is set, only path
    sample points within that distance of the nearest tumor cell are counted --
    i.e. the peritumoral shell the immune cell must cross to reach the tumor. This
    isolates the rim while preserving directionality (still along the immune->tumor
    path). None = legacy whole-path behavior.

    THE NULL
    --------
    Comparing the observed fraction to nothing is meaningless -- in a tissue that
    is 80% CAF, most paths hit a CAF by chance. The null here relabels a random
    subset of non-tumor, non-immune cells as pseudo-barrier, preserving both the
    tissue geometry and the barrier cell COUNT. So the score answers: are these
    specific cells positioned as a barrier more than an equally numerous random
    stromal population would be?

    LIMITATION
    ----------
    Straight-line paths are a simplification; real T cell migration follows ECM
    fiber tracks and vessels. This is a geometric proxy, not a migration model.
    """
    tum = coords[labels == tumor]
    bar = coords[labels == barrier]
    imm = coords[labels == immune]
    if len(tum) < 5 or len(bar) < 5 or len(imm) < 5:
        log(f"insufficient cells for barrier score "
            f"(tumor={len(tum)}, barrier={len(bar)}, immune={len(imm)})", "WARN")
        return None

    rng = np.random.default_rng(seed)
    if len(imm) > n_samples:
        imm = imm[rng.choice(len(imm), n_samples, replace=False)]

    tumor_tree = cKDTree(tum)
    _, ti = tumor_tree.query(imm, k=1)
    partners = tum[ti]

    def _corridor_density(barrier_pts):
        """Mean barrier cells per corridor sample point along tumor-immune paths.

        Continuous, not binary. An earlier binary 'is any point blocked'
        formulation saturated at ~1.0 in dense stroma and had zero
        discriminative power -- in an 80% stromal tissue every path touches
        some CAF. Counting barrier MASS in the corridor keeps the signal.
        With tumor_shell_um, only peritumoral-shell sample points are averaged.
        """
        tree = cKDTree(barrier_pts)
        vals = []
        for p0, p1 in zip(imm, partners):
            seg = p1 - p0
            L = np.linalg.norm(seg)
            if L < 1e-9:
                continue
            n_steps = max(3, int(L / (corridor_um / 2)))
            ts = np.linspace(0.15, 0.85, n_steps)   # exclude endpoints
            pts = p0 + ts[:, None] * seg
            if tumor_shell_um is not None:           # rim focus
                dts, _ = tumor_tree.query(pts, k=1)
                pts = pts[dts <= tumor_shell_um]
                if len(pts) == 0:
                    continue
            n_near = tree.query_ball_point(pts, r=corridor_um,
                                           return_length=True)
            vals.append(float(np.mean(n_near)))
        return float(np.mean(vals)) if vals else np.nan

    observed = _corridor_density(bar)

    # Null: equally many random stromal cells. Preserves tissue geometry and
    # barrier cell COUNT, so the test asks about POSITIONING, not abundance.
    other_mask = ~np.isin(labels, [tumor, immune])
    other = coords[other_mask]
    nulls = []
    if len(other) >= len(bar):
        for i in range(20):
            idx = rng.choice(len(other), len(bar), replace=False)
            nulls.append(_corridor_density(other[idx]))
    null_mean = float(np.mean(nulls)) if nulls else np.nan
    null_sd = float(np.std(nulls)) if nulls else np.nan
    z = (observed - null_mean) / null_sd if (nulls and null_sd > 0) else np.nan
    p = ((np.sum(np.array(nulls) >= observed) + 1) / (len(nulls) + 1)) if nulls else np.nan

    return {
        "barrier": barrier, "tumor": tumor, "immune": immune,
        "observed_corridor_density": float(observed),
        "null_mean": null_mean, "null_sd": null_sd,
        "enrichment_ratio": float(observed / null_mean) if null_mean else None,
        "z_score": float(z) if not np.isnan(z) else None,
        "p_value": float(p) if not np.isnan(p) else None,
        "n_immune_sampled": int(len(imm)), "n_barrier": int(len(bar)),
        "corridor_um": corridor_um, "tumor_shell_um": tumor_shell_um,
        "interpretation": (
            "barrier interposed beyond chance" if (not np.isnan(z) and z > 2)
            else "no evidence of interposition beyond stromal density"),
    }


def rim_enrichment(coords, labels, tumor="Tumor", barrier="myCAF",
                   shell_um=30.0, n_perms=1000, seed=0):
    """Is the barrier type enriched in the peritumoral SHELL (rim), vs a matched
    random-stroma null? A direct, path-free test of 'surrounds the tumor'.

    Rim = non-tumor cells within `shell_um` of the nearest tumor cell. We compare
    the barrier fraction inside the rim to a null that relabels an equal COUNT of
    random non-tumor cells (so it tests positioning, not abundance) -- and also to
    the barrier's overall fraction (enrichment ratio). Robust to thin rims and
    immune-cell distance, unlike the straight-line path metric.
    """
    tum = coords[labels == tumor]
    nontum_mask = labels != tumor
    other = coords[nontum_mask]
    is_bar = (labels[nontum_mask] == barrier)
    if len(tum) < 5 or is_bar.sum() < 5 or len(other) < 10:
        return None

    d, _ = cKDTree(tum).query(other, k=1)
    in_shell = d <= shell_um
    n_shell = int(in_shell.sum())
    if n_shell < 5:
        return None

    obs_barrier_in_shell = int((is_bar & in_shell).sum())
    obs_frac = obs_barrier_in_shell / n_shell
    overall_frac = float(is_bar.mean())          # barrier share of all non-tumor

    # Null: shuffle which non-tumor cells are 'barrier' (matched count), recount
    # how many land in the shell. Tests positioning of THESE cells in the rim.
    rng = np.random.default_rng(seed)
    nb = int(is_bar.sum())
    shell_idx = np.where(in_shell)[0]
    null = np.empty(n_perms)
    N = len(other)
    for i in range(n_perms):
        pick = rng.choice(N, nb, replace=False)
        null[i] = np.isin(pick, shell_idx).sum() / n_shell
    null_mean = float(null.mean()); null_sd = float(null.std())
    z = (obs_frac - null_mean) / null_sd if null_sd > 0 else np.nan
    p = float((np.sum(null >= obs_frac) + 1) / (n_perms + 1))

    return {
        "barrier": barrier, "tumor": tumor, "shell_um": shell_um,
        "n_shell": n_shell, "barrier_in_shell": obs_barrier_in_shell,
        "obs_rim_fraction": float(obs_frac), "overall_fraction": overall_frac,
        "null_mean": null_mean, "null_sd": null_sd,
        "enrichment_ratio": float(obs_frac / overall_frac) if overall_frac else None,
        "z_score": float(z) if not np.isnan(z) else None, "p_value": p,
        "interpretation": (
            "barrier enriched in tumor rim beyond chance"
            if (not np.isnan(z) and z > 2)
            else "no rim enrichment beyond stromal density"),
    }


# ==========================================================================
# Metric 4: co-occurrence across length scales
# ==========================================================================

def cooccurrence(adata, celltype_key=CELLTYPE_KEY, interval=None, max_dist=300):
    """Co-occurrence probability ratio vs radius. Reveals the length scale of
    association, which a single-radius enrichment test hides."""
    import squidpy as sq
    if interval is None:
        interval = np.linspace(10, max_dist, 25)
    sq.gr.co_occurrence(adata, cluster_key=celltype_key, interval=interval,
                        show_progress_bar=False)
    res = adata.uns[f"{celltype_key}_co_occurrence"]
    return res["occ"], res["interval"], list(adata.obs[celltype_key].cat.categories)


# ==========================================================================
# Orchestration
# ==========================================================================

def run_spatial_suite(adata, platform="xenium", celltype_key=CELLTYPE_KEY,
                      tumor="Tumor", immune="CD8_T",
                      barriers=("myCAF", "iCAF"), radius=50, n_perms=1000,
                      seed=0):
    """Run all four metrics. Returns a results dict."""
    spec = check_platform(platform)
    coords = np.asarray(adata.obsm[SPATIAL_KEY])
    labels = np.asarray(adata.obs[celltype_key])

    counts = pd.Series(labels).value_counts()
    log(f"cell types: {dict(counts)}")

    out = {"platform": platform, "platform_caveat": spec["caveat"],
           "single_cell": spec["single_cell"], "cell_counts": counts.to_dict()}

    log("--- metric 1: neighborhood enrichment ---")
    try:
        out["nhood_zscore"] = neighborhood_enrichment(
            adata, celltype_key, n_perms=n_perms,
            radius=radius if spec["single_cell"] else None, seed=seed)
    except Exception as e:
        log(f"neighborhood enrichment failed: {e}", "WARN")
        out["nhood_zscore"] = None

    log("--- metric 2: radial profiles + proximity ---")
    out["radial"] = {}
    for ct in list(barriers) + [immune]:
        prof = radial_profile(coords, labels, tumor, ct)
        if prof is not None:
            out["radial"][ct] = prof
    if len(barriers) >= 2:
        pt = proximity_test(coords, labels, tumor, barriers[0], barriers[1],
                            n_perms=n_perms, seed=seed)
        out["proximity_test"] = pt
        if pt:
            log(f"{barriers[0]} median dist to {tumor}: {pt['median_dist_a']:.1f}um | "
                f"{barriers[1]}: {pt['median_dist_b']:.1f}um | p={pt['p_permutation']:.3g}")
            log(f"positive control: {pt['interpretation']} "
                f"(literature expects {barriers[0]} closer)")

    log("--- metric 3: barrier score (containment test) ---")
    out["barrier"] = {}
    for b in barriers:
        bs = barrier_score(coords, labels, tumor=tumor, barrier=b,
                           immune=immune, seed=seed)
        if bs:
            out["barrier"][b] = bs
            zz = bs["z_score"]
            log(f"{b}: corridor={bs['observed_corridor_density']:.3f} "
                f"null={bs['null_mean']:.3f} z={zz if zz is None else round(zz,2)} "
                f"-> {bs['interpretation']}")

    log("--- metric 4: co-occurrence vs radius ---")
    try:
        occ, interval, cats = cooccurrence(adata, celltype_key)
        out["cooccurrence"] = {"occ": occ, "interval": interval, "categories": cats}
    except Exception as e:
        log(f"co-occurrence failed: {e}", "WARN")
        out["cooccurrence"] = None

    return out
