"""
Synthetic PDAC spatial tissue generator (import-safe).

Adapted from make_synthetic_spatial.py so it can be imported by the dashboard
without the original module's import-time os.chdir / data-dir side effects.

Two ground-truth tissues, IDENTICAL cell counts, different architecture:

  contained : tumor islets ringed by myCAF; CD8 pushed outside the ring.
              -> barrier score MUST be high, radial CD8 density MUST rise with
                 distance from tumor.
  diffuse   : same counts and same tumor islets, but myCAF scattered at random
              and CD8 free to sit anywhere.
              -> barrier score MUST be near null.

Never report numbers from synthetic tissue -- it exists only to prove the
metrics can separate architecture from abundance.
"""

import numpy as np
import pandas as pd

FIELD = 1500.0        # um
N_ISLETS = 7
ISLET_R = 90.0
RING_INNER, RING_OUTER = 95.0, 145.0

# Default population sizes (identical across modes -> abundance is not the signal)
N_MYCAF = 1100
N_ICAF = 700
N_CD8 = 450
N_MAC = 350


def _islet_centers(rng, n=N_ISLETS, field=FIELD, min_sep=330):
    pts, tries = [], 0
    while len(pts) < n and tries < 8000:
        tries += 1
        p = rng.uniform(180, field - 180, 2)
        if all(np.linalg.norm(p - q) > min_sep for q in pts):
            pts.append(p)
    return np.array(pts)


def _in_disc(rng, center, radius, n):
    t = rng.uniform(0, 2 * np.pi, n)
    r = radius * np.sqrt(rng.uniform(0, 1, n))
    return center + np.c_[r * np.cos(t), r * np.sin(t)]


def _in_annulus(rng, center, r0, r1, n):
    t = rng.uniform(0, 2 * np.pi, n)
    r = np.sqrt(rng.uniform(r0 ** 2, r1 ** 2, n))
    return center + np.c_[r * np.cos(t), r * np.sin(t)]


def make_tissue(mode="contained", seed=42):
    """Return (coords[N,2] in microns, labels[N] str, islet_centers[K,2])."""
    rng = np.random.default_rng(seed)
    centers = _islet_centers(rng)
    xy, lab = [], []

    # --- tumor: identical in both modes -------------------------------
    for c in centers:
        n = rng.integers(150, 220)
        xy.append(_in_disc(rng, c, ISLET_R, n))
        lab += ["Tumor"] * n

    if mode == "contained":
        # myCAF forms rings around islets
        per = N_MYCAF // len(centers)
        for c in centers:
            xy.append(_in_annulus(rng, c, RING_INNER, RING_OUTER, per))
            lab += ["myCAF"] * per
        # iCAF distal in stroma
        pts = rng.uniform(0, FIELD, (N_ICAF * 3, 2))
        d = np.min(np.linalg.norm(pts[:, None] - centers[None], axis=2), axis=1)
        pts = pts[d > RING_OUTER + 25][:N_ICAF]
        xy.append(pts); lab += ["iCAF"] * len(pts)
        # CD8 excluded: only outside the ring
        pts = rng.uniform(0, FIELD, (N_CD8 * 4, 2))
        d = np.min(np.linalg.norm(pts[:, None] - centers[None], axis=2), axis=1)
        pts = pts[d > RING_OUTER + 15][:N_CD8]
        xy.append(pts); lab += ["CD8_T"] * len(pts)

    elif mode == "diffuse":  # same counts, no architecture
        for name, n in [("myCAF", N_MYCAF), ("iCAF", N_ICAF), ("CD8_T", N_CD8)]:
            pts = rng.uniform(0, FIELD, (n * 3, 2))
            d = np.min(np.linalg.norm(pts[:, None] - centers[None], axis=2), axis=1)
            pts = pts[d > ISLET_R + 5][:n]      # only exclude tumor core
            xy.append(pts); lab += [name] * len(pts)
    else:
        raise ValueError("mode must be 'contained' or 'diffuse'")

    pts = rng.uniform(0, FIELD, (N_MAC * 3, 2))
    d = np.min(np.linalg.norm(pts[:, None] - centers[None], axis=2), axis=1)
    pts = pts[d > ISLET_R + 5][:N_MAC]
    xy.append(pts); lab += ["Macrophage"] * len(pts)

    coords = np.vstack(xy)
    labels = np.array(lab)
    jit = rng.normal(0, 4, coords.shape)     # measurement noise
    return coords + jit, labels, centers


# Marker panel used to build a mock expression matrix (only needed when an
# AnnData object is required, i.e. for the squidpy-based metrics 1 & 4).
GENES = ["ACTA2", "TAGLN", "POSTN", "LRRC15", "IL6", "CXCL12", "CFD",
         "KRT19", "EPCAM", "CD8A", "CD3E", "GZMB", "CD68", "CD14",
         "COL1A1", "DCN", "PDPN", "MSLN", "CD74", "HLA-DRA"]

_PROFILE = {
    "Tumor":      dict(KRT19=4, EPCAM=4),
    "myCAF":      dict(ACTA2=4, TAGLN=4, POSTN=3.5, LRRC15=3, COL1A1=3),
    "iCAF":       dict(IL6=4, CXCL12=3.5, CFD=3, PDPN=2.5, DCN=3),
    "CD8_T":      dict(CD8A=4, CD3E=4, GZMB=3),
    "Macrophage": dict(CD68=4, CD14=3.5, HLA_DRA=3),
}


def to_anndata(coords, labels, seed=1):
    """Build an AnnData with mock expression. Requires the `anndata` package."""
    import anndata as ad
    rng = np.random.default_rng(seed)
    X = rng.normal(0.5, 0.3, (len(labels), len(GENES))).clip(0)
    gi = {g: i for i, g in enumerate(GENES)}
    for ct, sig in _PROFILE.items():
        m = labels == ct
        for g, v in sig.items():
            g = g.replace("_", "-")
            if g in gi:
                X[m, gi[g]] += rng.normal(v, 0.5, m.sum()).clip(0)
    a = ad.AnnData(X=X.astype(np.float32))
    a.var_names = GENES
    a.obs["cell_type"] = pd.Categorical(labels)
    a.obsm["spatial"] = coords
    return a
