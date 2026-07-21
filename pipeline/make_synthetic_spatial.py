"""
SYNTHETIC spatial tissue. Validation only -- never report numbers from this.

Generates two tissues with KNOWN ground truth:

  CONTAINED:  tumor islets ringed by myCAF; CD8 pushed outside the ring.
              -> barrier score MUST be high, radial CD8 density MUST rise with
                 distance from tumor.
  DIFFUSE:    same cell counts and same tumor islets, but myCAF scattered at
              random through the stroma and CD8 free to sit anywhere.
              -> barrier score MUST be near null.

If the module cannot separate these two, it does not measure containment.
The DIFFUSE tissue is the important one: it has the SAME myCAF abundance, so
any method that only measures composition will call the two identical.
"""
import os
import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs("data", exist_ok=True)

FIELD = 1500.0        # um
N_ISLETS = 7
ISLET_R = 90.0
RING_INNER, RING_OUTER = 95.0, 145.0


def _islet_centers(rng, n=N_ISLETS, field=FIELD, min_sep=330):
    pts = []
    tries = 0
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


def make_tissue(mode="contained", seed=0):
    rng = np.random.default_rng(seed)
    centers = _islet_centers(rng)
    xy, lab = [], []

    # --- tumor: identical in both modes -------------------------------
    for c in centers:
        n = rng.integers(150, 220)
        xy.append(_in_disc(rng, c, ISLET_R, n)); lab += ["Tumor"] * n

    n_my = 1100
    n_ic = 700
    n_cd8 = 450
    n_mac = 350

    if mode == "contained":
        # myCAF forms rings around islets
        per = n_my // len(centers)
        for c in centers:
            xy.append(_in_annulus(rng, c, RING_INNER, RING_OUTER, per))
            lab += ["myCAF"] * per
        # iCAF distal in stroma
        pts = rng.uniform(0, FIELD, (n_ic * 3, 2))
        d = np.min(np.linalg.norm(pts[:, None] - centers[None], axis=2), axis=1)
        pts = pts[d > RING_OUTER + 25][:n_ic]
        xy.append(pts); lab += ["iCAF"] * len(pts)
        # CD8 excluded: only outside the ring
        pts = rng.uniform(0, FIELD, (n_cd8 * 4, 2))
        d = np.min(np.linalg.norm(pts[:, None] - centers[None], axis=2), axis=1)
        pts = pts[d > RING_OUTER + 15][:n_cd8]
        xy.append(pts); lab += ["CD8_T"] * len(pts)

    else:  # diffuse -- same counts, no architecture
        for name, n in [("myCAF", n_my), ("iCAF", n_ic), ("CD8_T", n_cd8)]:
            pts = rng.uniform(0, FIELD, (n * 3, 2))
            d = np.min(np.linalg.norm(pts[:, None] - centers[None], axis=2), axis=1)
            pts = pts[d > ISLET_R + 5][:n]      # only exclude tumor core
            xy.append(pts); lab += [name] * len(pts)

    pts = rng.uniform(0, FIELD, (n_mac * 3, 2))
    d = np.min(np.linalg.norm(pts[:, None] - centers[None], axis=2), axis=1)
    pts = pts[d > ISLET_R + 5][:n_mac]
    xy.append(pts); lab += ["Macrophage"] * len(pts)

    coords = np.vstack(xy)
    labels = np.array(lab)
    jit = rng.normal(0, 4, coords.shape)
    return coords + jit, labels, centers


def to_anndata(coords, labels):
    import anndata as ad
    rng = np.random.default_rng(1)
    genes = ["ACTA2", "TAGLN", "POSTN", "LRRC15", "IL6", "CXCL12", "CFD",
             "KRT19", "EPCAM", "CD8A", "CD3E", "GZMB", "CD68", "CD14",
             "COL1A1", "DCN", "PDPN", "MSLN", "CD74", "HLA-DRA"]
    prof = {
        "Tumor":      dict(KRT19=4, EPCAM=4),
        "myCAF":      dict(ACTA2=4, TAGLN=4, POSTN=3.5, LRRC15=3, COL1A1=3),
        "iCAF":       dict(IL6=4, CXCL12=3.5, CFD=3, PDPN=2.5, DCN=3),
        "CD8_T":      dict(CD8A=4, CD3E=4, GZMB=3),
        "Macrophage": dict(CD68=4, CD14=3.5, HLA_DRA=3),
    }
    X = rng.normal(0.5, 0.3, (len(labels), len(genes))).clip(0)
    gi = {g: i for i, g in enumerate(genes)}
    for ct, sig in prof.items():
        m = labels == ct
        for g, v in sig.items():
            g = g.replace("_", "-")
            if g in gi:
                X[m, gi[g]] += rng.normal(v, 0.5, m.sum()).clip(0)

    a = ad.AnnData(X=X.astype(np.float32))
    a.var_names = genes
    a.obs["cell_type"] = pd.Categorical(labels)
    a.obsm["spatial"] = coords
    return a


if __name__ == "__main__":
    for mode in ["contained", "diffuse"]:
        coords, labels, centers = make_tissue(mode, seed=42)
        a = to_anndata(coords, labels)
        a.write(f"data/SYN_SPATIAL_{mode}.h5ad")
        np.save(f"data/SYN_SPATIAL_{mode}_centers.npy", centers)
        vc = pd.Series(labels).value_counts().to_dict()
        print(f"{mode}: {len(labels)} cells  {vc}")
    print("\nGround truth: CONTAINED has myCAF rings + CD8 exclusion;")
    print("DIFFUSE has identical cell counts but no architecture.")
