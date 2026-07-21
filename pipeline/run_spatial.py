"""
Layer 4 spatial analysis runner.

Usage:
  python3 run_spatial.py --adata data/sample.h5ad --platform xenium --outdir results

Input: .h5ad with
  - adata.obsm['spatial']  : (n_cells, 2) coordinates IN MICRONS
  - adata.obs['cell_type'] : categorical cell/spot labels

For Visium, cell_type should be the argmax deconvolution label (RCTD /
cell2location / CARD) -- and see the platform caveat, because spot-level
argmax is a coarse approximation.
"""

import argparse, json, os
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, pd.DataFrame):
        return o.to_dict(orient="list")
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


def main():
    ap = argparse.ArgumentParser(description="Spatial CAF architecture analysis")
    ap.add_argument("--adata", required=True)
    ap.add_argument("--platform", default="xenium",
                    choices=["visium", "visium_hd", "xenium", "cosmx"])
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--celltype-key", default="cell_type")
    ap.add_argument("--tumor", default="Tumor")
    ap.add_argument("--immune", default="CD8_T")
    ap.add_argument("--barriers", default="myCAF,iCAF")
    ap.add_argument("--radius", type=float, default=50.0,
                    help="neighbor graph radius in microns")
    ap.add_argument("--n-perms", type=int, default=1000)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    import anndata as ad
    from spatial_core import run_spatial_suite, log
    from spatial_viz import plot_all_spatial

    os.makedirs(a.outdir, exist_ok=True)
    adata = ad.read_h5ad(a.adata)

    if "spatial" not in adata.obsm:
        raise ValueError("adata.obsm['spatial'] missing -- coordinates required")
    if a.celltype_key not in adata.obs:
        raise ValueError(f"adata.obs['{a.celltype_key}'] missing. "
                         "Run deconvolution/annotation first.")
    if not pd.api.types.is_categorical_dtype(adata.obs[a.celltype_key]):
        adata.obs[a.celltype_key] = adata.obs[a.celltype_key].astype("category")

    log(f"loaded {adata.n_obs} units x {adata.n_vars} genes")

    res = run_spatial_suite(
        adata, platform=a.platform, celltype_key=a.celltype_key,
        tumor=a.tumor, immune=a.immune,
        barriers=tuple(a.barriers.split(",")),
        radius=a.radius, n_perms=a.n_perms)

    # save
    if res.get("nhood_zscore") is not None:
        res["nhood_zscore"].to_csv(f"{a.outdir}/spatial_nhood_z{a.tag}.csv")
    for ct, prof in res.get("radial", {}).items():
        prof.to_csv(f"{a.outdir}/spatial_radial_{ct}{a.tag}.csv", index=False)

    summary = {k: v for k, v in res.items()
               if k not in ("nhood_zscore", "radial", "cooccurrence")}
    with open(f"{a.outdir}/spatial_results{a.tag}.json", "w") as f:
        json.dump(_jsonable(summary), f, indent=2)
    log(f"wrote {a.outdir}/spatial_results{a.tag}.json")

    plot_all_spatial(adata, res, outdir=a.outdir, tag=a.tag)


if __name__ == "__main__":
    main()
