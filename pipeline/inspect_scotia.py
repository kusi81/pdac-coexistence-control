"""SCOTIA raw_meta_data_final.h5ad 구조 조사 — 저자 주석/좌표/치료상태 컬럼 특정."""
import sys, os, warnings
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import anndata as ad
import numpy as np

H5 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                  "data", "scotia", "raw_meta_data_final.h5ad")


def main():
    a = ad.read_h5ad(H5, backed="r")   # 메타만, X는 디스크
    print("shape:", a.shape)
    print("\n=== obs columns ===")
    for c in a.obs.columns:
        s = a.obs[c]
        dt = str(s.dtype)
        info = ""
        if dt == "category" or s.nunique() <= 40:
            vc = s.value_counts().head(15)
            info = " | " + ", ".join(f"{k}:{v}" for k, v in vc.items())
        else:
            try:
                info = f" | numeric range [{np.nanmin(s.values):.2f}, {np.nanmax(s.values):.2f}]"
            except Exception:
                info = f" | nunique={s.nunique()}"
        print(f"  {c:28s} [{dt}]{info}")
    print("\n=== obsm keys ===", list(a.obsm.keys()))
    for k in a.obsm.keys():
        try:
            print(f"  {k}: shape {a.obsm[k].shape}")
        except Exception:
            pass
    print("\n=== obs head ===")
    print(a.obs.head(6).to_string())


if __name__ == "__main__":
    main()
