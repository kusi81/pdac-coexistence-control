"""SCOTIA 저자 주석으로 myCAF 양성대조 재검 (F3 부활).

우리 Xenium module-score 주석은 양성대조 실패 → SCOTIA 저자 주석(CosMx 1000-plex,
Pericyte/SMC를 CAF와 분리)으로 우리 지표를 그대로 적용.
양성대조: myCAF가 iCAF보다 Malignant에 가까워야(SCOTIA Nat Genet 2024).
시료별 proximity_test + 치료군(Untreated/CRT) 집계.
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree
from spatial_core import proximity_test, median_nn_distance

H5 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                  "data", "scotia", "raw_meta_data_final.h5ad")
PX_UM = 0.12028   # CosMx SMI 픽셀 크기(µm/px) — 스케일 참고용


def main():
    a = ad.read_h5ad(H5, backed="r")
    obs = a.obs
    # 라벨: 서브타입 기반 + Malignant는 majortype으로 확정
    sub = obs["annotation_subtypes"].astype(str).values
    maj = obs["annotation_majortypes"].astype(str).values
    lab = sub.copy()
    lab[maj == "Malignant"] = "Malignant"
    X = obs["CenterX_global_px"].astype(float).values
    Y = obs["CenterY_global_px"].astype(float).values
    samp = obs["sample_id"].astype(str).values
    treat = obs["treatment_status"].astype(str).values

    rows = []
    print(f"{'sample':<9}{'treat':<11}{'양성대조':<9}"
          f"{'myCAF':>8}{'iCAF':>8}{'diff':>7}{'p':>7}{'nnµm':>7}")
    print("-" * 66)
    for s in pd.unique(samp):
        m = samp == s
        coords = np.column_stack([X[m], Y[m]])
        labels = lab[m]
        tr = treat[m][0]
        pr = proximity_test(coords, labels, "Malignant", "myCAF", "iCAF",
                            n_perms=200)
        if not pr:
            continue
        # 스케일 참고: 전체 세포 최근접 이웃 거리(px→µm)
        d, _ = cKDTree(coords).query(coords, k=2)
        nn_um = float(np.median(d[:, 1])) * PX_UM
        ok = pr["difference"] < 0          # myCAF 더 가까우면 통과
        da_um = pr["median_dist_a"] * PX_UM
        db_um = pr["median_dist_b"] * PX_UM
        rows.append(dict(sample=s, treat=tr, ok=ok,
                         mycaf_um=da_um, icaf_um=db_um,
                         diff_um=(da_um - db_um), p=pr["p_permutation"]))
        print(f"{s:<9}{tr:<11}{'✅통과' if ok else '❌실패':<9}"
              f"{da_um:>8.1f}{db_um:>8.1f}{da_um-db_um:>+7.1f}"
              f"{pr['p_permutation']:>7.3f}{nn_um:>7.1f}")

    df = pd.DataFrame(rows)
    npass = int(df["ok"].sum())
    print(f"\n양성대조 통과: {npass}/{len(df)} 시료 "
          f"({'myCAF가 iCAF보다 종양에 가까움 = SCOTIA 재현' if npass > len(df)//2 else '미흡'})")
    print("\n치료군별 평균 (µm):")
    for tr, g in df.groupby("treat"):
        print(f"  {tr:<11} n={len(g)}  myCAF {g['mycaf_um'].mean():6.1f}  "
              f"iCAF {g['icaf_um'].mean():6.1f}  "
              f"diff {g['diff_um'].mean():+6.1f}  통과 {int(g['ok'].sum())}/{len(g)}")

    out = os.path.join(os.path.dirname(H5), "scotia_posctrl.csv")
    df.to_csv(out, index=False)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
