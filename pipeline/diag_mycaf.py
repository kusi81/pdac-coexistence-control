"""myCAF 양성대조 실패 진단.

용의자: ① 주석(패널 커버리지·서브타입 분리) ② Tumor 경계 ③ rim 지표.
측정:
  A. Xenium 패널에 각 모듈 마커가 존재하는가 (커버리지) — 핵심
  B. 라벨 조성 (Tumor/CAF/myCAF/iCAF/apCAF 비율)
  C. proximity_test(myCAF vs iCAF → Tumor) = 모듈의 공식 양성대조
  D. median_nn 거리: myCAF/iCAF/apCAF → 가장 가까운 Tumor
두 시료: 31076(myCAF rim +9, '되는 편') vs 28429(myCAF rim -40, '최악').
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np
import data_loader as dl
from data_loader import (PDAC_COARSE_MODULES, CAF_SUBTYPE_MODULES)
from spatial_core import proximity_test, median_nn_distance

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAMPLES = {
    "31076": "output-XETG00248__0014632__31076__20240529__214643",
    "28429": "output-XETG00248__0014635__28429__20240529__214643",
}


def coverage_report(adata):
    panel = {g.upper() for g in map(str, adata.var_names)}
    print(f"  패널 유전자 수: {len(panel)}")
    for name, mods in [("COARSE", PDAC_COARSE_MODULES),
                       ("CAF_SUBTYPE", CAF_SUBTYPE_MODULES)]:
        print(f"  [{name}]")
        for t, genes in mods.items():
            present = [g for g in genes if g.upper() in panel]
            absent = [g for g in genes if g.upper() not in panel]
            print(f"    {t:12s} {len(present)}/{len(genes)} present"
                  f"  ✅{present}  ❌{absent}")


def main():
    first = True
    for code, folder in SAMPLES.items():
        b = os.path.join(ROOT, "data", "xenium", folder)
        print(f"\n{'='*70}\n[{code}] {folder}")
        ad = dl.load_xenium_bundle(b)
        if first:
            print("\n--- A. 패널 커버리지 (시료 무관, 한 번만) ---")
            coverage_report(ad)
            first = False
        labels, cov = dl.annotate_pdac(ad)
        coords = np.asarray(ad.obsm["spatial"], float)

        print("\n--- B. 라벨 조성 ---")
        uniq, cnt = np.unique(labels, return_counts=True)
        tot = len(labels)
        for u, c in sorted(zip(uniq, cnt), key=lambda x: -x[1]):
            print(f"    {u:12s} {c:8d}  ({100*c/tot:5.1f}%)")

        print("\n--- C. proximity_test (양성대조: myCAF가 iCAF보다 Tumor에 가까워야) ---")
        pr = proximity_test(coords, labels, source="Tumor",
                            target_a="myCAF", target_b="iCAF", n_perms=300)
        if pr:
            print(f"    myCAF med dist→Tumor: {pr['median_dist_a']:.1f}")
            print(f"    iCAF  med dist→Tumor: {pr['median_dist_b']:.1f}")
            print(f"    diff(a-b): {pr['difference']:+.1f}  "
                  f"p_perm={pr['p_permutation']:.3f}  → {pr['interpretation']}")
            ok = pr['difference'] < 0
            print(f"    양성대조 {'✅통과' if ok else '❌실패(iCAF가 더 가까움)'}")
        else:
            print("    (세포 수 부족)")

        print("\n--- D. median_nn 거리 → 가장 가까운 Tumor (작을수록 인접) ---")
        for ct in ["myCAF", "iCAF", "apCAF", "CD8_T", "Macrophage"]:
            d, n = median_nn_distance(coords, labels, source="Tumor", target=ct)
            print(f"    {ct:12s} med={d:7.1f}µm  (n={n})")


if __name__ == "__main__":
    main()
