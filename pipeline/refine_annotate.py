"""빠른 트랙: 주석 정밀화 → myCAF 양성대조 재검 (다운로드 0).

수정 2가지:
  1) 특이 myCAF 마커 추가 — POSTN·LRRC15·CTHRC1(골드스탠다드)를 서브타입 모듈에.
  2) Perivascular(혈관주위/평활근) 대분류 신설 — RGS5·NOTCH3·MYH11·CSPG4·DES·PLN.
     ACTA2/TAGLN이 잘못 포섭하던 세포를 CAF에서 분리 → myCAF 특이도↑.

검증: proximity_test(myCAF vs iCAF → Tumor)가 통과(myCAF가 더 가까움)로 뒤집히는가.
성공 시 data_loader에 반영.
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
from spatial_core import proximity_test, median_nn_distance

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAMPLES = {
    "31076": ("naive", "output-XETG00248__0014632__31076__20240529__214643"),
    "28429": ("naive", "output-XETG00248__0014635__28429__20240529__214643"),
    "38245": ("CRT",   "output-XETG00248__0014632__38245__20240529__214643"),
    "35406": ("CRT",   "output-XETG00248__0014635__35406__20240529__214643"),
    "39928": ("naive", "output-XETG00248__0014632__39928__20240529__214643"),
}

# ── 정밀화 모듈 ──
REFINED_COARSE = {
    "Tumor":       ["EPCAM", "KRT19", "KRT8", "KRT18", "KRT17", "TFF1", "MSLN"],
    "CD8_T":       ["CD8A", "CD3D", "CD3E", "GZMB", "GZMK", "NKG7", "CCL5"],
    "Macrophage":  ["CD68", "CD14", "ITGAM", "AIF1", "C1QA", "C1QB", "CD163"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Perivascular": ["RGS5", "NOTCH3", "MYH11", "CSPG4", "DES", "PLN"],  # 신설
    "CAF":         ["FAP", "LUM", "DCN", "COL1A2", "COL3A1", "COL5A1", "PDGFRA"],
}
REFINED_SUBTYPE = {
    # 특이 ECM myCAF 마커 중심 (ACTA2/TAGLN 의존↓)
    "myCAF": ["POSTN", "LRRC15", "CTHRC1", "MMP11", "THBS2", "COL10A1",
              "COL11A1", "GREM1", "TAGLN"],
    "iCAF":  ["IL6", "CXCL12", "CXCL1", "CXCL2", "CCL2", "CFD", "DPT",
              "C3", "C7", "CLU", "APOD"],
    "apCAF": ["CD74", "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "SLPI", "HLA-DRB1"],
}


def refined_annotate(adata):
    cs, ctypes, _ = dl._module_scores(adata, REFINED_COARSE)
    coarse = np.array([ctypes[k] for k in np.argmax(cs, axis=1)])
    labels = coarse.copy()
    caf = coarse == "CAF"
    if caf.any():
        ss, stypes, _ = dl._module_scores(adata, REFINED_SUBTYPE)
        labels[caf] = np.array([stypes[k] for k in np.argmax(ss[caf], axis=1)])
    return labels


def main():
    print(f"{'시료':<8}{'군':<6}{'양성대조':<10}{'myCAF µm':>9}{'iCAF µm':>9}"
          f"{'차이':>7}{'p':>7}")
    print("-" * 60)
    results = {}
    for code, (treat, folder) in SAMPLES.items():
        b = os.path.join(ROOT, "data", "xenium", folder)
        if not os.path.isdir(b):
            continue
        ad = dl.load_xenium_bundle(b)
        labels = refined_annotate(ad)
        coords = np.asarray(ad.obsm["spatial"], float)
        pr = proximity_test(coords, labels, "Tumor", "myCAF", "iCAF", n_perms=300)
        if not pr:
            print(f"{code:<8}{treat:<6}(세포 부족)")
            continue
        ok = pr["difference"] < 0   # myCAF가 더 가까우면 통과
        # 조성
        uniq, cnt = np.unique(labels, return_counts=True)
        comp = {u: c for u, c in zip(uniq, cnt)}
        tot = len(labels)
        results[code] = dict(treat=treat, pr=pr, ok=ok, comp=comp, tot=tot)
        print(f"{code:<8}{treat:<6}{'✅통과' if ok else '❌실패':<10}"
              f"{pr['median_dist_a']:>9.1f}{pr['median_dist_b']:>9.1f}"
              f"{pr['difference']:>+7.1f}{pr['p_permutation']:>7.3f}")

    print("\n=== 조성(%) — Perivascular 분리 효과 확인 ===")
    cts = ["Tumor", "myCAF", "iCAF", "apCAF", "Perivascular", "Endothelial",
           "Macrophage", "CD8_T"]
    hdr = "시료    " + "".join(f"{c[:5]:>8}" for c in cts)
    print(hdr)
    for code, r in results.items():
        row = f"{code:<8}"
        for c in cts:
            pct = 100 * r["comp"].get(c, 0) / r["tot"]
            row += f"{pct:>8.1f}"
        print(row)

    npass = sum(1 for r in results.values() if r["ok"])
    print(f"\n양성대조 통과: {npass}/{len(results)} 시료")


if __name__ == "__main__":
    main()
