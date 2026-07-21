"""480 Xenium 패널에 CAF 서브타입/혈관주위 정밀화 후보 마커가 존재하는지 조사."""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import data_loader as dl

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B = os.path.join(ROOT, "data", "xenium",
                 "output-XETG00248__0014632__31076__20240529__214643")

CANDIDATES = {
    "myCAF_specific": ["POSTN", "LRRC15", "CTHRC1", "MMP11", "TAGLN", "ACTA2",
                       "THBS2", "INHBA", "COL10A1", "COL11A1", "GREM1", "SULF1",
                       "FN1", "COMP", "MMP14", "TNC"],
    "iCAF_specific": ["IL6", "CXCL12", "CXCL1", "CXCL2", "CCL2", "PDGFRA", "CFD",
                      "DPT", "LMNA", "HAS1", "HAS2", "C3", "C7", "PLA2G2A",
                      "CLU", "APOD"],
    "apCAF": ["CD74", "HLA-DRA", "HLA-DPA1", "HLA-DQA1", "SLPI", "HLA-DRB1"],
    "perivascular_SMC": ["RGS5", "NOTCH3", "MYH11", "PDGFRB", "CSPG4", "DES",
                         "MCAM", "KCNJ8", "ACTA2", "MYL9", "PLN", "CNN1"],
    "panCAF": ["PDGFRA", "PDGFRB", "LUM", "DCN", "FAP", "COL1A1", "COL1A2",
               "COL3A1", "VIM", "COL6A1", "COL6A2", "COL5A1"],
}


def main():
    ad = dl.load_xenium_bundle(B)
    panel = {g.upper() for g in map(str, ad.var_names)}
    print(f"패널 유전자 수: {len(panel)}\n")
    for grp, genes in CANDIDATES.items():
        present = [g for g in genes if g.upper() in panel]
        absent = [g for g in genes if g.upper() not in panel]
        print(f"[{grp}]  {len(present)}/{len(genes)}")
        print(f"   ✅ {present}")
        print(f"   ❌ {absent}\n")


if __name__ == "__main__":
    main()
