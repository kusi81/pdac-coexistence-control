"""개선 지표 검증: 기존 whole-path barrier vs rim-focus barrier vs rim-enrichment.

① 합성(contained/diffuse): rim 지표가 여전히 contained myCAF를 잡고 diffuse는 null인가.
② 실제 Xenium(가진 2시료): rim 초점이 실데이터 myCAF containment를 살리는가.
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from synthetic import make_tissue
from spatial_core import barrier_score, rim_enrichment
import data_loader as dl


def _z(bs):
    return None if not bs or bs.get("z_score") is None else round(bs["z_score"], 1)


def row(coords, labels, tag, barrier, tumor, immune, shell):
    old = barrier_score(coords, labels, tumor=tumor, barrier=barrier, immune=immune,
                        corridor_um=shell, seed=0)
    rim = barrier_score(coords, labels, tumor=tumor, barrier=barrier, immune=immune,
                        corridor_um=shell, tumor_shell_um=shell * 1.3, seed=0)
    enr = rim_enrichment(coords, labels, tumor=tumor, barrier=barrier,
                         shell_um=shell, n_perms=500, seed=0)
    print(f"  {tag:<26} {barrier:<6}  whole-path z={_z(old)!s:>6}  "
          f"rim-focus z={_z(rim)!s:>6}  rim-enrich z={_z(enr)!s:>6}")


print("=== ① 합성 검증 (shell=30µm) ===")
for mode in ["contained", "diffuse"]:
    c, l, _ = make_tissue(mode, seed=42)
    row(c, l, f"synth-{mode}", "myCAF", "Tumor", "CD8_T", 30)
    row(c, l, f"synth-{mode}", "iCAF", "Tumor", "CD8_T", 30)

print("\n=== ② 실제 Xenium PDAC (shell=30µm, 저자 주석) ===")
bundles = sorted(glob.glob(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "xenium", "output-*")))
bundles = [b for b in bundles if os.path.isdir(b)]
for b in bundles:
    ad = dl.load_xenium_bundle(b)
    labels, _ = dl.annotate_pdac(ad)
    coords = np.asarray(ad.obsm["spatial"], float)
    tag = os.path.basename(b)[14:26]
    row(coords, labels, f"xen-{tag}", "myCAF", "Tumor", "CD8_T", 30)
    row(coords, labels, f"xen-{tag}", "iCAF", "Tumor", "CD8_T", 30)
print("\n기대: contained myCAF는 세 지표 모두 양수, diffuse·iCAF는 낮음.")
print("실데이터: rim 지표에서 myCAF가 whole-path보다 개선되면 성공.")
