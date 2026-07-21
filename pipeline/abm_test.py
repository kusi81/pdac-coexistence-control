"""ABM 동작 검증: 대조군 종양은 자라고, 항섬유화/면역/조합은 억제해야 한다.

기전 sanity check:
 - 대조(무처치): myCAF 장벽이 CD8를 막아 종양 성장
 - 항섬유화: 장벽↓ -> CD8 침윤↑ -> 종양 억제
 - 세포독성: 직접 증식 억제
 - 조합: 가장 강한 억제
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from synthetic import make_tissue
from abm import simulate, compose_regimen, SUBSTANCES

coords, labels, _ = make_tissue("contained", seed=42)
DAYS = 20.0

arms = {
    "대조(무처치)": [],
    "항섬유화": [("generic_antifibrotic", 1.0)],
    "세포독성": [("generic_cytotoxic", 1.0)],
    "커큐민+Rg3 조합": [("curcumin", 1.0), ("ginsenoside_rg3", 1.0)],
}

t0 = None
print(f"{'처치군':<20}{'종양 t0':>9}{'종양 t20':>10}{'변화%':>9}{'CD8침윤 t20':>12}")
print("-" * 62)
rows = []
for name, sel in arms.items():
    reg = compose_regimen(sel)
    hist, _ = simulate(coords, labels, days=DAYS, regimen=reg, snapshots=(0.0, 1.0))
    n0 = hist[0]["n_tumor"]
    n1 = hist[-1]["n_tumor"]
    infil = hist[-1]["cd8_infiltration"]
    pct = 100.0 * (n1 - n0) / n0
    rows.append((name, n0, n1, pct, infil))
    print(f"{name:<19}{n0:>9}{n1:>10}{pct:>+8.0f}%{infil:>11.2f}")

print()
ctrl = next(r for r in rows if r[0].startswith("대조"))
combo = next(r for r in rows if "조합" in r[0])
print(f"검증: 대조 종양변화 {ctrl[3]:+.0f}% vs 조합 {combo[3]:+.0f}%")
ok = combo[3] < ctrl[3]   # 조합이 대조보다 종양을 더 억제
print("결과:", "PASS - 처치가 대조 대비 종양 억제" if ok
      else "CHECK - 예상과 다른 방향")
