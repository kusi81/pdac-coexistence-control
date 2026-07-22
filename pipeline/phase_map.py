"""위상도 (L1) — '기질=자원' 가설이 성립하는 조건 지도.

리뷰어 방어: myCAF containment를 모델에 구현(confinement+압력+약물차단) 후, 이 가설이
'항상' 성립한다고 주장하지 않고 **성립 조건을 지도화**한다.

두 힘의 트레이드오프:
  x축 = confinement 강도(caf_pressure) — 기질의 이득(종양 가둠, 경계효과)
  y축 = 면역배제 강도(cd8_barrier_alpha) — 기질의 비용(CD8 차단, 전면효과)
각 (x,y)에서 myCAF 수준(k_caf_activate)을 스윕 → 최적 myCAF 수준(종양 최소)을 찾는다.
  최적 myCAF > 0  → '기질을 유지·조절'이 최선 = 자원 가설 성립
  최적 myCAF = 0  → '기질 제거'가 최선 = 가설 불성립(면역억제 우세)

고정 치료 하 escape 체제. 텐서는 data/phase_map.npz에 저장(재작도용).
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np
from synthetic import make_tissue
from abm import simulate, control_metrics

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = dict(k_prolif=0.17, cd8_recruit=6, k_kill=0.40,
            caf_confine=0.9, caf_confine_ref=2.5)
DRUG = [("generic_cytotoxic", 0.30)]
DAYS = 90
SEEDS = [1, 2, 3]
PRESSURE = [0.0, 0.9, 1.8, 2.8]        # x: confinement 강도
ALPHA = [0.2, 0.6, 1.2, 2.0]           # y: 면역배제 강도
MYCAF = [0.0, 0.05, 0.11, 0.20, 0.30]  # myCAF 수준 스윕


def main():
    c, l, _ = make_tissue("contained", seed=42)
    n0 = (l == "Tumor").sum()
    # tumor[xi, yi, mi] = 평균 최종 종양/초기
    T = np.zeros((len(PRESSURE), len(ALPHA), len(MYCAF)))
    for xi, pres in enumerate(PRESSURE):
        for yi, alpha in enumerate(ALPHA):
            for mi, act in enumerate(MYCAF):
                vals = []
                for s in SEEDS:
                    p = dict(BASE); p.update(caf_pressure=pres,
                                             cd8_barrier_alpha=alpha,
                                             k_caf_activate=act, seed=s)
                    h, _ = simulate(c, l, days=DAYS, params=p,
                                    regimen_subs=DRUG, snapshots=(0.0, 1.0))
                    vals.append(control_metrics(h, n0=n0)["final_frac"])
                T[xi, yi, mi] = float(np.mean(vals))
            opt = MYCAF[int(np.argmin(T[xi, yi]))]
            print(f"pressure={pres:.1f} alpha={alpha:.1f} → "
                  + " ".join(f"{a}:{v:.2f}" for a, v in zip(MYCAF, T[xi, yi]))
                  + f"  최적myCAF={opt}")
    np.savez(os.path.join(ROOT, "data", "phase_map.npz"),
             T=T, pressure=PRESSURE, alpha=ALPHA, mycaf=MYCAF)
    print("\nsaved data/phase_map.npz — plot with pipeline/plot_phase_map.py")


if __name__ == "__main__":
    main()
