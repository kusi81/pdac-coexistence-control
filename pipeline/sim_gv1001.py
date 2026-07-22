"""GV1001 순차 병용 — '장벽을 연 뒤 면역 부스트' 가설 시뮬레이션.

개념: 조밀 myCAF는 CD8를 배제하므로 면역 프라이밍(GV1001) 단독은 무력하다. 항섬유화로
장벽을 연 창(window)에서 GV1001을 투여하면 T세포가 종양에 도달해 효과를 낸다.
→ 5개 팔 비교로 '무의미가 아니라 타이밍의 문제'임을 보인다.

  1 무처치
  2 GV1001 단독(연속)       — 장벽 배제로 약할 것으로 예상
  3 항섬유화 단독(연속)      — 장벽은 열되 면역 부스트 없음
  4 동시병용(항섬유+GV1001)  — 장벽 열림 + 면역 부스트
  5 순차(개방→열린창 GV1001→리셋 사이클)

면역배제가 유효한 체제(cd8_barrier_alpha 높음)에서, containment 기본 ON.
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt
matplotlib.rcParams["axes.unicode_minus"] = False
from matplotlib.lines import Line2D
from synthetic import make_tissue
from abm import simulate, control_metrics

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 강한 면역배제 체제 — GV1001 단독은 장벽에 막혀 실패, 개방(항섬유화) 필요
PARAMS = dict(k_prolif=0.15, cd8_recruit=8, k_kill=0.6, cd8_barrier_alpha=3.0)
DAYS = 120
SEEDS = [1, 2, 3]
ANTIFIB = [("danshen", 1.0), ("astragaloside", 1.0)]      # 순수 항섬유화(장벽 개방)
GV = [("gv1001", 1.0)]
# 순차 스케줄: 개방 → 열린 창에서 GV1001 → 휴지(장벽 재형성) 반복
SEQ = [
    dict(label="open barrier", days=6, substances=ANTIFIB),
    dict(label="immune strike (open window)", days=6, substances=ANTIFIB + GV),
    dict(label="rest (barrier re-forms)", days=6, substances=[]),
]
ARMS = [
    ("Untreated", "none", dict(), "#7F8C8D"),
    ("GV1001 alone", "subs", dict(regimen_subs=GV), "#2980B9"),
    ("Anti-fibrotic alone", "subs", dict(regimen_subs=ANTIFIB), "#E67E22"),
    ("Concurrent (anti-fib + GV1001)", "subs", dict(regimen_subs=ANTIFIB + GV), "#8E44AD"),
    ("Sequential (open → GV1001 → reset)", "sched", dict(schedule=SEQ), "#27AE60"),
]


def run(kind, kw, seed):
    p = dict(PARAMS); p["seed"] = seed
    h, _ = simulate(c, l, days=DAYS, params=p, snapshots=(0.0, 1.0), **kw)
    return h, control_metrics(h, n0=n0)


c, l, _ = make_tissue("contained", seed=42)
n0 = (l == "Tumor").sum()


def main():
    print(f"{'arm':<36}{'final':>7}{'resist':>8}{'toxic':>7}{'ctrl':>7}")
    results = {}
    for name, kind, kw, col in ARMS:
        hs, ms = [], []
        for s in SEEDS:
            h, m = run(kind, kw, s)
            hs.append(h); ms.append(m)
        # seed 평균 궤적(공통 시점 가정) + 지표 평균
        agg = {k: float(np.mean([m[k] for m in ms]))
               for k in ("final_frac", "final_resistant_frac", "cum_toxicity",
                         "control_score")}
        results[name] = (hs[0], agg, col)
        print(f"{name:<36}{agg['final_frac']:>6.2f}x{agg['final_resistant_frac']:>8.2f}"
              f"{agg['cum_toxicity']:>7.0f}{agg['control_score']:>7.1f}")

    # 그림: 궤적 + 요약 막대
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 5.4))
    for name, (h, agg, col) in results.items():
        t = [x["t"] for x in h]; y = [x["n_tumor"] / n0 for x in h]
        axA.plot(t, y, lw=2, color=col, label=name)
    axA.axhline(1.5, color="#E74C3C", ls=":", lw=1)
    axA.text(2, 1.53, "progression threshold 1.5x", fontsize=8, color="#E74C3C")
    axA.set_xlabel("Time (days)"); axA.set_ylabel("Tumor / initial")
    axA.set_title("a  Tumor trajectory — GV1001 needs an open barrier to act",
                  fontsize=10.5, fontweight="bold")
    axA.legend(frameon=False, fontsize=7.6)

    names = list(results)
    finals = [results[n][1]["final_frac"] for n in names]
    cols = [results[n][2] for n in names]
    y = np.arange(len(names))
    axB.barh(y, finals, color=cols, edgecolor="black", linewidth=0.7)
    for i, f in enumerate(finals):
        axB.text(f + 0.02, i, f"{f:.2f}x", va="center", fontsize=8.5)
    axB.set_yticks(y); axB.set_yticklabels(names, fontsize=8.5); axB.invert_yaxis()
    axB.set_xlabel("Final tumor burden (fold vs initial)  →  lower is better")
    axB.set_title("b  GV1001 alone ≈ futile (barrier excludes T cells);\n"
                  "opening the barrier first unlocks it", fontsize=10.5, fontweight="bold")
    fig.suptitle("GV1001 (telomerase vaccine) sequential combination — timing, not futility "
                 "(in-silico)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(ROOT, "assets", "gv1001_sequential.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
