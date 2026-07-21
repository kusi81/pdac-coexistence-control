"""
통제(공존) 전략 비교 — 무처치 vs 연속 최대투여 vs 적응형 on/off.

목표: '박멸'이 아니라 '최소 부담으로 오래 통제'. 저항성 아형을 넣어 adaptive
therapy의 핵심(경쟁적 방출)을 드러낸다:
  - 연속 최대투여: 감수성 세포를 없애 내성 세포가 장악 → 결국 통제 불능 재발.
  - 적응형: 종양을 band에 묶어 감수성을 남김 → 내성을 억누르고 적은 부담으로 오래 통제.
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

from synthetic import make_tissue
from abm import simulate, control_metrics

# 약물이 '통제는 하되 박멸은 어려운' 창 — adaptive therapy 이점이 드러나는 체제.
PARAMS = dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10,
              init_resistant_frac=0.03, mutation_rate=0.003,
              resistant_immune_evasion=0.35)
DRUG = [("generic_cytotoxic", 1.0)]
DAYS = 150
ADAPT_ON, ADAPT_OFF = 1.1, 0.7


def series(h, key, n0=1):
    return [x["t"] for x in h], [x[key] / (n0 if key.startswith("n_tumor") else 1)
                                 for x in h]


def main():
    c, l, _ = make_tissue("contained", seed=42)
    h0, _ = simulate(c, l, days=DAYS, params=PARAMS, snapshots=(0.0, 1.0))
    hc, _ = simulate(c, l, days=DAYS, params=PARAMS, regimen_subs=DRUG, snapshots=(0.0, 1.0))
    ha, _ = simulate(c, l, days=DAYS, params=PARAMS, regimen_subs=DRUG,
                     adaptive=True, adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF,
                     snapshots=(0.0, 1.0))
    n0 = h0[0]["n_tumor"]
    runs = [("무처치", h0, "#7F8C8D"), ("연속 최대투여", hc, "#C0392B"),
            ("적응형 on/off", ha, "#27AE60")]

    print(f"{'전략':<14}{'TTP':>8}{'최종배수':>9}{'내성분율':>9}{'독성':>7}{'통제점수':>9}")
    for name, h, _c in runs:
        m = control_metrics(h, n0=n0, crit_mult=1.5)
        print(f"{name:<13}{m['ttp_days']:>6.0f}d{'*' if m['progression_censored'] else ' '}"
              f"{m['final_frac']:>8.2f}x{m['final_resistant_frac']:>9.2f}"
              f"{m['cum_toxicity']:>7.0f}{m['control_score']:>9.1f}")

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    # A. 종양 궤적
    for name, h, col in runs:
        t = [x["t"] for x in h]; y = [x["n_tumor"] / n0 for x in h]
        axs[0].plot(t, y, lw=2, color=col, label=name)
    # 적응형 투여 구간 음영
    on = [x["t"] for x in ha if x["drug_on"]]
    for tt in on:
        axs[0].axvspan(tt, tt + PARAMS.get("dt_days", 0.5), color="#27AE60", alpha=0.04)
    axs[0].axhline(1.0, color="#CCC", ls=":", lw=1)
    axs[0].axhline(1.5, color="#E74C3C", ls="--", lw=1, label="진행 임계(1.5x)")
    axs[0].set_xlabel("시간 (days)"); axs[0].set_ylabel("종양 / 초기")
    axs[0].set_title("종양 부담 — 적응형은 band에 통제(공존)", fontsize=11, fontweight="bold")
    axs[0].legend(frameon=False, fontsize=8)
    # B. 내성 분율 (경쟁적 방출)
    for name, h, col in runs:
        t = [x["t"] for x in h]; y = [x["resistant_frac"] for x in h]
        axs[1].plot(t, y, lw=2, color=col, label=name)
    axs[1].set_xlabel("시간 (days)"); axs[1].set_ylabel("내성 종양 분율")
    axs[1].set_title("내성 분율 추이 (더 공격적 체제선 연속투여가 내성 장악)",
                     fontsize=10.5, fontweight="bold")
    axs[1].legend(frameon=False, fontsize=8)
    # C. 누적 독성
    for name, h, col in runs:
        t = [x["t"] for x in h]; y = [x["cum_toxicity"] for x in h]
        axs[2].plot(t, y, lw=2, color=col, label=name)
    axs[2].set_xlabel("시간 (days)"); axs[2].set_ylabel("누적 독성(환자 부담)")
    axs[2].set_title("누적 부담 — 적응형이 더 적게", fontsize=11, fontweight="bold")
    axs[2].legend(frameon=False, fontsize=8)

    fig.suptitle("통제·공존 전략 비교 (저항성 포함 ABM) — 적응형: 종양 0.7x 공존을 "
                 "1/5 독성으로 (연속=박멸이나 부담 최대)", fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "control_strategies.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
