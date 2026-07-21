"""민감도 그림 재작성 (정직 프레이밍) — CSV에서만 읽음, 재시뮬 없음.

정정: 적응형 저독성 이점은 resistance_cost 전 범위서 강건(cost 비의존).
resistance_cost의 실제 역할 = 내성 억제(경쟁적 방출): cost 0→내성 0.20, cost≥0.1→<0.05.
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt
matplotlib.rcParams["axes.unicode_minus"] = False

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_CS = 5.57
sw = pd.read_csv(os.path.join(ROOT, "data", "sens_rcost.csv"))
tor = pd.read_csv(os.path.join(ROOT, "data", "sens_tornado.csv")).sort_values("span")

fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.6))

# ── a: resistance_cost — 이점 강건 + 내성 억제가 진짜 역할 ──
ax = axs[0]
ax.plot(sw.resistance_cost, sw.adapt_cs, "o-", color="#27AE60", lw=2,
        label="적응형 control_score")
ax.plot(sw.resistance_cost, sw.cont_cs, "s-", color="#C0392B", lw=2,
        label="연속 control_score")
ax.axvspan(0.20, 0.30, color="#3498DB", alpha=0.10)
ax.text(0.25, 5.0, "문헌범위\n0.2-0.3", ha="center", fontsize=8, color="#2471A3")
ax.set_xlabel("resistance_cost (내성 적합도 비용)")
ax.set_ylabel("control_score (통제기간/(독성+1))")
ax.set_ylim(0, 10)
ax2 = ax.twinx()
ax2.plot(sw.resistance_cost, sw.adapt_res, "^--", color="#8E44AD", lw=1.8,
         label="적응형 최종 내성분율")
ax2.set_ylabel("최종 내성 분율", color="#8E44AD")
ax2.tick_params(axis="y", labelcolor="#8E44AD")
ax2.set_ylim(0, 0.25)
l1, la = ax.get_legend_handles_labels()
l2, lb = ax2.get_legend_handles_labels()
ax.legend(l1 + l2, la + lb, frameon=False, fontsize=8, loc="center right")
ax.set_title("a  적응형 저독성 통제는 cost 전범위서 강건(연속 대비 CS 4.5~5.9 vs 1.2)\n"
             "resistance_cost의 역할 = 내성 억제: cost 0→0.20, cost≥0.1→<0.05",
             fontsize=9.8, fontweight="bold")

# ── b: 토네이도 ──
ax = axs[1]
y = np.arange(len(tor))
for i, (_, r) in enumerate(tor.iterrows()):
    ax.plot([r.cs_lo, r.cs_hi], [i, i], color="#95A5A6", lw=2, zorder=1)
    ax.scatter([r.cs_lo], [i], color="#2980B9", s=45, zorder=2,
               label="0.5×" if i == 0 else None)
    ax.scatter([r.cs_hi], [i], color="#E67E22", s=45, zorder=2,
               label="2×" if i == 0 else None)
ax.axvline(BASE_CS, color="black", ls="--", lw=1.2, label="baseline")
ax.set_yticks(y); ax.set_yticklabels(tor.param)
ax.set_xlabel("적응형 control_score")
ax.set_title("b  OAT 민감도(±50%): 결과는 종양-면역 균형에 지배\n"
             "(k_prolif·cd8_barrier_alpha·k_kill) · 내성 파라미터는 부차적",
             fontsize=9.8, fontweight="bold")
ax.legend(frameon=False, fontsize=8, loc="lower right")

fig.suptitle("Figure S3 — 민감도분석: 적응형 저독성 통제는 파라미터에 강건; "
             "결과 민감도는 종양-면역 균형에 집중(→ 우선 실측 대상)",
             fontsize=11.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(ROOT, "assets", "sensitivity.png")
fig.savefig(out, dpi=118, bbox_inches="tight")
print("wrote", out)
