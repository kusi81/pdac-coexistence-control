"""Fig 3 재작도 (글리프 수정) — scotia_rim.csv에서만 읽음, 재계산 없음."""
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
import fonts  # noqa (sets axes.unicode_minus=False)
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
df = pd.read_csv(os.path.join(ROOT, "data", "scotia", "scotia_rim.csv"))
CTS = ["myCAF", "iCAF", "apCAF", "CD8+ T", "CD4+ T", "Treg", "Macrophage",
       "Plasma", "Endothelial", "Pericyte"]
CTS = [c for c in CTS if c in df.columns]
CMAP = {"myCAF": "#C0392B", "iCAF": "#E67E22", "apCAF": "#8E44AD",
        "CD8+ T": "#2980B9", "CD4+ T": "#5DADE2", "Treg": "#85C1E9",
        "Macrophage": "#16A085", "Plasma": "#7F8C8D", "Endothelial": "#27AE60",
        "Pericyte": "#95A5A6"}
naive = df[df.treat == "Untreated"]
crt = df[df.treat == "CRT"]

fig, ax = plt.subplots(figsize=(10, 6))
y = np.arange(len(CTS)); h = 0.38
for grp, off, alpha in [(naive, +h/2, 0.55), (crt, -h/2, 1.0)]:
    means = [grp[c].mean() for c in CTS]
    ax.barh(y + off, means, height=h, alpha=alpha,
            color=[CMAP.get(c, "#999") for c in CTS],
            edgecolor="black", linewidth=0.8)
    for _, r in grp.iterrows():
        ax.scatter([r[c] for c in CTS], y + off, s=14, color="black",
                   zorder=3, alpha=0.55)
ax.set_yticks(y); ax.set_yticklabels(CTS); ax.invert_yaxis()
ax.axvline(0, color="black", lw=1)
ax.axvline(2, color="#888", ls="--", lw=1); ax.axvline(-2, color="#888", ls="--", lw=1)
ax.set_xlabel("종양 rim(30µm) 농축 z  (>2 인접, <-2 배제) · 점=개별 시료")
ax.set_title("SCOTIA 저자주석 — 종양 주변 세포 농축 (Untreated vs CRT)\n"
             "양성대조 15/16 통과 · 대형패널 CosMx", fontsize=12.5, fontweight="bold")
ax.legend(handles=[Patch(facecolor="#888", alpha=0.55, label="치료 전(Untreated)"),
                   Patch(facecolor="#888", alpha=1.0, label="CRT 치료후")],
          frameon=False, fontsize=9, loc="lower right")
fig.tight_layout()
out = os.path.join(ROOT, "assets", "rim_scotia.png")
fig.savefig(out, dpi=115, bbox_inches="tight")
print("wrote", out)
