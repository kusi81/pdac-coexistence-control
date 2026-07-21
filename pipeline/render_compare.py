"""contained vs diffuse 조직 지도 + 장벽 점수를 한 장의 PNG로 렌더."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa: F401  (import 부작용으로 한글 폰트 등록)
import matplotlib.pyplot as plt

from synthetic import make_tissue
from analysis import run_core_metrics
from plots import COLORS

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "assets", "compare_modes.png")

fig, axes = plt.subplots(2, 2, figsize=(12, 11))
plt.rcParams.update({"font.size": 9})

for col, mode in enumerate(["contained", "diffuse"]):
    coords, labels, _ = make_tissue(mode=mode, seed=42)
    res = run_core_metrics(coords, labels, n_perms=1000, seed=0)

    # --- top row: tissue map ---
    ax = axes[0, col]
    order = ["Macrophage", "iCAF", "myCAF", "CD8_T", "Tumor"]
    for ct in order:
        m = labels == ct
        ax.scatter(coords[m, 0], coords[m, 1], s=5, c=COLORS.get(ct, "#999"),
                   alpha=0.75, linewidths=0, label=f"{ct} ({int(m.sum())})")
    ax.set_aspect("equal")
    ax.set_title(f"{mode.upper()} tissue map", fontsize=12, fontweight="bold")
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    if col == 1:
        ax.legend(frameon=False, fontsize=7.5, markerscale=2.4,
                  loc="center left", bbox_to_anchor=(1.01, 0.5))

    # --- bottom row: barrier z ---
    ax = axes[1, col]
    names = list(res["barrier"].keys())
    zs = [res["barrier"][n]["z_score"] for n in names]
    cols = ["#C0392B" if z > 2 else "#7F8C8D" for z in zs]
    bars = ax.barh(names, zs, color=cols, height=0.5)
    ax.axvline(2, color="black", ls="--", lw=1)
    ax.set_xlim(-20, 25)
    ax.set_xlabel("barrier z  (>2 = interposed beyond chance)")
    ax.set_title(f"{mode.upper()} barrier score", fontsize=12, fontweight="bold")
    for b, z in zip(bars, zs):
        ax.text(z + (0.6 if z >= 0 else -0.6), b.get_y() + b.get_height() / 2,
                f"{z:+.1f}", va="center", ha="left" if z >= 0 else "right",
                fontsize=10, fontweight="bold")

fig.suptitle("동일 세포 수 · 다른 구조 — 구성비로는 구분 불가, 위치(장벽 점수)로만 구분됨",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig(OUT, dpi=110, bbox_inches="tight")
print("wrote", OUT)
