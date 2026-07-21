"""ABM 시뮬레이션 결과 시각화 (Streamlit용 Figure 반환)."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa: F401  한글 폰트 등록
import matplotlib.pyplot as plt

from plots import COLORS

plt.rcParams.update({"figure.dpi": 130, "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})


def _series(hist, key):
    return np.array([h["t"] for h in hist]), np.array([h[key] for h in hist])


def fig_tumor_trajectory(hist_ctrl, hist_treated, label_treated="처치"):
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    t, y = _series(hist_ctrl, "n_tumor")
    ax.plot(t, y, "-o", ms=3, lw=1.8, color="#7F8C8D", label="대조(무처치)")
    t2, y2 = _series(hist_treated, "n_tumor")
    ax.plot(t2, y2, "-o", ms=3, lw=1.8, color="#C0392B", label=label_treated)
    ax.axhline(y[0], color="#BBB", ls=":", lw=1)
    ax.set_xlabel("시간 (days)")
    ax.set_ylabel("종양 세포 수")
    ax.set_title("종양 성장 궤적 — 대조 vs 처치", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    return fig


def fig_mechanism_trajectory(hist_ctrl, hist_treated, label_treated="처치"):
    """CD8 침윤 + myCAF 장벽 세포수 궤적 (기전 판독)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.0))

    for hist, name, col, ls in [(hist_ctrl, "대조", "#7F8C8D", "-"),
                                (hist_treated, label_treated, "#27AE60", "-")]:
        t, inf = _series(hist, "cd8_infiltration")
        ax1.plot(t, inf, ls, marker="o", ms=3, lw=1.7, color=col, label=name)
    ax1.set_xlabel("시간 (days)"); ax1.set_ylabel("CD8 침윤도 (종양 25µm 내 비율)")
    ax1.set_title("면역 침투 — 장벽이 풀리면 상승", fontsize=10)
    ax1.legend(frameon=False, fontsize=8.5)

    for hist, name, col in [(hist_ctrl, "대조", "#7F8C8D"),
                            (hist_treated, label_treated, "#C0392B")]:
        t, my = _series(hist, "n_myCAF")
        ax2.plot(t, my, "-o", ms=3, lw=1.7, color=col, label=name)
    ax2.set_xlabel("시간 (days)"); ax2.set_ylabel("myCAF 장벽 세포 수")
    ax2.set_title("장벽 질량 — 항섬유화 처치 시 하락", fontsize=10)
    ax2.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    return fig


def fig_snapshot_single(coords, labels, t, size=4.2):
    """단일 시점 조직 스냅샷 (정사각) — Streamlit 컬럼에 넣어 반응형으로 배치."""
    fig, ax = plt.subplots(figsize=(size, size))
    order = ["Macrophage", "iCAF", "myCAF", "CD8_T", "Tumor"]
    for ct in order:
        m = labels == ct
        if not m.any():
            continue
        ax.scatter(coords[m, 0], coords[m, 1], s=4, c=COLORS.get(ct, "#999"),
                   alpha=0.75, linewidths=0)
    n_tum = int((labels == "Tumor").sum())
    ax.set_aspect("equal")
    ax.set_title(f"t = {t:.0f}d  (종양 {n_tum})", fontsize=10, fontweight="bold")
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    fig.tight_layout()
    return fig


def fig_snapshots(snaps, title_prefix=""):
    """조직 스냅샷 여러 시점을 가로로 나열."""
    times = sorted(snaps.keys())
    n = len(times)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 4.6))
    if n == 1:
        axes = [axes]
    order = ["Macrophage", "iCAF", "myCAF", "CD8_T", "Tumor"]
    for ax, t in zip(axes, times):
        coords, labels = snaps[t]
        for ct in order:
            m = labels == ct
            if not m.any():
                continue
            ax.scatter(coords[m, 0], coords[m, 1], s=4, c=COLORS.get(ct, "#999"),
                       alpha=0.75, linewidths=0,
                       label=f"{ct} ({int(m.sum())})" if ax is axes[-1] else None)
        n_tum = int((labels == "Tumor").sum())
        ax.set_aspect("equal")
        ax.set_title(f"{title_prefix}t = {t:.0f}d  (종양 {n_tum})", fontsize=10)
        ax.set_xlabel("x (µm)")
    axes[0].set_ylabel("y (µm)")
    axes[-1].legend(frameon=False, fontsize=7, markerscale=2.2,
                    loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    return fig


def fig_arms_bar(arm_results):
    """여러 처치군의 최종 종양 변화% 막대 비교.
    arm_results: [(name, pct_change, infiltration), ...]"""
    if not arm_results:
        return None
    names = [a[0] for a in arm_results]
    pcts = [a[1] for a in arm_results]
    fig, ax = plt.subplots(figsize=(7.2, 0.6 * len(names) + 1.6))
    cols = ["#7F8C8D" if n.startswith("대조") else
            ("#C0392B" if p >= 0 else "#27AE60") for n, p in zip(names, pcts)]
    ax.barh(names, pcts, color=cols, height=0.6)
    ax.axvline(0, color="black", lw=1)
    for i, p in enumerate(pcts):
        ax.text(p + (1 if p >= 0 else -1), i, f"{p:+.0f}%", va="center",
                ha="left" if p >= 0 else "right", fontsize=9, fontweight="bold")
    ax.set_xlabel("종양 세포 수 변화 % (음수 = 억제)")
    ax.set_title("처치군별 종양 성장 억제 비교", fontsize=11, fontweight="bold")
    fig.tight_layout()
    return fig
