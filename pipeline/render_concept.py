"""개념 다이어그램: 섬유화 공통 축 — 성상세포 → TGF-β → 섬유화 → 이중 운명."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa  한글 폰트
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(13, 8.4))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")


def box(x, y, w, h, text, fc, ec, fs=10, tc="#111", bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=2",
                                fc=fc, ec=ec, lw=1.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal", wrap=True)


def arrow(x1, y1, x2, y2, color="#555", lw=2.2, style="-|>", label=None, lc="#333"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=18, lw=lw, color=color,
                                 shrinkA=2, shrinkB=2))
    if label:
        ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 2.4, label, ha="center",
                fontsize=9.5, color=lc, fontweight="bold")


ax.text(50, 96, "섬유화 공통 축 — 성상세포에서 종양까지",
        ha="center", fontsize=15, fontweight="bold")
ax.text(50, 91.5, "간경화도 종양의 폭증도 같은 세포·같은 신호에서 출발한다",
        ha="center", fontsize=10.5, color="#555")

# --- 상단 흐름: 손상 → 성상세포 → (TGF-β) → 활성 근섬유아세포 ---
box(3, 74, 17, 11, "만성 손상·염증\n(간염·췌장염·종양)", "#F5EEE6", "#B08968", 10)
box(28, 74, 18, 11, "성상세포 (휴지)\nHSC · PSC\n비타민A 저장", "#EAF2F8", "#5499C7", 10)
box(56, 74, 19, 11, "활성 근섬유아세포\n= myCAF\nα-SMA+ 콜라겐 분비", "#FDEBD0", "#E67E22", 10, bold=True)
arrow(20, 79.5, 28, 79.5)
arrow(46, 79.5, 56, 79.5, color="#C0392B", label="TGF-β", lc="#C0392B")

# --- 이중 운명 분기 ---
ax.text(65.5, 69.5, "이중 운명 (context-dependent)", ha="center", fontsize=9.5,
        color="#7F8C8D", style="italic")
arrow(63, 74, 40, 60, color="#27AE60")     # to restraining (left-down)
arrow(68, 74, 88, 60, color="#C0392B")     # to promoting (right-down)

# 왼쪽: 억제 (처음 의도)
box(19, 45, 30, 14,
    "① 억제 (처음 의도)\n종양 포위·기계적 압박\n→ 증식 제지 / 면역 배제",
    "#E9F7EF", "#27AE60", 10, bold=True)
# 오른쪽: 촉진 (포섭 후)
box(65, 45, 31, 14,
    "② 촉진 (포섭 후)\nECM 리모델링·성장인자\n→ 증식·전이 촉진",
    "#FDEDEC", "#C0392B", 10, bold=True)

# 포섭 화살표 (억제→촉진 전환)
arrow(49, 52, 65, 52, color="#8E44AD", lw=2.6, label="암에게 '포섭'", lc="#8E44AD")
ax.text(57, 48.4, "LOX–FAK–β-catenin\n촉진 상태로 잠김(lock)", ha="center",
        fontsize=8.2, color="#8E44AD")
# 가역성 (되돌리기)
arrow(65, 56, 49, 56, color="#16A085", lw=1.8, style="-|>",
      label="항섬유화로 되돌리기", lc="#16A085")

# --- 장기 맥락 (하단) ---
box(8, 22, 38, 15,
    "췌장 (PDAC)\n장벽 우세 — desmoplasia가 면역배제\nCAF 제거가 오히려 악화(역설)\nABM: caf_protumor = 0",
    "#EAF2F8", "#2E86AB", 9.5)
box(54, 22, 38, 15,
    "간 (HCC) / 간경화\n전암 토양 우세 — 섬유화가 증식 촉진\n경화된 간에서 HCC 80–90% 발생\nABM: caf_protumor > 0",
    "#FDEDEC", "#C0392B", 9.5)
arrow(34, 45, 27, 37, color="#2E86AB", lw=1.8)
arrow(80, 45, 73, 37, color="#C0392B", lw=1.8)

# --- 치료 논리 (맨 아래 배너) ---
box(15, 4, 70, 11,
    "치료 논리: 가소성(plasticity) — 촉진→억제 표현형으로 되돌리기\n"
    "커큐민·쑥(artesunate)·산삼·단삼·황기 = TGF-β 축 차단 → myCAF/HSC 재프로그래밍",
    "#F4ECF7", "#8E44AD", 9.8, bold=True)

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "assets", "concept_fibrosis_axis.png")
fig.savefig(out, dpi=115, bbox_inches="tight")
print("wrote", out)
