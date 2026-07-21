"""
에를로티닙 · 젬시타빈 vs 커큐민 3D 구조 비교 — '실험 공결정 vs 도킹 예측'.

  1M17 (EGFR + 에를로티닙, AQ4)   : 실험 — 약물이 ATP 주머니에 실측 (초록)
  1P62 (dCK + 젬시타빈, GEO)      : 실험 — 약물이 활성화효소에 실측 (초록)
  1BG1 (STAT3, 커큐민 도킹)       : 예측 — 약물 없음, 예측 잔기만 (빨강)

핵심: 같은 '3D 표적+결합부위'라도 근거 수준이 결론을 가른다.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
import numpy as np
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt

PDB_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "data", "pdb")


def parse_pdb(path, lig_resn=None):
    """반환: ca(dict resSeq->xyz, 주 사슬), lig(list xyz, 지정 리간드 원자)."""
    by_chain = {}     # chain -> {resSeq: xyz(CA)}
    lig = []
    for line in open(path):
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        ch = line[21]; resn = line[17:20].strip(); atom = line[12:16].strip()
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        if rec == "ATOM" and atom == "CA":
            try:
                by_chain.setdefault(ch, {})[int(line[22:26])] = xyz
            except ValueError:
                pass
        elif rec == "HETATM" and lig_resn and resn == lig_resn and atom[0] != "H":
            lig.append(xyz)
    ca = max(by_chain.values(), key=len) if by_chain else {}   # 최대 사슬
    return ca, lig


PANELS = [
    ("1M17.pdb", "AQ4", "#7FB3D5", "에를로티닙",
     "[실험] EGFR + 에를로티닙 (1M17)\n약물이 ATP 주머니에 실측", "#1E8449", None),
    ("1P62.pdb", "GEO", "#82C4B0", "젬시타빈",
     "[실험] dCK + 젬시타빈 (1P62)\n약물이 활성화효소에 실측", "#1E8449", None),
    ("1BG1.pdb", None, "#F5B7B1", None,
     "[도킹 예측] STAT3 + 커큐민 (1BG1)\n약물 없음, 예측 잔기만", "#B9770E",
     (474, 325, 324, 252, 258, 250)),
]


def _clean(ax):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass


def main():
    fig = plt.figure(figsize=(16, 5.6))
    for k, (pdb, lig_resn, pcol, drug, title, tcol, hi_res) in enumerate(PANELS):
        ca, lig = parse_pdb(os.path.join(PDB_DIR, pdb), lig_resn)
        ax = fig.add_subplot(1, 3, k + 1, projection="3d")
        seq = sorted(ca); P = np.array([ca[i] for i in seq])
        ax.plot(P[:, 0], P[:, 1], P[:, 2], color=pcol, lw=1.0, alpha=0.6)
        if lig:
            L = np.array(lig)
            ax.scatter(L[:, 0], L[:, 1], L[:, 2], color="#27AE60", s=55,
                       edgecolor="black", linewidth=0.3,
                       label=f"{drug} (실측, {len(L)}원자)")
            print(f"{pdb}: CA {len(ca)}, {drug}({lig_resn}) {len(lig)}원자")
        if hi_res:
            hs = [r for r in hi_res if r in ca]
            if hs:
                H = np.array([ca[r] for r in hs])
                ax.scatter(H[:, 0], H[:, 1], H[:, 2], color="#C0392B", s=75,
                           edgecolor="black", linewidth=0.4,
                           label=f"도킹 예측 잔기 ({len(hs)})")
            print(f"{pdb}: CA {len(ca)}, 예측잔기 {hs}")
        ax.set_title(title, fontsize=10.5, fontweight="bold", color=tcol)
        ax.legend(fontsize=8, loc="upper left")
        _clean(ax)

    fig.suptitle("에를로티닙·젬시타빈(실험) vs 커큐민(도킹) — 초록 약물이 보이면 자세 신뢰, "
                 "빨강 예측 잔기만이면 가설", fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.normpath(os.path.join(PDB_DIR, "..", "..", "assets",
                                        "drug_structures_3d.png"))
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
