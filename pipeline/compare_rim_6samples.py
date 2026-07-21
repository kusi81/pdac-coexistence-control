"""
6개 PDAC Xenium 시료 rim 패널 비교 — 치료 전(naive) vs 후(CRT).

가설: CRT가 면역 배제 장벽을 풀어 CD8를 종양 rim으로 끌어온다.
      → naive는 CD8 rim 배제(면역 냉각), CRT후는 CD8 rim 농축.
각 시료를 annotate_pdac(저자 module-score) 후 run_rim_panel로 종양주변 농축 계산,
치료군별로 세포타입 rim z를 비교.

보유 시료만으로도 돌아감(부분 미리보기 → 6개 완성 시 전체).
"""
import sys, os, glob, warnings
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

import data_loader as dl
from analysis import run_rim_panel
from plots import COLORS

# 부위코드 → 치료상태 (GEO GSE274673 메타데이터)
TREATMENT = {
    "31076": ("naive", "GSM8454446"), "39928": ("naive", "GSM8454449"),
    "28429": ("naive", "GSM8454450"), "37080": ("CRT", "GSM8454447"),
    "38245": ("CRT", "GSM8454448"), "35406": ("CRT", "GSM8454451"),
}
CTS = ["myCAF", "iCAF", "apCAF", "Macrophage", "CD8_T", "Endothelial"]
SHELL = 30.0   # µm


def _region(bundle):
    for code in TREATMENT:
        if code in bundle:
            return code
    return None


def main():
    bundles = sorted(glob.glob(os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "data", "xenium", "output-*")))
    bundles = [b for b in bundles if os.path.isdir(b)]
    data = {}   # region -> {treat, gsm, z per ct}
    for b in bundles:
        code = _region(b)
        if code is None:
            continue
        treat, gsm = TREATMENT[code]
        ad = dl.load_xenium_bundle(b)
        labels, _ = dl.annotate_pdac(ad)
        coords = np.asarray(ad.obsm["spatial"], float)
        rows = run_rim_panel(coords, labels, tumor="Tumor", shell_um=SHELL, n_perms=500)
        zmap = {r["cell_type"]: r["z"] for r in rows}
        data[code] = dict(treat=treat, gsm=gsm, z=zmap, n=len(labels))
        print(f"[{gsm} {code} {treat:5}] n={len(labels):6}  " +
              "  ".join(f"{ct}:{zmap.get(ct, float('nan')):+.0f}" for ct in CTS))

    naive = [d for d in data.values() if d["treat"] == "naive"]
    crt = [d for d in data.values() if d["treat"] == "CRT"]
    print(f"\n보유: naive {len(naive)}개, CRT {len(crt)}개")
    print(f"{'세포타입':<12}{'naive 평균 rim z':>16}{'CRT 평균 rim z':>16}")
    for ct in CTS:
        nz = [d["z"][ct] for d in naive if ct in d["z"]]
        cz = [d["z"][ct] for d in crt if ct in d["z"]]
        nm = np.mean(nz) if nz else float("nan")
        cm = np.mean(cz) if cz else float("nan")
        print(f"  {ct:<12}{nm:>14.1f}  {cm:>14.1f}")

    _figure(data, naive, crt)
    return data


def _figure(data, naive, crt):
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    y = np.arange(len(CTS))
    h = 0.36
    for grp, off, hatch, alpha, lab in [(naive, +h/2, None, 0.55, "치료 전(naive)"),
                                        (crt, -h/2, None, 1.0, "CRT 치료후")]:
        means = [np.mean([d["z"][ct] for d in grp if ct in d["z"]] or [np.nan])
                 for ct in CTS]
        ax.barh(y + off, means, height=h, alpha=alpha,
                color=[COLORS.get(ct, "#999") for ct in CTS],
                edgecolor="black", linewidth=0.8, label=lab)
        # 개별 시료 점
        for d in grp:
            xs = [d["z"].get(ct, np.nan) for ct in CTS]
            ax.scatter(xs, y + off, s=16, color="black", zorder=3, alpha=0.6)
    ax.set_yticks(y); ax.set_yticklabels(CTS)
    ax.axvline(0, color="black", lw=1)
    ax.axvline(2, color="#888", ls="--", lw=1); ax.axvline(-2, color="#888", ls="--", lw=1)
    ax.set_xlabel("종양 rim(30µm) 농축 z  (>2 인접, <−2 배제) · 점=개별 시료")
    ax.set_title("치료 전 vs 후 — 종양 주변 세포 농축\n"
                 "(연한=naive, 진한=CRT후; CRT가 CD8를 rim으로 끌어오는가)",
                 fontsize=12, fontweight="bold")
    # 커스텀 범례(투명도)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="#888", alpha=0.55, label="치료 전(naive)"),
                       Patch(facecolor="#888", alpha=1.0, label="CRT 치료후")],
              frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "assets", "rim_naive_vs_crt.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
