"""완성 F3 — SCOTIA 저자주석으로 종양-rim 세포조성 (Untreated vs CRT).

Xenium 주석 실패로 깨졌던 F3를, 저자주석 CosMx(양성대조 15/16 통과)로 재구성.
각 시료: run_rim_panel(tumor=Malignant, 30µm) → 세포타입별 rim z.
치료군(Untreated/CRT)별 평균 + 개별 시료 점. 그림 assets/rim_scotia.png.
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from analysis import run_rim_panel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5 = os.path.join(ROOT, "data", "scotia", "raw_meta_data_final.h5ad")
PX_UM = 0.12028
SHELL = 30.0
# 표시할 세포타입(존재하는 것만)
CTS = ["myCAF", "iCAF", "apCAF", "CD8+ T", "CD4+ T", "Treg", "Macrophage",
       "Plasma", "Endothelial", "Pericyte"]
CMAP = {"myCAF": "#C0392B", "iCAF": "#E67E22", "apCAF": "#8E44AD",
        "CD8+ T": "#2980B9", "CD4+ T": "#5DADE2", "Treg": "#85C1E9",
        "Macrophage": "#16A085", "Plasma": "#7F8C8D", "Endothelial": "#27AE60",
        "Pericyte": "#95A5A6"}


def build_labels(obs):
    sub = obs["annotation_subtypes"].astype(str).values
    maj = obs["annotation_majortypes"].astype(str).values
    lab = maj.copy()                       # 대분류 기본
    caf = maj == "CAF"
    lab[caf] = sub[caf]                    # CAF는 서브타입으로
    lab[maj == "Malignant"] = "Malignant"
    return lab


def main():
    a = ad.read_h5ad(H5, backed="r")
    obs = a.obs
    lab = build_labels(obs)
    X = obs["CenterX_global_px"].astype(float).values * PX_UM
    Y = obs["CenterY_global_px"].astype(float).values * PX_UM
    samp = obs["sample_id"].astype(str).values
    treat = obs["treatment_status"].astype(str).values

    data = {}
    for s in pd.unique(samp):
        m = samp == s
        tr = treat[m][0]
        if tr not in ("Untreated", "CRT"):     # CRTL(n=1) 제외
            continue
        coords = np.column_stack([X[m], Y[m]])
        labels = lab[m]
        rows = run_rim_panel(coords, labels, tumor="Malignant",
                             shell_um=SHELL, n_perms=300)
        z = {r["cell_type"]: r["z"] for r in rows}
        data[s] = dict(treat=tr, z=z, n=int(m.sum()))
        print(f"[{s:8} {tr:9}] n={m.sum():6}  " +
              "  ".join(f"{c}:{z.get(c, float('nan')):+.0f}" for c in CTS
                        if c in z))

    naive = [d for d in data.values() if d["treat"] == "Untreated"]
    crt = [d for d in data.values() if d["treat"] == "CRT"]
    present = [c for c in CTS if any(c in d["z"] for d in data.values())]
    print(f"\n보유: Untreated {len(naive)}, CRT {len(crt)}")
    print(f"{'세포타입':<12}{'Untreated 평균 z':>16}{'CRT 평균 z':>14}")
    summ = {}
    for c in present:
        nz = [d["z"][c] for d in naive if c in d["z"]]
        cz = [d["z"][c] for d in crt if c in d["z"]]
        nm = float(np.mean(nz)) if nz else float("nan")
        cm = float(np.mean(cz)) if cz else float("nan")
        summ[c] = (nm, cm)
        print(f"  {c:<12}{nm:>14.1f}  {cm:>14.1f}")

    _figure(present, data, naive, crt)
    pd.DataFrame([{"sample": s, **{"treat": d["treat"]},
                   **{c: d["z"].get(c) for c in present}}
                  for s, d in data.items()]).to_csv(
        os.path.join(ROOT, "data", "scotia", "scotia_rim.csv"), index=False)
    print("wrote scotia_rim.csv")


def _figure(present, data, naive, crt):
    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(present))
    h = 0.38
    for grp, off, alpha, lab in [(naive, +h/2, 0.55, "치료 전(Untreated)"),
                                 (crt, -h/2, 1.0, "CRT 치료후")]:
        means = [np.mean([d["z"][c] for d in grp if c in d["z"]] or [np.nan])
                 for c in present]
        ax.barh(y + off, means, height=h, alpha=alpha,
                color=[CMAP.get(c, "#999") for c in present],
                edgecolor="black", linewidth=0.8)
        for d in grp:
            xs = [d["z"].get(c, np.nan) for c in present]
            ax.scatter(xs, y + off, s=14, color="black", zorder=3, alpha=0.55)
    ax.set_yticks(y); ax.set_yticklabels(present)
    ax.invert_yaxis()
    ax.axvline(0, color="black", lw=1)
    ax.axvline(2, color="#888", ls="--", lw=1)
    ax.axvline(-2, color="#888", ls="--", lw=1)
    ax.set_xlabel("종양 rim(30µm) 농축 z  (>2 인접, <-2 배제) · 점=개별 시료")
    ax.set_title("SCOTIA 저자주석 — 종양 주변 세포 농축 (Untreated vs CRT)\n"
                 "양성대조 15/16 통과 · 대형패널 CosMx", fontsize=12.5,
                 fontweight="bold")
    ax.legend(handles=[Patch(facecolor="#888", alpha=0.55, label="치료 전(Untreated)"),
                       Patch(facecolor="#888", alpha=1.0, label="CRT 치료후")],
              frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    out = os.path.join(ROOT, "assets", "rim_scotia.png")
    fig.savefig(out, dpi=115, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
