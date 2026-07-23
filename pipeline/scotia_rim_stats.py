"""CRT rim 분석 — patient-level 통계 재작성 (리뷰어 #6).

기존 Fig 3은 세포를 pooled한 'mean rim z' 위주였다. 여기서는 **환자(sample)를 분석
단위**로 하여, 세포타입별 Untreated(n=9) vs CRT(n=6) 그룹 비교를 정식 통계로 제시:
  - 각 환자의 rim z (양성대조 통과 CosMx 저자주석; run_rim_panel)
  - 그룹 평균 ± 95% CI (환자 단위 부트스트랩)
  - Mann-Whitney U (소표본 비모수), 효과크기 Cliff's delta
  - Benjamini-Hochberg 다중검정 보정(세포타입 수)
인과적 표현 없이 '치료상태와 rim 조성의 횡단면 연관'으로 보고.

결과: data/scotia/scotia_rim_stats.csv, assets/rim_scotia_stats.png (English).
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
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu
from analysis import run_rim_panel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5 = os.path.join(ROOT, "data", "scotia", "raw_meta_data_final.h5ad")
PX_UM = 0.12028
SHELL = 30.0
CTS = ["myCAF", "iCAF", "apCAF", "CD8+ T", "CD4+ T", "Treg", "Macrophage",
       "Plasma", "Endothelial", "Pericyte"]
CMAP = {"myCAF": "#C0392B", "iCAF": "#E67E22", "apCAF": "#8E44AD",
        "CD8+ T": "#2980B9", "CD4+ T": "#5DADE2", "Treg": "#85C1E9",
        "Macrophage": "#16A085", "Plasma": "#7F8C8D", "Endothelial": "#27AE60",
        "Pericyte": "#95A5A6"}


def build_labels(obs):
    sub = obs["annotation_subtypes"].astype(str).values
    maj = obs["annotation_majortypes"].astype(str).values
    lab = maj.copy()
    caf = maj == "CAF"
    lab[caf] = sub[caf]
    lab[maj == "Malignant"] = "Malignant"
    return lab


def cliffs_delta(a, b):
    """Cliff's delta 효과크기 (a=CRT, b=Untreated). 범위 [-1,1]."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    gt = sum((x > y) for x in a for y in b)
    lt = sum((x < y) for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def boot_ci(x, n=5000, seed=0):
    x = np.asarray(x, float)
    if len(x) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    bs = [np.mean(rng.choice(x, len(x), replace=True)) for _ in range(n)]
    return tuple(np.percentile(bs, [2.5, 97.5]))


def main():
    a = ad.read_h5ad(H5, backed="r")
    obs = a.obs
    lab = build_labels(obs)
    X = obs["CenterX_global_px"].astype(float).values * PX_UM
    Y = obs["CenterY_global_px"].astype(float).values * PX_UM
    samp = obs["sample_id"].astype(str).values
    treat = obs["treatment_status"].astype(str).values

    per_sample = {}
    for s in pd.unique(samp):
        m = samp == s
        tr = treat[m][0]
        if tr not in ("Untreated", "CRT"):
            continue
        rows = run_rim_panel(np.column_stack([X[m], Y[m]]), lab[m],
                             tumor="Malignant", shell_um=SHELL, n_perms=300)
        per_sample[s] = dict(treat=tr, z={r["cell_type"]: r["z"] for r in rows})
        print(f"[{s:8} {tr:9}] " + "  ".join(
            f"{c}:{per_sample[s]['z'].get(c, float('nan')):+.0f}"
            for c in CTS if c in per_sample[s]["z"]), flush=True)

    naive = [d for d in per_sample.values() if d["treat"] == "Untreated"]
    crt = [d for d in per_sample.values() if d["treat"] == "CRT"]
    present = [c for c in CTS if sum(c in d["z"] for d in per_sample.values()) >= 6]
    print(f"\nAnalysis unit = patient. Untreated n={len(naive)}, CRT n={len(crt)}")

    # 세포타입별 patient-level 통계
    stat = []
    for c in present:
        nz = [d["z"][c] for d in naive if c in d["z"]]
        cz = [d["z"][c] for d in crt if c in d["z"]]
        if len(nz) < 3 or len(cz) < 3:
            continue
        try:
            U, p = mannwhitneyu(cz, nz, alternative="two-sided")
        except ValueError:
            p = np.nan
        stat.append(dict(cell_type=c,
                         untr_mean=np.mean(nz), untr_ci=boot_ci(nz),
                         crt_mean=np.mean(cz), crt_ci=boot_ci(cz),
                         delta=cliffs_delta(cz, nz), p=p,
                         n_untr=len(nz), n_crt=len(cz)))
    # Benjamini-Hochberg
    ps = np.array([s["p"] for s in stat], float)
    order = np.argsort(np.where(np.isnan(ps), 1.0, ps))
    m = len(ps)
    padj = np.full(m, np.nan)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        k = m - rank
        val = min(prev, ps[idx] * m / k) if not np.isnan(ps[idx]) else np.nan
        padj[idx] = val
        if not np.isnan(val):
            prev = val
    for s, pa in zip(stat, padj):
        s["p_bh"] = pa

    print(f"\n{'cell type':<12}{'Untr mean[95%CI]':>22}{'CRT mean[95%CI]':>22}"
          f"{'delta':>8}{'p':>8}{'p_BH':>8}")
    for s in stat:
        print(f"  {s['cell_type']:<12}"
              f"{s['untr_mean']:>7.1f}[{s['untr_ci'][0]:.0f},{s['untr_ci'][1]:.0f}]"
              f"{s['crt_mean']:>9.1f}[{s['crt_ci'][0]:.0f},{s['crt_ci'][1]:.0f}]"
              f"{s['delta']:>+8.2f}{s['p']:>8.3f}{s['p_bh']:>8.3f}")
    df = pd.DataFrame([{**{k: v for k, v in s.items() if k not in ('untr_ci', 'crt_ci')},
                        "untr_lo": s["untr_ci"][0], "untr_hi": s["untr_ci"][1],
                        "crt_lo": s["crt_ci"][0], "crt_hi": s["crt_ci"][1]}
                       for s in stat])
    os.makedirs(os.path.join(ROOT, "data", "scotia"), exist_ok=True)
    df.to_csv(os.path.join(ROOT, "data", "scotia", "scotia_rim_stats.csv"), index=False)
    figure(stat, naive, crt)


def figure(stat, naive, crt):
    present = [s["cell_type"] for s in stat]
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    y = np.arange(len(present)); h = 0.34
    for s, yy in zip(stat, y):
        c = CMAP.get(s["cell_type"], "#999")
        # group means with 95% CI, patients as points
        ax.errorbar(s["untr_mean"], yy + h/2,
                    xerr=[[s["untr_mean"] - s["untr_ci"][0]],
                          [s["untr_ci"][1] - s["untr_mean"]]],
                    fmt="o", color=c, alpha=0.55, ms=7, capsize=3, mec="black", mew=0.6)
        ax.errorbar(s["crt_mean"], yy - h/2,
                    xerr=[[s["crt_mean"] - s["crt_ci"][0]],
                          [s["crt_ci"][1] - s["crt_mean"]]],
                    fmt="s", color=c, alpha=1.0, ms=7, capsize=3, mec="black", mew=0.6)
    # patient points
    for grp, off, al in [(naive, +h/2, 0.4), (crt, -h/2, 0.5)]:
        for d in grp:
            for c, yy in zip(present, y):
                if c in d["z"]:
                    ax.scatter(d["z"][c], yy + off + np.random.default_rng(
                        abs(hash(c)) % 999).uniform(-0.05, 0.05), s=9,
                        color="black", alpha=al, zorder=1)
    # significance stars (BH-adjusted)
    for s, yy in zip(stat, y):
        if not np.isnan(s.get("p_bh", np.nan)):
            mark = "**" if s["p_bh"] < 0.01 else ("*" if s["p_bh"] < 0.05 else "")
            if mark:
                xm = max(s["untr_ci"][1], s["crt_ci"][1])
                ax.text(xm + 0.6, yy, mark, fontsize=13, va="center", color="black")
    ax.set_yticks(y); ax.set_yticklabels(present); ax.invert_yaxis()
    ax.axvline(0, color="black", lw=1)
    ax.axvline(2, color="#888", ls="--", lw=1); ax.axvline(-2, color="#888", ls="--", lw=1)
    ax.set_xlabel("Peritumoral rim (30 µm) enrichment z  (>2 adjacent, <-2 excluded)  ·  "
                  "points = individual patients")
    ax.set_title("Patient-level peritumoral composition: Untreated (n=9) vs CRT (n=6)\n"
                 "group mean ± 95% CI (patient bootstrap); "
                 "* BH-adjusted p<0.05, ** p<0.01; cross-sectional association",
                 fontsize=11, fontweight="bold")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#888",
               markeredgecolor="k", alpha=0.55, markersize=9, label="Untreated (n=9)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#888",
               markeredgecolor="k", markersize=9, label="CRT-treated (n=6)")],
        frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    out = os.path.join(ROOT, "assets", "rim_scotia_stats.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
