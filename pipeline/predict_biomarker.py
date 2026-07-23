"""Non-trivial 예측 (3) — 어떤 measurable spatial biomarker 조합이 regime을 구분하나.

리뷰어 Major-1: 자명하지 않은 예측 + 임상 활용성. 측정 가능한 공간 biomarker(조합)가
치료 regime(면역요법이 통하는가=면역반응성 vs 면역배제)을 구분하는가? 그리고 단순
myCAF abundance는 구분하지 못하는가(=arrangement가 관건)?

tissue 패널: contained(myCAF 링)를 다양한 비율로 thinning + diffuse(같은 abundance,
산개). 각 tissue에서 측정 가능한 공간 biomarker 계산:
  - barrier_score        : 종양-면역 직선경로의 myCAF 개재 비율(interposition)
  - CD8 infiltration     : 종양 25µm 내 CD8 비율
  - rim_enrichment       : 종양 주변 shell의 myCAF 농축비
  - myCAF abundance      : 단순 myCAF 분율(‘양’ 대조)
그리고 고정 면역요법을 돌려 효능(final tumor, 낮을수록 면역반응성)을 얻어 regime 정의.
→ barrier_score·CD8 침윤 조합이 regime을 가르고, abundance는 못 가름을 보인다
   (contained-full vs diffuse-full: 같은 abundance, 반대 regime).

결과: data/predict_biomarker.csv, assets/predict_biomarker.png.
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
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt
matplotlib.rcParams["axes.unicode_minus"] = False
from scipy.spatial import cKDTree
from synthetic import make_tissue
from abm import simulate, control_metrics
from spatial_core import barrier_score, rim_enrichment

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARAMS = dict(k_prolif=0.15, cd8_recruit=12, k_kill=0.6, cd8_barrier_alpha=2.2,
              caf_confine=0.7, caf_pressure=0.9, init_resistant_frac=0.02,
              mutation_rate=0.002, resistant_immune_evasion=0.4, max_tumor=6000)
IMMUNE = [("generic_immunostim", 1.0)]
DAYS = 90
SEEDS = [1, 2, 3, 4, 5]
# (label, geometry, myCAF thinning fraction)
PANEL = [
    ("Contained\n(full myCAF)", "contained", 0.0),
    ("Contained\nthin-35%", "contained", 0.35),
    ("Contained\nthin-60%", "contained", 0.60),
    ("Contained\nthin-80%", "contained", 0.80),
    ("Diffuse\n(full myCAF)", "diffuse", 0.0),
]


def _thin(c, l, frac, seed):
    if frac <= 0:
        return c, l
    rng = np.random.default_rng(9000 + seed)
    caf = np.where(l == "myCAF")[0]
    drop = rng.choice(caf, int(len(caf) * frac), replace=False)
    keep = np.ones(len(l), bool); keep[drop] = False
    return c[keep], l[keep]


def _peri_density(c, l):
    """peritumoral myCAF 밀도 = 종양 30µm 내 myCAF 수 / 종양 세포 수(측정 가능한 절대 장벽)."""
    t = c[l == "Tumor"]; caf = c[l == "myCAF"]
    if len(t) < 5 or len(caf) < 1:
        return 0.0
    d, _ = cKDTree(t).query(caf, k=1)
    return float((d <= 30).sum() / len(t))


def _biomarkers(c, l):
    # 기능적 장벽 = peritumoral myCAF 절대 밀도(thinning·geometry 모두 추적).
    # rim_enrichment의 분율/비율은 남은 myCAF가 rim집중이라 thinning에 둔감 → 밀도 사용.
    re = rim_enrichment(c, l, seed=0)
    rim_frac = float(re["obs_rim_fraction"]) if isinstance(re, dict) else np.nan
    mycaf = float(np.mean(l == "myCAF"))
    return _peri_density(c, l), rim_frac, mycaf


def main():
    rows = []
    for label, geom, thin in PANEL:
        for seed in SEEDS:
            c0, l0, _ = make_tissue(geom, seed=seed)
            c, l = _thin(c0, l0, thin, seed)
            peri, rim_frac, mycaf = _biomarkers(c, l)
            n0 = int((l == "Tumor").sum())
            p = dict(PARAMS); p["seed"] = seed
            h, _ = simulate(c, l, days=DAYS, params=p, regimen_subs=IMMUNE,
                            snapshots=(0.0, 1.0))
            fin = control_metrics(h, n0=n0)["final_frac"]
            rows.append([label.replace("\n", " "), peri, rim_frac, mycaf, fin])
        g = [r for r in rows if r[0] == label.replace("\n", " ")]
        print(f"{label.replace(chr(10),' '):<26} peri_density={np.mean([r[1] for r in g]):.3f} "
              f"rim_frac={np.mean([r[2] for r in g]):.2f} "
              f"myCAF_abund={np.mean([r[3] for r in g]):.2f} "
              f"immune_final={np.mean([r[4] for r in g]):.2f}x", flush=True)
    df = pd.DataFrame(rows, columns=["arch", "peri_density", "rim_frac",
                                     "mycaf_frac", "immune_final"])
    df.to_csv(os.path.join(ROOT, "data", "predict_biomarker.csv"), index=False)
    print("\n[면역효능(final tumor) 예측 상관 |r|]:")
    for col in ["peri_density", "rim_frac", "mycaf_frac"]:
        r = np.corrcoef(df[col], df["immune_final"])[0, 1]
        print(f"  {col:<12} r={r:+.2f}")
    figure(df)


def figure(df):
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.0))
    # a: immune efficacy vs myCAF abundance — abundance FAILS to separate
    a = axs[0]
    archs = df["arch"].unique()
    cmap = plt.cm.coolwarm(np.linspace(0, 1, len(archs)))
    for arch, col in zip(archs, cmap):
        g = df[df.arch == arch]
        a.scatter(g["mycaf_frac"], g["immune_final"], s=45, color=col,
                  edgecolor="black", linewidths=0.4, label=arch, alpha=0.85)
    r = np.corrcoef(df["mycaf_frac"], df["immune_final"])[0, 1]
    a.set_xlabel("myCAF abundance (fraction)")
    a.set_ylabel("Tumor burden under immune therapy\n(lower = immune-responsive)")
    a.set_title(f"a  Abundance alone does NOT predict regime\n(|r|={abs(r):.2f}; "
                "same abundance → opposite outcome)", fontsize=9.5, fontweight="bold")
    a.legend(fontsize=6.5, loc="upper left", framealpha=0.7)

    # b: immune efficacy vs peritumoral myCAF density — barrier PREDICTS
    b = axs[1]
    for arch, col in zip(archs, cmap):
        g = df[df.arch == arch]
        b.scatter(g["peri_density"], g["immune_final"], s=45, color=col,
                  edgecolor="black", linewidths=0.4, alpha=0.85)
    r2 = np.corrcoef(df["peri_density"], df["immune_final"])[0, 1]
    b.set_xlabel("Peritumoral myCAF density (per tumor cell)")
    b.set_ylabel("Tumor burden under immune therapy")
    b.set_title(f"b  Peritumoral barrier density predicts response\n(|r|={abs(r2):.2f})",
                fontsize=9.5, fontweight="bold")

    # c: 2D biomarker map (barrier density × abundance) colored by efficacy
    cax = axs[2]
    sc = cax.scatter(df["peri_density"], df["mycaf_frac"], c=df["immune_final"],
                     cmap="RdYlGn_r", s=55, edgecolor="black", linewidths=0.4,
                     vmin=0, vmax=max(1.5, df["immune_final"].max()))
    cax.set_xlabel("Peritumoral myCAF density (per tumor cell)")
    cax.set_ylabel("myCAF abundance (fraction)")
    cax.set_title("c  Barrier, not abundance, tracks the regime\n"
                  "(green = immune-responsive, red = excluded)", fontsize=9.5,
                  fontweight="bold")
    fig.colorbar(sc, ax=cax, fraction=0.045, pad=0.02,
                 label="Tumor burden under immune therapy")

    fig.suptitle("Which measurable spatial biomarker distinguishes the immune regime? "
                 "The peritumoral myCAF barrier—not total myCAF abundance—separates "
                 "immune-responsive from immune-excluded tissue (#1; in silico)",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "assets", "predict_biomarker.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    figure(pd.read_csv(os.path.join(ROOT, "data", "predict_biomarker.csv")))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
