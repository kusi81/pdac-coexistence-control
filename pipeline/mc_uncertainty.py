"""Monte Carlo epistemic 불확실성 — compound 가정 자체를 샘플링.

리뷰어(Major 3)/PM 지적: 30-seed는 architecture/stochastic(aleatoric) 불확실성만
개선한다. 더 중요한 것은 epistemic 불확실성 — compound effect coefficient·exposure
weight·bioavailability·synergy가 assigned value라면, 그 값이 틀리면 30번 반복해도
결론은 같게 틀린다. 이 스크립트는 새로 정의된 API 기준으로 이 가정들을 분포에서
샘플링(Monte Carlo)해 regimen ranking의 강건성을 검정한다.

draw마다 샘플:
  - effect coefficient factor : 각 API 효능계수의 편차를 로그정규(±)로 스케일
  - bioavailability           : API별 사전분포(Table S2 접지: SAC 높음…curcumin 매우낮음)
                                → 효능을 곱으로 약화(생체이용률 불확실성)
  - exposure weight factor    : 독성/노출 가중치를 로그정규로 섭동(독립)
  - synergy                   : Uniform[0, 0.3] (0 또는 불확실성 분포)
+ architecture seed = draw id (aleatoric도 일부 포함)

각 draw에서 모든 regimen을 동일 가정으로 돌려 (control, exposure) 랭크 → ensemble
전체의 rank 분포·progression-free 비율로 ranking stability 보고.
결과: data/mc_uncertainty.csv, assets/mc_uncertainty.png.
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
import abm
from abm import simulate, control_metrics
from synthetic import make_tissue

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARAMS = dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10,
              init_resistant_frac=0.03, mutation_rate=0.003,
              resistant_immune_evasion=0.35, max_tumor=2800)
DAYS = 120
ADAPT_ON, ADAPT_OFF = 1.1, 0.7
CAP = 1500
NDRAW = 100
WORKERS = 6
EFF_CV = 0.35            # 효능계수 불확실성(로그정규 CV)
WEIGHT_CV = 0.40         # 노출가중치 불확실성

# API-기반 regimen (새로 정의된 활성성분)
REGIMENS = [
    ("SAC", [("sac", 1.0)]),
    ("Eupatilin", [("eupatilin", 1.0)]),
    ("20(S)-Rg3", [("rg3_20s", 1.0)]),
    ("Curcumin", [("curcumin", 1.0)]),
    ("SAC+Eupatilin", [("sac", 1.0), ("eupatilin", 1.0)]),
    ("Eupatilin+Rg3", [("eupatilin", 1.0), ("rg3_20s", 1.0)]),
    ("SAC+Rg3+Curcumin", [("sac", 1.0), ("rg3_20s", 1.0), ("curcumin", 1.0)]),
    ("Gemcitabine", [("gemcitabine", 1.0)]),
]
COMPOUNDS = ["sac", "eupatilin", "rg3_20s", "curcumin", "gemcitabine"]
# 생체이용률 사전분포 Beta(a,b) — Table S2 접지(경구): SAC 높음 … curcumin 매우 낮음
BIOAV = {"sac": (9.0, 1.0), "eupatilin": (5.0, 3.0), "rg3_20s": (2.0, 6.0),
         "curcumin": (1.5, 10.0), "gemcitabine": (50.0, 1.0)}  # gem: IV≈1
_TIS = {}


def _tissue(seed):
    if seed not in _TIS:
        c, l, _ = make_tissue("contained", seed=seed)
        idx = np.random.default_rng(1000 + seed).choice(len(l), CAP, replace=False)
        _TIS[seed] = (c[idx], l[idx], int((l[idx] == "Tumor").sum()))
    return _TIS[seed]


def _lognorm(rng, cv):
    s = np.sqrt(np.log(1 + cv * cv))
    return float(np.exp(rng.normal(-0.5 * s * s, s)))     # 중앙값 1


def evaluate(draw):
    """draw별로 compound 가정을 샘플→모든 regimen 실행→per-regimen 지표 리스트."""
    rng = np.random.default_rng(10_000 + draw)
    # 1) 이 draw의 compound 가정 샘플
    eff_factor, weight_factor, bioav = {}, {}, {}
    for name in COMPOUNDS:
        eff_factor[name] = _lognorm(rng, EFF_CV)
        weight_factor[name] = _lognorm(rng, WEIGHT_CV)
        a, b = BIOAV[name]
        bioav[name] = float(rng.beta(a, b))
    synergy = float(rng.uniform(0.0, 0.30))
    # 2) 모듈 dict 임시 변경(효능=편차×eff×bioav, 노출가중치×weight)
    saved_eff = {n: dict(abm.SUBSTANCES[n]["effects"]) for n in COMPOUNDS}
    saved_tox = {n: abm.TOXICITY[n] for n in COMPOUNDS}
    try:
        for name in COMPOUNDS:
            fac = eff_factor[name] * bioav[name]     # 효능 스케일(생체이용률 반영)
            abm.SUBSTANCES[name]["effects"] = {
                p: 1.0 + (m - 1.0) * fac for p, m in saved_eff[name].items()}
            abm.TOXICITY[name] = saved_tox[name] * weight_factor[name]
        seed = 1 + (draw % 8)                         # architecture도 일부 변동
        c, l, n0 = _tissue(seed)
        out = []
        for label, subs in REGIMENS:
            p = dict(PARAMS); p["seed"] = seed
            h, _ = simulate(c, l, days=DAYS, params=p, regimen_subs=subs,
                            synergy=synergy, adaptive=True,
                            adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF,
                            snapshots=(0.0, 1.0))
            m = control_metrics(h, n0=n0)
            out.append([label, float(m["cum_toxicity"]),
                        float(m["progression_censored"]), float(m["final_frac"])])
    finally:
        for name in COMPOUNDS:                        # 복원
            abm.SUBSTANCES[name]["effects"] = saved_eff[name]
            abm.TOXICITY[name] = saved_tox[name]
    # 3) draw 내 랭크(통제 우선, 그다음 노출)
    order = sorted(range(len(out)),
                   key=lambda i: (-out[i][2], out[i][1]))
    rank = {out[i][0]: r + 1 for r, i in enumerate(order)}
    return [(draw, lab, expo, cens, fin, rank[lab]) for lab, expo, cens, fin in out]


def main():
    from concurrent.futures import ProcessPoolExecutor
    print(f"{NDRAW} draws × {len(REGIMENS)} regimens = {NDRAW*len(REGIMENS)} sims, "
          f"{WORKERS}워커 병렬 (Monte Carlo epistemic)", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for k, res in enumerate(ex.map(evaluate, range(NDRAW), chunksize=2)):
            rows.extend(res)
            if (k + 1) % 20 == 0 or k == NDRAW - 1:
                print(f"  {k+1}/{NDRAW} draws", flush=True)
    df = pd.DataFrame(rows, columns=["draw", "regimen", "exposure",
                                     "censored", "final", "rank"])
    df.to_csv(os.path.join(ROOT, "data", "mc_uncertainty.csv"), index=False)

    print(f"\n{'regimen':<20}{'ctrl%':>7}{'expo med[95%]':>20}{'rank med[95%]':>16}{'top3%':>7}")
    labels = [r[0] for r in REGIMENS]
    for lab in labels:
        g = df[df.regimen == lab]
        el, em, eh = np.percentile(g.exposure, [2.5, 50, 97.5])
        rl, rm, rh = np.percentile(g["rank"], [2.5, 50, 97.5])
        print(f"{lab:<20}{g.censored.mean()*100:>6.0f}%"
              f"{em:>9.0f}[{el:.0f}-{eh:.0f}]"
              f"{rm:>8.0f}[{rl:.0f}-{rh:.0f}]{(g['rank']<=3).mean()*100:>6.0f}%")
    figure(df)


def figure(df):
    labels = [r[0] for r in REGIMENS]
    order = sorted(labels, key=lambda L: np.median(df[df.regimen == L]["rank"]))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14.5, 5.8))
    cols = plt.cm.viridis(np.linspace(0.1, 0.9, len(order)))

    # a: rank stability across MC ensemble (epistemic)
    data = [df[df.regimen == L]["rank"].values for L in order]
    bp = a1.boxplot(data, vert=False, patch_artist=True, widths=0.6,
                    medianprops=dict(color="black", lw=1.4),
                    flierprops=dict(marker=".", ms=3, alpha=0.3))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    a1.set_yticklabels(order, fontsize=8.5); a1.invert_yaxis()
    a1.set_xlabel("Per-draw rank (1 = best: controlled at lowest exposure)")
    a1.set_title("a  Ranking stability under compound-assumption uncertainty\n"
                 f"(Monte Carlo over effect, exposure weight, bioavailability, synergy; "
                 f"n={NDRAW})", fontsize=10, fontweight="bold")

    # b: progression-free fraction (control robustness) across ensemble
    ctrl = [df[df.regimen == L]["censored"].mean() * 100 for L in order]
    y = np.arange(len(order))
    a2.barh(y, ctrl, color=cols, edgecolor="black", linewidth=0.6, alpha=0.85)
    for i, v in enumerate(ctrl):
        a2.text(v + 1, i, f"{v:.0f}%", va="center", fontsize=8)
    a2.set_yticks(y); a2.set_yticklabels(order, fontsize=8.5); a2.invert_yaxis()
    a2.set_xlim(0, 108)
    a2.set_xlabel("Draws controlling the tumor (progression-free), %")
    a2.set_title("b  Control robustness to compound uncertainty\n"
                 "(fraction of Monte Carlo draws progression-free)",
                 fontsize=10, fontweight="bold")

    fig.suptitle("Epistemic (compound-assumption) uncertainty by Monte Carlo over the "
                 "API-resolved regimens — sampling effect coefficients, exposure weights, "
                 "bioavailability, and synergy (in silico)", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "assets", "mc_uncertainty.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    figure(pd.read_csv(os.path.join(ROOT, "data", "mc_uncertainty.csv")))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
