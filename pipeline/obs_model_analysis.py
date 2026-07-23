"""관측모델(observation model) — 적응형 controller의 임상 실현가능성 검정.

리뷰어(임상) 지적: 현재 적응형은 실제 종양부담을 지연·오차 없이 즉시 관측하는 이상적
controller다. 임상에서는 CA19-9류 혈청 biomarker로만, 그것도 이산 측정·측정오차·
biomarker-부담 시간지연·non-secretor(~10%)·최소지속·안전규칙 하에서 관측한다.

이 스크립트는 동일 종양에서 네 arm을 비교:
  1) continuous            — 연속 최대투여 참조
  2) ideal adaptive        — 실제 부담 즉시관측(기존)
  3) observed CA19-9        — 지연·잡음·이산측정·최소지속·안전규칙 관측
  4) observed (non-secretor)— CA19-9 무정보 → 안전측 연속치료로 폴백
+ 측정주기(interval) 민감도. 핵심 질문: 노이즈/지연에도 적응형의 저노출 이점이 유지되는가.

20 seed, 6워커 병렬, max_tumor 상한. 결과: data/obs_model.csv, assets/obs_model.png.
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
from synthetic import make_tissue
from abm import simulate, control_metrics

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARAMS = dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10,
              init_resistant_frac=0.03, mutation_rate=0.003,
              resistant_immune_evasion=0.35, max_tumor=6000)
DRUG = [("generic_cytotoxic", 1.0)]
DAYS = 150
ADAPT_ON, ADAPT_OFF = 1.1, 0.7
CAP = 1800
SEEDS = list(range(1, 21))       # 20 seeds
WORKERS = 6
OBS = dict(obs_model=True, obs_interval=28.0, obs_noise_cv=0.25,
           obs_lag_days=14.0, min_on_days=14.0, min_off_days=21.0,
           obs_confirm=2, obs_safety_mult=1.4)
# arm: (label, adaptive?, extra params)
ARMS = [
    ("Continuous", False, {}),
    ("Ideal adaptive", True, {}),
    ("Observed CA19-9", True, dict(OBS)),
    ("Observed (non-secretor)", True, dict(OBS, nonsecretor=True)),
]
INTERVALS = [14.0, 28.0, 56.0]   # 측정주기 민감도(일)


def _tissue(seed):
    c, l, _ = make_tissue("contained", seed=seed)
    idx = np.random.default_rng(1000 + seed).choice(len(l), CAP, replace=False)
    return c[idx], l[idx], int((l[idx] == "Tumor").sum())


def _run(seed, adaptive, extra):
    c, l, n0 = _tissue(seed)
    p = dict(PARAMS); p.update(extra); p["seed"] = seed
    kw = dict(days=DAYS, params=p, regimen_subs=DRUG, snapshots=(0.0, 1.0))
    if adaptive:
        kw.update(adaptive=True, adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF)
    h, _ = simulate(c, l, **kw)
    return h, control_metrics(h, n0=n0)


def evaluate(task):
    kind = task[0]
    if kind == "arm":
        _, ai, seed = task
        label, adaptive, extra = ARMS[ai]
        _, m = _run(seed, adaptive, extra)
        return ("arm", ai, seed, m["cum_toxicity"], float(m["progression_censored"]),
                m["final_frac"])
    else:                             # interval sensitivity (observed, secretor)
        _, iv, seed = task
        _, m = _run(seed, True, dict(OBS, obs_interval=INTERVALS[iv]))
        return ("iv", iv, seed, m["cum_toxicity"], float(m["progression_censored"]),
                m["final_frac"])


def main():
    from concurrent.futures import ProcessPoolExecutor
    tasks = [("arm", ai, s) for ai in range(len(ARMS)) for s in SEEDS]
    tasks += [("iv", iv, s) for iv in range(len(INTERVALS)) for s in SEEDS
              if INTERVALS[iv] != 28.0]     # 28일은 arm2와 동일 → 중복 제외
    print(f"{len(tasks)} runs, {WORKERS}워커 병렬", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for k, r in enumerate(ex.map(evaluate, tasks, chunksize=4)):
            rows.append(r)
            if (k + 1) % 40 == 0 or k == len(tasks) - 1:
                print(f"  {k+1}/{len(tasks)}", flush=True)
    df = pd.DataFrame(rows, columns=["kind", "idx", "seed", "exposure",
                                     "censored", "final"])
    df.to_csv(os.path.join(ROOT, "data", "obs_model.csv"), index=False)

    print(f"\n{'arm':<26}{'expo med[IQR]':>18}{'control%':>10}")
    for ai, (label, _a, _e) in enumerate(ARMS):
        g = df[(df.kind == "arm") & (df.idx == ai)]
        lo, md, hi = np.percentile(g.exposure, [25, 50, 75])
        print(f"{label:<26}{md:>7.0f}[{lo:.0f}-{hi:.0f}]{g.censored.mean()*100:>9.0f}%")
    print("측정주기 민감도(관측 적응형, secretor):")
    for iv, val in enumerate(INTERVALS):
        if val == 28.0:
            g = df[(df.kind == "arm") & (df.idx == 2)]
        else:
            g = df[(df.kind == "iv") & (df.idx == iv)]
        print(f"  interval={val:>4.0f}d: expo med={np.median(g.exposure):.0f} "
              f"control={g.censored.mean()*100:.0f}%")

    # 대표 궤적(관측모델 기전 시연) — seed 1
    h, _ = _run(1, True, dict(OBS))
    figure(df, h)


def figure(df, h):
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.0))

    # ── Panel a: observation-model mechanics (one representative patient) ──
    a = axs[0]
    t = np.array([x["t"] for x in h])
    true = np.array([x["n_tumor"] for x in h]); true = true / max(true[0], 1)
    bio = np.array([x["biomarker"] for x in h])
    obsr = np.array([x["obs_ratio"] for x in h])
    on = np.array([x["drug_on"] for x in h])
    a.plot(t, true, color="#34495E", lw=2, label="True tumor burden")
    a.plot(t, bio, color="#2980B9", lw=1.6, ls="--", label="Latent CA19-9 (lagged)")
    # 측정 시점(obs_ratio 변화 지점) 산점
    chg = np.where(np.abs(np.diff(obsr)) > 1e-9)[0] + 1
    a.scatter(t[chg], obsr[chg], s=26, color="#E67E22", zorder=5,
              label="Measurements (noisy)", edgecolor="black", linewidths=0.4)
    # 투여 on 구간 음영
    for i in range(1, len(t)):
        if on[i]:
            a.axvspan(t[i-1], t[i], color="#C0392B", alpha=0.06, lw=0)
    a.axhline(ADAPT_ON, color="#E74C3C", ls=":", lw=1)
    a.axhline(ADAPT_OFF, color="#27AE60", ls=":", lw=1)
    a.text(2, ADAPT_ON + 0.02, "on-band", fontsize=7, color="#E74C3C")
    a.text(2, ADAPT_OFF - 0.06, "off-band", fontsize=7, color="#1E8449")
    a.set_xlabel("Time (days)"); a.set_ylabel("Burden / initial")
    a.set_title("a  Observation-model mechanics\n(lagged, noisy, interval-sampled "
                "CA19-9; red = drug on)", fontsize=10, fontweight="bold")
    a.legend(frameon=False, fontsize=7.5, loc="upper right")

    # ── Panel b: exposure by arm (does the advantage survive?) ──
    b = axs[1]
    labels = [x[0] for x in ARMS]
    cols = ["#C0392B", "#27AE60", "#2980B9", "#8E44AD"]
    data = [df[(df.kind == "arm") & (df.idx == ai)]["exposure"].values
            for ai in range(len(ARMS))]
    bp = b.boxplot(data, patch_artist=True, widths=0.6,
                   medianprops=dict(color="black", lw=1.4))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    ctrl = [df[(df.kind == "arm") & (df.idx == ai)]["censored"].mean() * 100
            for ai in range(len(ARMS))]
    for i, (d, cc) in enumerate(zip(data, ctrl)):
        b.text(i + 1, np.median(d) + 4, f"ctrl {cc:.0f}%", ha="center", fontsize=7.5)
    b.set_xticklabels(["Continuous", "Ideal\nadaptive", "Observed\nCA19-9",
                       "Observed\nnon-secretor"], fontsize=8)
    b.set_ylabel("Cumulative modeled exposure  →  lower better")
    b.set_title("b  Adaptive advantage under realistic observation\n"
                "(20 seeds; survives but attenuates; non-secretor loses it)",
                fontsize=10, fontweight="bold")

    # ── Panel c: measurement-interval sensitivity ──
    cax = axs[2]
    xs, meds, los, his = [], [], [], []
    for iv, val in enumerate(INTERVALS):
        if val == 28.0:
            g = df[(df.kind == "arm") & (df.idx == 2)]["exposure"].values
        else:
            g = df[(df.kind == "iv") & (df.idx == iv)]["exposure"].values
        lo, md, hi = np.percentile(g, [25, 50, 75])
        xs.append(val); meds.append(md); los.append(md - lo); his.append(hi - md)
    cax.errorbar(xs, meds, yerr=[los, his], fmt="o-", color="#2980B9", lw=2,
                 capsize=3, label="Observed CA19-9 adaptive")
    ideal = df[(df.kind == "arm") & (df.idx == 1)]["exposure"].median()
    cont = df[(df.kind == "arm") & (df.idx == 0)]["exposure"].median()
    cax.axhline(ideal, color="#27AE60", ls="--", lw=1.4, label="Ideal adaptive")
    cax.axhline(cont, color="#C0392B", ls="--", lw=1.4, label="Continuous")
    cax.set_xlabel("CA19-9 measurement interval (days)")
    cax.set_ylabel("Median exposure")
    cax.set_title("c  Longer intervals erode the benefit\n"
                  "(exposure rises toward continuous)", fontsize=10,
                  fontweight="bold")
    cax.set_xticks(INTERVALS)
    cax.legend(frameon=False, fontsize=8, loc="center right")

    fig.suptitle("Adaptive control under a realistic CA19-9 observation model "
                 "(delay, noise, interval sampling, non-secretors, safety rules) — in silico",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "assets", "obs_model.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    df = pd.read_csv(os.path.join(ROOT, "data", "obs_model.csv"))
    h, _ = _run(1, True, dict(OBS))
    figure(df, h)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
