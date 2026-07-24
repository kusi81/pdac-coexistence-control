"""전역 민감도분석 (Sobol) — OAT를 넘어서는 리뷰어 방어.

리뷰어 지적: ±50% OAT(한 번에 하나)는 상호작용·비선형을 놓친다. 전역(Sobol)
분산기반 지표로 (a) 각 파라미터가 결과 분산에 기여하는 1차 지수 S1,
(b) 상호작용 포함 총효과 ST 를 구한다. ST≫S1 이면 그 파라미터는 상호작용을
통해 작동. 8개 파라미터에는 L1 장벽의 두 축을 명시적으로 포함:
  - cd8_barrier_alpha : 면역배제(pervasive) 강도
  - caf_confine       : 물리적 봉쇄(boundary) 강도
출력: 적응형 arm 의 control_score(주), final_frac(보조). 2시드 평균으로
확률성 완화, CAP 서브샘플로 계산량 축소. 결과: CSV + Fig S7.
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
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze

from synthetic import make_tissue
from abm import simulate, control_metrics

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# OAT와 동일한 검증된 baseline 체제 (적응형 이점이 드러남)
BASE = dict(k_prolif=0.15, cd8_recruit=10, k_kill=0.5, k_caf_activate=0.10,
            init_resistant_frac=0.03, mutation_rate=0.003,
            resistant_immune_evasion=0.35, resistance_cost=0.24,
            max_tumor=4500)      # 폭주 regime 세포 상한(런타임 bound; progression 불변)
DRUG = [("generic_cytotoxic", 1.0)]
DAYS = 130
ADAPT_ON, ADAPT_OFF = 1.1, 0.7
SEEDS = [1]            # 단일시드(N↑ 우선); CI로 강건성 보고
CAP = 1800
N = 32                 # Saltelli base sample (2의 거듭제곱)

# 8개 파라미터: 종양동역학·면역·내성 + 장벽 2종(면역배제/물리봉쇄)
PROBLEM = {
    "num_vars": 8,
    "names": ["k_prolif", "k_kill", "cd8_recruit", "k_caf_activate",
              "resistance_cost", "resistant_immune_evasion",
              "cd8_barrier_alpha", "caf_confine"],
    "bounds": [[0.075, 0.30], [0.25, 1.0], [5.0, 20.0], [0.05, 0.20],
               [0.10, 0.40], [0.20, 0.70], [0.5, 3.0], [0.2, 1.5]],
}
PRETTY = {
    "k_prolif": "tumor proliferation", "k_kill": "T-cell killing",
    "cd8_recruit": "CD8 recruitment", "k_caf_activate": "CAF activation",
    "resistance_cost": "resistance cost",
    "resistant_immune_evasion": "resistant immune evasion",
    "cd8_barrier_alpha": "immune exclusion (barrier α)",
    "caf_confine": "physical confinement",
}

_TISSUE = {}


def tissue(seed):
    """seed별 tissue를 CAP까지 서브샘플(계산량 축소)."""
    if seed not in _TISSUE:
        c, l, _ = make_tissue("contained", seed=42)
        idx = np.random.default_rng(seed).choice(len(l), CAP, replace=False)
        _TISSUE[seed] = (c[idx], l[idx])
    return _TISSUE[seed]


def evaluate(x):
    """파라미터 벡터 x → (control_score, final_frac), 적응형 arm, 2시드 평균."""
    p = dict(BASE)
    for name, val in zip(PROBLEM["names"], x):
        p[name] = float(val)
    cs, ff = [], []
    for seed in SEEDS:
        c, l = tissue(seed)
        pp = dict(p); pp["seed"] = seed
        h, _ = simulate(c, l, days=DAYS, params=pp, regimen_subs=DRUG,
                        adaptive=True, adapt_on=ADAPT_ON, adapt_off=ADAPT_OFF,
                        snapshots=(0.0, 1.0))
        m = control_metrics(h, crit_mult=1.5)
        cs.append(m["control_score"]); ff.append(m["final_frac"])
    return float(np.mean(cs)), float(np.mean(ff))


WORKERS = 6            # 8코어 중 6워커 병렬(teardown 전 완결 위해 단축)


def main():
    from concurrent.futures import ProcessPoolExecutor
    X = sobol_sample.sample(PROBLEM, N, calc_second_order=False)
    ntot = X.shape[0]
    print(f"Saltelli 표본: {ntot} 평가 × {len(SEEDS)} seeds "
          f"= {ntot*len(SEEDS)} runs (D={PROBLEM['num_vars']}, N={N}), "
          f"{WORKERS}워커 병렬", flush=True)
    Ycs = np.zeros(ntot); Yff = np.zeros(ntot)
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for i, (cs, ff) in enumerate(ex.map(evaluate, list(X), chunksize=4)):
            Ycs[i], Yff[i] = cs, ff
            if (i + 1) % 80 == 0 or i == ntot - 1:
                print(f"  {i+1}/{ntot} 완료", flush=True)

    rows = []
    results = {}
    for label, Y in [("control_score", Ycs), ("final_frac", Yff)]:
        Si = sobol_analyze.analyze(PROBLEM, Y, calc_second_order=False,
                                   print_to_console=False)
        results[label] = Si
        for j, name in enumerate(PROBLEM["names"]):
            rows.append(dict(output=label, param=name,
                             S1=Si["S1"][j], S1_conf=Si["S1_conf"][j],
                             ST=Si["ST"][j], ST_conf=Si["ST_conf"][j]))
        print(f"\n[{label}] (S1 / ST)")
        order = np.argsort(-Si["ST"])
        for j in order:
            print(f"  {PROBLEM['names'][j]:26s} "
                  f"S1={Si['S1'][j]:+.3f}±{Si['S1_conf'][j]:.3f}  "
                  f"ST={Si['ST'][j]:+.3f}±{Si['ST_conf'][j]:.3f}")

    df = pd.DataFrame(rows)
    os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
    df.to_csv(os.path.join(ROOT, "data", "sobol_indices.csv"), index=False)
    print("\nwrote data/sobol_indices.csv")
    figure(results)
    # 판정: control_score 총효과 상위 파라미터
    st = results["control_score"]["ST"]
    top = PROBLEM["names"][int(np.argmax(st))]
    inter = [PROBLEM["names"][j] for j in range(len(st))
             if results["control_score"]["ST"][j]
             - results["control_score"]["S1"][j] > 0.10]
    print(f"[판정] control_score 최대 총효과 파라미터: {top}")
    print(f"[판정] 상호작용 뚜렷(ST-S1>0.10): {', '.join(inter) or '없음'}")


def figure(results):
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.6))
    for ax, label, title in [
            (axs[0], "control_score", "a  Control score (higher = better)"),
            (axs[1], "final_frac", "b  Final tumor burden (lower = better)")]:
        Si = results[label]
        names = PROBLEM["names"]
        order = np.argsort(Si["ST"])          # 아래→위 증가
        y = np.arange(len(names))
        s1 = np.clip(np.array(Si["S1"])[order], 0, None)
        st = np.clip(np.array(Si["ST"])[order], 0, None)
        s1c = np.array(Si["S1_conf"])[order]
        stc = np.array(Si["ST_conf"])[order]
        lab = [PRETTY[names[j]] for j in order]
        ax.barh(y + 0.18, st, height=0.34, color="#E67E22",
                xerr=stc, error_kw=dict(lw=0.8, ecolor="#7f5320"),
                label="Total-order $S_T$ (incl. interactions)")
        ax.barh(y - 0.18, s1, height=0.34, color="#2980B9",
                xerr=s1c, error_kw=dict(lw=0.8, ecolor="#1b4f72"),
                label="First-order $S_1$ (main effect)")
        ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=9)
        ax.set_xlabel("Sobol sensitivity index")
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.legend(frameon=False, fontsize=8.5, loc="lower right")
        ax.set_xlim(0, max(0.05, ax.get_xlim()[1]))
    fig.suptitle("Figure S7 — Global (Sobol) sensitivity of adaptive-therapy control to "
                 "eight parameters, incl. the two barrier axes\n"
                 "(immune exclusion vs physical confinement); "
                 f"D=8, N={N}, {len(SEEDS)}-seed mean, in silico",
                 fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(ROOT, "assets", "sobol.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


def replot():
    """저장된 지수 CSV에서 재작도 — 주 출력은 final tumor burden 하나만.
    (composite control_score는 구조적 결함으로 본문에서 삭제했으므로 보고하지 않음)."""
    import pandas as pd
    df = pd.read_csv(os.path.join(ROOT, "data", "sobol_indices.csv"))
    g = df[df.output == "final_frac"].set_index("param").loc[PROBLEM["names"]]
    names = PROBLEM["names"]
    order = np.argsort(g["ST"].values)
    y = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    s1 = np.clip(g["S1"].values[order], 0, None)
    st = np.clip(g["ST"].values[order], 0, None)
    ax.barh(y + 0.18, st, height=0.34, color="#E67E22",
            xerr=g["ST_conf"].values[order],
            error_kw=dict(lw=0.8, ecolor="#7f5320"),
            label="Total-order $S_T$ (incl. interactions)")
    ax.barh(y - 0.18, s1, height=0.34, color="#2980B9",
            xerr=g["S1_conf"].values[order],
            error_kw=dict(lw=0.8, ecolor="#1b4f72"),
            label="First-order $S_1$ (main effect)")
    ax.set_yticks(y)
    ax.set_yticklabels([PRETTY[names[j]] for j in order], fontsize=9)
    ax.set_xlabel("Sobol sensitivity index")
    ax.set_title("Global (Sobol) sensitivity of final tumor burden\n"
                 f"D=8 parameters, Saltelli N={N}; bars are bootstrap 95% CIs",
                 fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    out = os.path.join(ROOT, "assets", "sobol.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "replot":
        replot()
    else:
        main()
