"""L2 환자접지 — 실제 환자 공간데이터를 ABM 초기조건으로.

리뷰어 지적: 치료 시뮬은 합성 조직에서 돌렸으므로 'spatially grounded'는 과장.
이 스크립트는 SCOTIA(CosMx) 실제 환자 시료에서 종양-포함 1500µm 창을 잘라
세포 위치·타입을 그대로 ABM 초기조건으로 사용 → 프레임워크가 실제 환자 조직 위에서
돌아감을 보인다. naive vs CRT 두 환자의 서로 다른 구조에서 서로 다른 동역학.
"""
import sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
warnings.filterwarnings("ignore")
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import fonts  # noqa
import matplotlib.pyplot as plt
matplotlib.rcParams["axes.unicode_minus"] = False
from abm import simulate, control_metrics

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5 = os.path.join(ROOT, "data", "scotia", "raw_meta_data_final.h5ad")
PX = 0.12028
WIN = 750          # ±750µm → 1500µm 창
CAP = 2500         # ABM 세포 수 상한(성능)
CORE = {"Malignant": "Tumor", "myCAF": "myCAF", "iCAF": "iCAF",
        "CD8+ T": "CD8_T", "Macrophage": "Macrophage"}
COL = {"Tumor": "#34495E", "myCAF": "#C0392B", "iCAF": "#E67E22",
       "CD8_T": "#2980B9", "Macrophage": "#16A085"}
NAT = [("curcumin", 1.0), ("garlic", 1.0), ("ginsenoside_rg3", 1.0)]
PATIENTS = [("U7-a", "Untreated (naive)"), ("T4-a", "CRT-treated")]


def patient_tissue(obs, sid):
    m = (obs["sample_id"].astype(str) == sid).values
    maj = obs["annotation_majortypes"].astype(str).values[m]
    sub = obs["annotation_subtypes"].astype(str).values[m]
    X = obs["CenterX_global_px"].astype(float).values[m] * PX
    Y = obs["CenterY_global_px"].astype(float).values[m] * PX
    lab = np.where(maj == "Malignant", "Tumor",
                   np.where(maj == "CAF", sub, maj))
    # 종양 중심 창
    tum = lab == "Tumor"
    cx, cy = np.median(X[tum]), np.median(Y[tum])
    win = (np.abs(X - cx) < WIN) & (np.abs(Y - cy) < WIN)
    # 핵심 5종만
    keep = win & np.isin(lab, list(CORE.values()))
    coords = np.column_stack([X[keep] - (cx - WIN), Y[keep] - (cy - WIN)])
    labels = lab[keep]
    if len(labels) > CAP:                      # 비례 서브샘플
        idx = np.random.default_rng(0).choice(len(labels), CAP, replace=False)
        coords, labels = coords[idx], labels[idx]
    return coords, labels


def main():
    a = ad.read_h5ad(H5, backed="r")
    obs = a.obs
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 9.5))
    for r, (sid, title) in enumerate(PATIENTS):
        coords, labels = patient_tissue(obs, sid)
        n0 = int((labels == "Tumor").sum())
        comp = {t: int((labels == t).sum()) for t in CORE.values()}
        print(f"[{sid} {title}] {len(labels)}세포, 종양 n0={n0}, "
              + " ".join(f"{t}:{comp[t]}" for t in CORE.values()))
        # (좌) 초기 환자 조직 산점도
        axL = axs[r, 0]
        for t in CORE.values():
            mm = labels == t
            axL.scatter(coords[mm, 0], coords[mm, 1], s=3, c=COL[t],
                        label=f"{t} ({comp[t]})", alpha=0.7, linewidths=0)
        axL.set_aspect("equal"); axL.set_xticks([]); axL.set_yticks([])
        axL.set_title(f"{sid} — {title}\npatient tissue as ABM initial condition (1500µm)",
                      fontsize=10, fontweight="bold")
        axL.legend(fontsize=7, markerscale=2.5, loc="upper right", framealpha=0.7)
        # (우) 무처치 vs 적응형 치료 궤적 (patient tissue 위)
        axR = axs[r, 1]
        for arm, kw, col in [("Untreated", dict(), "#7F8C8D"),
                             ("Adaptive natural combo", dict(
                                 regimen_subs=NAT, adaptive=True,
                                 adapt_on=1.1, adapt_off=0.7), "#27AE60")]:
            ys = []
            for seed in (1, 2):
                p = dict(seed=seed)
                h, _ = simulate(coords, labels, days=90, params=p,
                                snapshots=(0.0, 1.0), **kw)
                ys.append([x["n_tumor"] / n0 for x in h])
                ts = [x["t"] for x in h]
            ym = np.mean([np.interp(ts, [x["t"] for x in h], y)
                          if False else y for y in ys], axis=0)
            axR.plot(ts, ym, lw=2, color=col, label=arm)
        axR.axhline(1.5, color="#E74C3C", ls=":", lw=1)
        axR.set_xlabel("Time (days)"); axR.set_ylabel("Tumor / initial")
        axR.set_title(f"{sid} — dynamics on patient tissue", fontsize=10,
                      fontweight="bold")
        axR.legend(frameon=False, fontsize=8)
    fig.suptitle("Patient-grounded initial conditions — the control framework runs on real "
                 "SCOTIA CosMx tissue\n(naive vs CRT patients have distinct architectures "
                 "and dynamics; in-silico)", fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(ROOT, "assets", "patient_grounded.png")
    fig.savefig(out, dpi=118, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
