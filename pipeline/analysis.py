"""
Analysis orchestration for the dashboard.

Splits the four metrics by dependency so the app degrades gracefully:

  ALWAYS available (numpy/scipy only):
    - metric 2  radial density profile + median NN distance
    - metric 2b proximity test (positive control: myCAF vs iCAF -> tumor)
    - metric 3  BARRIER SCORE -- the actual containment test

  OPTIONAL (needs anndata + squidpy):
    - metric 1  neighborhood enrichment
    - metric 4  co-occurrence vs radius

The most important claim (containment) lives entirely in the always-available
set, so the dashboard is useful even without the heavy spatial-omics stack.
"""

from __future__ import annotations
import numpy as np
import pandas as pd

from spatial_core import (
    PLATFORM_SPECS, check_platform,
    radial_profile, median_nn_distance, proximity_test, barrier_score,
    rim_enrichment,
)


def squidpy_available():
    try:
        import anndata  # noqa: F401
        import squidpy   # noqa: F401
        return True
    except Exception:
        return False


def platform_ceiling(platform):
    """Return (level, caveat) describing the strongest claim the platform allows."""
    spec = PLATFORM_SPECS[platform]
    if platform == "visium":
        level = "weak"
    elif platform == "visium_hd":
        level = "moderate"
    else:
        level = "supported"
    return level, spec["caveat"]


def run_core_metrics(coords, labels, tumor="Tumor", immune="CD8_T",
                     barriers=("myCAF", "iCAF"), corridor_um=30,
                     n_perms=1000, seed=0):
    """Metrics 2, 2b, 3 -- no squidpy required. Returns a results dict."""
    coords = np.asarray(coords, dtype=float)
    labels = np.asarray(labels)
    counts = pd.Series(labels).value_counts()

    out = {
        "cell_counts": counts.to_dict(),
        "radial": {},
        "median_nn": {},
        "proximity_test": None,
        "barrier": {},
    }

    # metric 2: radial profile + median NN distance for each barrier and immune
    for ct in list(barriers) + [immune]:
        prof = radial_profile(coords, labels, tumor, ct)
        if prof is not None:
            out["radial"][ct] = prof
        med, n = median_nn_distance(coords, labels, tumor, ct)
        out["median_nn"][ct] = {"median_um": med, "n": int(n)}

    # metric 2b: proximity positive control (myCAF should be closer than iCAF)
    if len(barriers) >= 2:
        out["proximity_test"] = proximity_test(
            coords, labels, tumor, barriers[0], barriers[1],
            n_perms=n_perms, seed=seed)

    # metric 3: barrier score -- the containment test
    for b in barriers:
        bs = barrier_score(coords, labels, tumor=tumor, barrier=b,
                           immune=immune, corridor_um=corridor_um, seed=seed)
        if bs:
            out["barrier"][b] = bs

    return out


def run_squidpy_metrics(coords, labels, platform="xenium",
                        radius=50, n_perms=1000, seed=0):
    """Metrics 1 & 4. Returns {} if the optional stack is unavailable/fails."""
    if not squidpy_available():
        return {}
    try:
        from synthetic import to_anndata
        from spatial_core import neighborhood_enrichment, cooccurrence
        adata = to_anndata(coords, labels)
        spec = PLATFORM_SPECS[platform]
        out = {}
        try:
            out["nhood_zscore"] = neighborhood_enrichment(
                adata, n_perms=n_perms,
                radius=radius if spec["single_cell"] else None, seed=seed)
        except Exception as e:  # noqa: BLE001
            out["nhood_error"] = str(e)
        try:
            occ, interval, cats = cooccurrence(adata)
            out["cooccurrence"] = {"occ": occ, "interval": interval,
                                   "categories": cats}
        except Exception as e:  # noqa: BLE001
            out["cooccurrence_error"] = str(e)
        return out
    except Exception as e:  # noqa: BLE001
        return {"error": str(e)}


def run_rim_panel(coords, labels, tumor="Tumor", shell_um=30.0, n_perms=1000,
                  seed=0):
    """모든 (종양 아닌) 세포타입의 종양-rim 농축을 계산. '누가 juxtatumoral인가'.

    반환: z 내림차순 정렬된 [{cell_type, z, enrichment_ratio, p, n_in_shell, verdict}].
    """
    import numpy as np
    labels = np.asarray(labels)
    cts = [c for c in pd.unique(labels) if c != tumor]
    rows = []
    for ct in cts:
        r = rim_enrichment(coords, labels, tumor=tumor, barrier=ct,
                           shell_um=shell_um, n_perms=n_perms, seed=seed)
        if not r or r.get("z_score") is None:
            continue
        z = r["z_score"]
        verdict = ("종양 인접(농축)" if z > 2 else
                   "종양에서 배제" if z < -2 else "중립")
        rows.append({"cell_type": ct, "z": z,
                     "enrichment_ratio": r["enrichment_ratio"],
                     "p": r["p_value"], "n_in_shell": r["barrier_in_shell"],
                     "verdict": verdict})
    rows.sort(key=lambda d: d["z"], reverse=True)
    return rows


def barrier_verdict(bs):
    """Human verdict from a barrier_score result dict."""
    if not bs or bs.get("z_score") is None:
        return "insufficient", "세포 수 부족 — 판정 불가"
    z = bs["z_score"]
    if z > 2:
        return "interposed", "장벽이 우연 이상으로 개재됨 (containment 지지)"
    if z < -2:
        return "anti", "오히려 경로에서 배제됨 (장벽 아님)"
    return "none", "기질 밀도 이상의 개재 근거 없음"
