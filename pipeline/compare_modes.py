"""contained vs diffuse 대조: 동일 세포 수, 다른 구조. 핵심 지표를 나란히 출력."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from synthetic import make_tissue
from analysis import run_core_metrics, barrier_verdict


def run(mode):
    coords, labels, _ = make_tissue(mode=mode, seed=42)
    return run_core_metrics(coords, labels, n_perms=1000, seed=0)


c, d = run("contained"), run("diffuse")


def cell(res, key):
    if key == "n":
        return sum(res["cell_counts"].values())
    if key == "myz":
        return res["barrier"]["myCAF"]["z_score"]
    if key == "icz":
        return res["barrier"]["iCAF"]["z_score"]
    if key == "myv":
        return barrier_verdict(res["barrier"]["myCAF"])[0]
    if key == "prox_a":
        return res["proximity_test"]["median_dist_a"]
    if key == "prox_b":
        return res["proximity_test"]["median_dist_b"]
    if key == "prox_p":
        return res["proximity_test"]["p_permutation"]
    if key == "cd8_near":
        r = res["radial"]["CD8_T"]
        return r["density_norm"].iloc[0]
    if key == "cd8_far":
        r = res["radial"]["CD8_T"]
        return r["density_norm"].iloc[-1]


print(f"{'지표':<34}{'contained':>14}{'diffuse':>14}")
print("-" * 62)
print(f"{'전체 세포 수':<32}{cell(c,'n'):>14}{cell(d,'n'):>14}")
print(f"{'myCAF 장벽 z':<33}{cell(c,'myz'):>+14.1f}{cell(d,'myz'):>+14.1f}")
print(f"{'  판정':<34}{cell(c,'myv'):>14}{cell(d,'myv'):>14}")
print(f"{'iCAF 장벽 z':<34}{cell(c,'icz'):>+14.1f}{cell(d,'icz'):>+14.1f}")
print(f"{'myCAF→종양 거리(µm)':<28}{cell(c,'prox_a'):>14.1f}{cell(d,'prox_a'):>14.1f}")
print(f"{'iCAF→종양 거리(µm)':<29}{cell(c,'prox_b'):>14.1f}{cell(d,'prox_b'):>14.1f}")
print(f"{'근접성 순열 p':<32}{cell(c,'prox_p'):>14.3g}{cell(d,'prox_p'):>14.3g}")
print(f"{'CD8 밀도 근거리(0µm)':<28}{cell(c,'cd8_near'):>14.2f}{cell(d,'cd8_near'):>14.2f}")
print(f"{'CD8 밀도 원거리(끝)':<28}{cell(c,'cd8_far'):>14.2f}{cell(d,'cd8_far'):>14.2f}")
