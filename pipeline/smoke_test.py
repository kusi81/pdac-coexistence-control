"""Headless smoke test: prove the core metrics separate contained vs diffuse.

Run:  python pipeline/smoke_test.py
Expects: contained myCAF z high & positive, diffuse myCAF z near 0.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    sys.stdout.reconfigure(encoding="utf-8")   # Korean Windows console is cp949
except Exception:
    pass

from synthetic import make_tissue
from analysis import run_core_metrics, barrier_verdict, squidpy_available


def summarize(mode):
    coords, labels, _ = make_tissue(mode=mode, seed=42)
    res = run_core_metrics(coords, labels, n_perms=500)
    print(f"\n=== {mode} ({len(labels)} cells) ===")
    pt = res["proximity_test"]
    if pt:
        print(f"  proximity: myCAF {pt['median_dist_a']:.1f}um vs "
              f"iCAF {pt['median_dist_b']:.1f}um  p={pt['p_permutation']:.3g} "
              f"-> {pt['interpretation']}")
    for b, bs in res["barrier"].items():
        v, txt = barrier_verdict(bs)
        z = bs["z_score"]
        print(f"  barrier {b}: z={z:+.1f}  obs={bs['observed_corridor_density']:.2f} "
              f"null={bs['null_mean']:.2f}  -> {v}")
    return res


if __name__ == "__main__":
    print(f"squidpy available: {squidpy_available()}")
    c = summarize("contained")
    d = summarize("diffuse")
    cz = c["barrier"].get("myCAF", {}).get("z_score", 0) or 0
    dz = d["barrier"].get("myCAF", {}).get("z_score", 0) or 0
    print(f"\nVALIDATION: contained myCAF z={cz:+.1f}, diffuse myCAF z={dz:+.1f}")
    ok = cz > 2 and abs(dz) < cz / 2
    print("RESULT:", "PASS ✅ metrics separate architecture from abundance"
          if ok else "CHECK ⚠️ unexpected separation")
