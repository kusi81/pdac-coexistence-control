# Spatial Layer (Layer 4) — CAF Architecture & Containment Test

Tests whether CAFs are **physically interposed** between tumor and T cells —
the question bulk composition (Layer 3b) cannot answer.

**Status:** validated on synthetic tissue with known ground truth. Not yet run
on real data (no network access to data hosts in the build environment).

## Install

```bash
pip install squidpy scanpy anndata   # squidpy 1.8.3 verified
```

## Run

```bash
python3 run_spatial.py \
  --adata data/sample.h5ad \
  --platform xenium \
  --outdir results
```

Input `.h5ad` needs:
- `adata.obsm['spatial']` — coordinates **in microns**
- `adata.obs['cell_type']` — categorical labels

Options: `--platform visium|visium_hd|xenium|cosmx`, `--tumor`, `--immune`,
`--barriers myCAF,iCAF`, `--radius 50`, `--n-perms 1000`.

## The four metrics, and why there are four

| # | Metric | Question | Directional? |
|---|---|---|---|
| 1 | Neighborhood enrichment | Who sits next to whom | No |
| 2 | Radial profile + proximity | How does density vary with distance | Partly |
| 3 | **Barrier score** | **Is CAF mass BETWEEN tumor and T cells** | **Yes** |
| 4 | Co-occurrence vs radius | At what length scale | No |

**Only metric 3 tests containment.** Enrichment and co-occurrence are symmetric —
they say "these are near each other" with no notion of *interposition*. Reporting
metric 1 and calling it containment would overclaim. This is the single most
important design point in the layer.

### How the barrier score works

For each (T cell → nearest tumor cell) pair, walk the connecting segment and
count barrier cells within a 30 µm corridor. Compare to a null that **relabels
an equal number of random stromal cells** as pseudo-barrier.

Matching the *count* is what makes it a test of **positioning rather than
abundance**. In a tissue that is 80–90% stroma, any test not matched on count
will call everything a barrier.

## Validation performed

Two synthetic tissues, **identical cell counts**, different architecture:

| Tissue | myCAF n | corridor obs/null | z | Verdict |
|---|---|---|---|---|
| Contained (rings + CD8 excluded) | 1099 | 2.86 / 2.10 | **+22.3** | interposed ✓ |
| Diffuse (random placement) | 1100 | 1.48 / 1.45 | +0.9 | no evidence ✓ |

- iCAF correctly shows **no** barrier in both (z = −16.7, −1.4)
- Positive control reproduced: myCAF 34.7 µm vs iCAF 191.9 µm from tumor
  (p=0.005), matching the literature expectation that myCAFs are more proximal
  to malignant cells than iCAFs
- Radial CD8 exclusion: 18× far/near in contained vs 0.16× in diffuse
- Graceful degradation on missing cell types and Visium platform confirmed

**One bug found and fixed during validation.** The first barrier implementation
used a binary "does any point on the path touch a CAF" rule. In dense stroma this
saturates at ~0.97 with a null of ~1.00, producing *negative* z-scores and zero
discriminative power — it would have reported "no containment" in a tissue built
with explicit rings. Replaced with continuous corridor density.

## Platform determines what you may claim

| Platform | Resolution | Containment claim |
|---|---|---|
| Visium | 55 µm spots, 1–10 cells | **Weak** — sub-spot geometry invisible |
| Visium HD | 8 µm bins | Moderate |
| Xenium | single cell | **Supported** (panel-limited, ~300–500 genes) |
| CosMx | single cell | **Supported** (check FOV edge effects) |

The module prints the ceiling for the chosen platform. For Visium, `cell_type`
must come from deconvolution (RCTD / cell2location / CARD / SPOTlight) and is an
argmax approximation.

## Interpreting real-data output

| Observation | Meaning |
|---|---|
| myCAF z > 2, iCAF z ≈ 0 | Containment supported; matches literature |
| Both z ≈ 0 | No architecture beyond stromal density |
| myCAF not closer than iCAF | **Check cell typing before believing it** — this contradicts a well-replicated finding |
| Everything z > 2 | Null likely mis-specified; check count matching |

## Honest limitations

1. **Straight-line paths.** Real T-cell migration follows ECM fiber tracks and
   vasculature. This is a geometric proxy, not a migration model.
2. **Static snapshot.** Containment is dynamic; fixed tissue cannot show flux.
3. **Correlation.** A barrier association does not prove the barrier *causes*
   exclusion — depletion experiments (which historically *worsened* outcomes)
   are the causal test.
4. **Cell typing dominates everything.** Every metric is downstream of the
   labels. Bad annotation produces confident, wrong geometry.

## Connecting to earlier layers

The full chain: Moffitt subtype (Layer 3) → CAF composition (Layer 3b) →
spatial architecture (Layer 4). The cross-layer question worth asking is whether
basal-like tumors show higher myCAF barrier scores, which would link the
transcriptional subtype to a physical mechanism of immune exclusion.
