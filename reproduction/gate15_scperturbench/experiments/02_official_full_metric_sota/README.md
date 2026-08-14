# Experiment 2 — official full-metric SOTA gate

## Question

After inserting WitnessCell into the released scPerturBench genetic-combination
results, does it lead under the same metrics, panels, within-operation ranks,
and final aggregation used by the official benchmark?

## Frozen design

- Inputs: the 12 target-free deployment predictions from Experiment 1.
- Primary panel: top 100 perturbation-affected genes.
- Primary metrics: PCC-delta, MSE, E-distance, symmetric KL, Wasserstein, and
  Common-DEGs.
- Robustness panel: up to 5,000 genes with the four applicable metrics.
- Primary aggregation: mean per-operation `Rank` across 654 condition-seed
  units, matching the official Fig. 4 code.
- Comparators: 15 published methods with eligible genetic-combination results.
- Sensitivity: dataset-equal weighting and a paired bootstrap clustered by
  dataset and perturbation.

## Run

The scorer is CPU-only. It requires an existing JAX/OTT scoring overlay but
does not install or touch CUDA, PyTorch, GEARS, or CPA.

```bash
export SCORER_PYTHONPATH=/path/to/scorer_vendor:/path/to/base_vendor
bash run.sh /path/to/official_h5ad_directory
```

Completed seed checkpoints are reused.

## Audit

```bash
python audit.py --package-root ../..
```

The audit verifies all 6,540 raw records, panel-specific metric contracts,
finite values, 654 converged Wasserstein solves, exact reconstruction of the
published ranks, the primary leaderboard, and bootstrap outputs.

## Frozen result

WitnessCell is the official-primary point-estimate leader:

| position | method | mean rank |
|---:|---|---:|
| 1 | WitnessCell | 5.1437 |
| 2 | scouter | 5.3578 |
| 3 | linearModel | 5.5581 |

It is also first on the 5,000-gene and both-panel operation-weighted
robustness aggregations. The dataset-equal top-100 sensitivity analysis ranks
it second, and the perturbation-cluster bootstrap does not resolve the small
lead over scouter at 95% confidence. This is therefore a formal benchmark SOTA
point estimate, not evidence of uniform dominance.
