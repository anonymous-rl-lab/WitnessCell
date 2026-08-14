# WitnessCell in the scPerturBench genetic-combination track

WitnessCell is the 28th entrant in the benchmark-wide inventory and the 16th
eligible entrant in scPerturBench's genetic-combination comparison. The other
15 combination-track methods are the published baselines shipped by the
benchmark; methods evaluated only in incompatible tasks are not counted as
competitors here.

The model has two layers:

1. a factorized single-perturbation backbone with a training-only population
   fallback and validation-selected saturation map;
2. a regularized incidence-kernel correction learned only from residuals of
   training double perturbations, with validation-only amplitude calibration.

For a combination with no witnessed endpoint, the correction is exactly zero.
The model therefore falls back to its factorized prediction instead of
inventing an unsupported interaction.

## Frozen formal result

On Norman, Wessels, Schmidt, and Replogle_exp6 with official split seeds 1--3,
WitnessCell ranks first by the official primary aggregation: top 100
perturbation-affected genes, six metrics, and mean per-operation rank over all
654 condition-seed units.

| position | method | mean rank |
|---:|---|---:|
| 1 | WitnessCell | 5.1437 |
| 2 | scouter | 5.3578 |
| 3 | linearModel | 5.5581 |

This is a point-ranking SOTA result, not uniform dominance. WitnessCell is
first on Norman and Wessels, second on Replogle_exp6, and ninth on Schmidt in
the top-100 dataset-wise sensitivity analysis. A perturbation-cluster bootstrap
does not resolve the small margin over scouter at 95% confidence. Read
`FORMAL_AUDIT_REPORT.md` before using a broader claim.

## Frozen workflow

1. Reproduce the GEARS 0.1.0 split, including its GO-support filter, and match
   every published test condition.
2. Run all four datasets and seeds with target-free CPU prediction.
3. Score top-100 with PCC-delta, MSE, E-distance, symmetric KL, Wasserstein,
   and Common-DEGs; score up to 5,000 genes with the four applicable metrics.
4. Reconstruct the published aggregation exactly before inserting WitnessCell.
5. Report the official operation-weighted primary rank, the all-gene robustness
   rank, dataset-equal-weight sensitivity, and perturbation-cluster bootstrap.

## Two formal experiment units

The evidence used for the SOTA decision is organized into two independently
auditable units:

1. `experiments/01_four_dataset_cpu_prediction`: 12 CPU fits, exact split
   matching, target-free deployment archives, frozen per-fit settings, and an
   independent prediction audit.
2. `experiments/02_official_full_metric_sota`: official six-metric scoring,
   exact published-rank reconstruction, primary and robustness leaderboards,
   clustered uncertainty, and an independent SOTA audit.

Each directory contains its own question, protocol, one-command entry point,
machine-readable frozen result, and audit script.

## Commands

No script installs CUDA, PyTorch, GEARS, CPA, or a model environment.

```bash
bash 00_probe.sh
bash 01_download_official_combo_data.sh /path/to/data
bash 02_run_combo_predictions.sh /path/to/data /path/to/predictions data/gears_assets

# The official distributional scorer needs JAX + ott-jax in an existing CPU
# overlay. If they live outside the active environment, point PYTHONPATH at it.
export SCORER_PYTHONPATH=/path/to/scorer_vendor:/path/to/base_vendor
bash 04_score_official_full.sh /path/to/data /path/to/predictions /path/to/scores
bash 05_aggregate_official_full.sh /path/to/scores /path/to/aggregate
```

`04_score_official_full.sh` writes one checkpoint per seed and reuses completed
seeds. E-distance uses the exact double-precision norm/dot-product identity,
which matches SciPy `cdist` to numerical precision while using multithreaded
BLAS. `--limit-conditions 1` on `score_official_full.py` is an engineering
smoke mode and is explicitly marked non-scientific.

## Important files

- `protocol.json`: frozen claim and leakage rules.
- `FORMAL_AUDIT_REPORT.md`: independent evidence and claim boundary.
- `results/formal_score/full`: 6,540 raw WitnessCell metric records and
  Wasserstein convergence audits.
- `results/formal_score/aggregate/formal_top100_leaderboard.csv`: official
  primary leaderboard.
- `results/formal_score/aggregate/formal_both_panels_leaderboard.csv`: all-gene
  robustness leaderboard.
- `results/formal_score/aggregate/paired_cluster_bootstrap.json`: uncertainty.
- `figures/figure_main_scperturbench.png`: frozen main figure.
- `integration`: adapter for the official repository layout.
- `experiments`: the two formal evidence units supporting the SOTA result.

The directional scorer remains a diagnostic only. It does not use the original
precomputed DEG panels and must never be presented as the formal SOTA result.
