# v1.2 two-gate formal SOTA release

This version organizes the completed SOTA evidence into two independent,
re-runnable experiment units without changing any model output or scientific
number.

## Experiment 1: four-dataset CPU prediction

- 4 datasets × 3 official split seeds = 12 fits.
- 654 held-out condition-seed predictions.
- Exact agreement with published test conditions for every dataset and seed.
- Target-free, pickle-free deployment archives.
- Frozen hyperparameters and validation diagnostics for every fit.
- Explicit zero-witness degradation to the factorized backbone.

## Experiment 2: official full-metric SOTA gate

- 6,540 raw metric records.
- Six official top-100 metrics and four all-gene metrics.
- 654/654 converged Wasserstein solves.
- Exact reconstruction of all published baseline ranks before insertion.
- WitnessCell official-primary position 1 with mean rank 5.1437.
- All-gene and both-panel operation-weighted robustness position 1.
- Dataset-equal top-100 sensitivity position 2.
- Paired perturbation-cluster bootstrap and explicit claim boundary.

Each experiment contains a README, runner, independent audit, and frozen
machine-readable result.
