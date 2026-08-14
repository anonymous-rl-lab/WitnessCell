# Frozen probe log

## Wessels official-split pilot

- Official GEARS/scPerturBench test conditions matched exactly for all seeds:
  seed 1 = 65, seed 2 = 55, seed 3 = 67.
- Across 187 seed-condition instances, the global-amplitude Witness correction
  reduced MSE versus the factorized-only ablation by about 30% on average.
- Partial leaderboard against the 15 published genetic methods:
  - top-100 Pearson + MSE: position 1, mean rank 3.39;
  - top-5000 Pearson + MSE: position 2, mean rank 2.36;
  - top-100 Pearson + MSE + energy + symmetric-KL: position 1, mean rank 4.30;
  - top-5000 four-metric pilot: position 2, mean rank 4.73.

These are directional gates, not a formal SOTA claim.  The formal result still
requires all four combination datasets, the frozen DEG artifact, Wasserstein,
DEG recovery, and the published aggregate transformation.

## Formal four-dataset result

- All 12 CPU predictions completed on Norman, Wessels, Schmidt, and
  Replogle_exp6 using official split seeds 1--3.
- All 6,540 expected metric rows are present and finite; every top-100
  Wasserstein solve converged.
- The published top-100 and top-5000 aggregate tables were reconstructed with
  exact ranks before WitnessCell was inserted.
- Official top-100 six-metric mean per-operation rank: WitnessCell 5.1437
  (position 1), scouter 5.3578 (position 2), linearModel 5.5581 (position 3).
- WitnessCell remains position 1 in the 5,000-gene and both-panel robustness
  aggregations. The dataset-equal top-100 sensitivity analysis ranks it second.
- Paired perturbation-cluster bootstrap: improvement over linearModel is stable
  at 95%; the smaller lead over scouter is not statistically resolved.

## Rejected probe: witness-count-stratified amplitude

A validation-only amplitude fitted separately for targets with one versus two
witnessed endpoints was tested once.  It worsened top-100 two-metric position
from 1 to 2 and did not improve top-5000 position.  The branch was rejected and
the simpler pooled closed-form amplitude remains frozen.
