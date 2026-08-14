# WitnessCell scPerturBench formal audit

## Frozen conclusion

WitnessCell is the rank-1 point estimate under the official primary genetic-
combination aggregation: top 100 perturbation-affected genes, six metrics, and
the mean per-operation rank over all four eligible combination datasets and
three official split seeds.

This is a track-specific benchmark result, not a claim against every method in
every scPerturBench task. The published benchmark contains 27 methods overall;
15 published methods have eligible results in this genetic-combination track.
Adding WitnessCell makes 16 entrants in the comparison.

## Packaged evidence units

The formal evidence is separated into two independently executable gates.
`experiments/01_four_dataset_cpu_prediction` establishes the exact splits,
12 target-isolated CPU predictions, and zero-witness fallback behavior.
`experiments/02_official_full_metric_sota` consumes those frozen predictions
and establishes the official full-metric ranking, robustness results, and
clustered uncertainty. This separation prevents engineering success from being
misreported as SOTA before official scoring is complete.

## Reproducibility identity

- Official repository commit:
  `6e24e7a9827e55d4567d2139427be9af0d1e7a6c`.
- Datasets: Norman, Wessels, Schmidt, and Replogle_exp6.
- Split seeds: 1, 2, and 3.
- Test condition-seed units: 654.
- Raw WitnessCell metric records: 6,540.
- Panels: top 100 genes (six metrics) and up to 5,000 genes (four metrics).
- Official split reproduction includes the GEARS 0.1.0 GO-support filter. All
  four datasets and all three seeds match the published test-condition table.
- Submission-side `deploy_predictions.npz` archives are target-free and
  pickle-free. Diagnostic archives attach test outcomes only after inference;
  test outcomes are never passed to the fit/predict interface. The formal
  scorer accepts the target-free deployment archive directly.

## Scoring audit

The scorer implements PCC-delta, MSE, Euclidean E-distance, log-transformed
symmetric Gaussian KL, OTT-JAX entropic Wasserstein, and top-100 Common-DEG
overlap. Distributional prediction follows the benchmark's Gaussian/control-
variance baseline contract.

The published aggregation was independently reconstructed before adding
WitnessCell. Across 9,810 published rows per panel, all official ranks were
reproduced exactly. Composite-score error was numerical only (maximum below
`1e-15`). The audit also identified that the downloadable raw table floors
negative PCC values at zero; signed PCC values must therefore be recovered
from the official aggregate table.

The exact Euclidean E-distance was accelerated with the identity
`||x-y|| = sqrt(||x||^2 + ||y||^2 - 2 x^T y)` in double precision. A real
Wessels condition produced identical four-decimal records for all 10 reported
panel/metric entries relative to SciPy `cdist`.

All 654 top-100 Wasserstein solves converged. Every expected row was present,
and no formal metric value was NaN or infinite.

## Official primary leaderboard

Lower mean rank is better.

| position | method | mean per-operation rank |
|---:|---|---:|
| 1 | **WitnessCell** | **5.1437** |
| 2 | scouter | 5.3578 |
| 3 | linearModel | 5.5581 |
| 4 | baseReg | 6.5642 |
| 5 | trainMean | 6.6070 |
| 6 | GEARS | 7.5245 |
| 7 | scGPT | 7.7278 |
| 8 | scELMo | 7.8287 |
| 9 | AttentionPert | 7.9740 |
| 10 | GeneCompass | 8.2783 |
| 11 | GenePert | 8.4511 |
| 12 | CPA | 9.2294 |
| 13 | biolord | 9.7416 |
| 14 | scFoundation | 11.6575 |
| 15 | baseMLP | 13.0031 |
| 16 | baseControl | 13.6529 |

WitnessCell's metric positions are E-distance 1, symmetric KL 1, MSE 2,
PCC-delta 2, Wasserstein 9, and Common-DEGs 9. Its lead is therefore a balanced
multimetric result rather than domination of every individual metric.

## Robustness and uncertainty

- The operation-weighted 5,000-gene panel also ranks WitnessCell first.
- Combining both panels operation-wise also ranks WitnessCell first (mean rank
  4.4725).
- Giving every dataset and panel equal weight ranks WitnessCell first overall,
  but giving every dataset equal weight on top-100 alone ranks it second.
- Top-100 dataset positions are Norman 1, Wessels 1, Replogle_exp6 2, and
  Schmidt 9. The official aggregate is consequently not evidence of uniform
  dominance in every experimental environment.
- A paired 20,000-replicate bootstrap clustered by dataset and perturbation
  estimates a WitnessCell-minus-scouter rank difference of -0.2141, with 95%
  interval [-0.6574, 0.2312] and 82.2% bootstrap probability that WitnessCell
  is better. Against linearModel, the difference is -0.4144 with interval
  [-0.7896, -0.0316] and 98.3% probability of improvement.

## Claim boundary

The defensible claim is: **WitnessCell attains the best point estimate in the
official scPerturBench genetic-combination leaderboard and remains first in the
all-gene robustness panel.** The margin over scouter is not statistically
decisive under perturbation-cluster resampling, and Schmidt remains a clear
failure environment. “Universal AIVC SOTA” or “better on every dataset and
metric” is not supported.

Machine-readable evidence is in `results/formal_score/full` and
`results/formal_score/aggregate`.
