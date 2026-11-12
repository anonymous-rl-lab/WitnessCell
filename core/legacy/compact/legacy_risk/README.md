# WitnessCell E20–E23 compact evidence replay

This module preserves the frozen results and protocols for preregistered computational replication (E20), selective prediction (E21), metric calibration (E22), and joint prediction/risk assessment (E23). It does not train a new model or revise a scientific verdict.

Run from any directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run.py --out /absolute/path/to/new-output
```

Requires Python 3.10+, NumPy, pandas, SciPy and scikit-learn. `REPLAY_REPORT.json` contains regenerated scores, confidence intervals, selection checks and the maximum numerical difference. Source-loss moments are float64; no lossy rounding or float16 conversion is used.

## What is reproducible

| Experiment | Compact input | Executed verification |
| --- | --- | --- |
| E20 | 654 condition-level paired loss rows, frozen protocols and 156-array identity receipt | Recompute all/dataset/subgroup MSE and PCC summaries; retain active/fallback outcomes. The original tar's 156 arrays were compared byte-for-byte to canonical E19 predictions during compaction. |
| E21 | Unchanged small experiment: calibration/risk/outcome rows, thresholds, protocols and inference arrays | Recompute selective risks and 20,000 pair-cluster bootstrap / within-seed permutation draws with original seeds 20260811/20260812. |
| E22 | 2,616 method-by-condition loss/correlation moment rows, metric-control rows, Gate21 weighted moments, frozen split identities | Reconstruct all five per-operation scalar metrics; recompute metric-control summaries, original 20,000 condition-cluster bootstrap (seed 20260822), and Gate21's original weighted-risk inference. |
| E23 | 2,334 query/method/scope loss-moment rows; frozen risks, OOF records, selections, all five comparator loss rows, metadata | Reconstruct full/non-target MSE, weighted MSE, relative loss, cosine loss, top-100 loss and weighted delta R². Recompute all 5,000-draw canonical-pair bootstrap distributions (master seed 1729), baseline comparisons, and all 6,224 frozen acceptance flags. |

E22's current-paper combination-only result uses 540 combination/seed units. It is kept distinct from the older 654-unit mixed single/combination report. The `trainMean` source label denotes the source-locked measured-single additive/mean-imputation comparator; it must not be relabelled as a constant mean predictor.

E21's unweighted 213-row result and E22's weighted 207-row complete-case result remain distinct. E22's complete-case corrections A1/A2 are retained; full-population evaluability is not silently declared successful. E23's unresolved Scouter comparison and pending full prediction-competitiveness gate are retained.

## Precision and source validation

`BUILD_VERIFICATION.json` records the audit performed against decoded v24 expression, prediction and weight arrays. E22's original frozen vectors were independently rescored, including NIR. E23's frozen prediction rows were exactly matched to the canonical E19 predictor, and all loss fields were regenerated with the original gene masks and top-100 tie ordering.

`SOURCE_INPUTS_SHA256.json` identifies the original files consumed to derive compact inputs. `ORIGINAL_RESAMPLING_ARRAYS_SHA256.json` records the original complete bootstrap/permutation arrays. The replay checks exact hashes when unchanged. E22's original resampling reference arrays are retained to permit element-wise tolerance checks for floating-point differences introduced by the original CSV round trip; the tolerance is 1e-12 and actual differences are reported. No bootstrap draws or confidence intervals are replaced with selected or favourable samples.

## Explicit scope limits

This is exact frozen-metric and uncertainty replay from compact sufficient statistics. It is not raw-cell preprocessing, model retraining, or independent biological replication. NIR identities and continuous DRF scores are retained at the per-condition level; current replay reaggregates them. E22 NIR was independently rescored from decoded vectors during the build audit. E20 array-identity verification is retained as a provenance receipt; redoing that byte identity comparison requires the external original result tar and canonical E19 vectors.

Original scientific source scripts and frozen protocol files remain for inspection. Their historical full-pipeline entry points can require external raw/contracts/training-cache assets that were removed from the compact package. Use this module's `run.py` for the self-contained replay. Historical manifest hashes identify the historical releases; they are not a manifest of this compact directory.

The producer `build_from_v24.py` is provided for audit. It requires the external original v24 repository and its reproduction store; it is not part of the ordinary replay. The later compact experiments E25–E32 are not changed by this module. E25's canonical Wessels dependencies are retained separately by the package integrator.
