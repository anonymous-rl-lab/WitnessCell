# WitnessCell partner-control supplement

This package supplies a component-specific adjustment of existing response predictions. It includes all 55 eligible functional configurations and 265 archived RNA candidate/fold records, four comparison policies, repeat consistency, query-panel sensitivity, and an outcome-free prediction interface. The analysis is retrospective; it does not replace the original functional evaluations.

## Run

Python 3 with NumPy and pandas is sufficient for numerical replay. Matplotlib is needed only for plots.

```bash
python run_analysis.py
python replay_rna.py
python independent_check.py
python reproduce.py
python software/test_partner_control.py
python software/witnesscell_partner_control.py fit --input example/training_oof.npz --output example/calibration.json
python software/witnesscell_partner_control.py predict --input example/query_predictions.npz --calibration example/calibration.json --output example/controlled_predictions.npz
```

The example uses saved HeLa internal held-combination predictions to fit a new scalar and applies it to 523 saved query predictions. It does not retrain the biological predictor. `evaluation_truth.npz` is separate and is never consumed by `fit` or `predict`.

## Supply your own arrays

The fit NPZ has exactly `y`, `baseline`, `prediction`, `edges`, `folds`, and `nodes`. Response arrays have shape N by R. Edges are unique unordered pairs of integer indices in `[0,nodes)`. Each pair and all its measurements must belong to one inner fold. The caller must generate held-combination baseline and correction predictions without using that fold's measurements. The interface cannot infer fitting provenance from numeric arrays.

The predict NPZ has exactly `baseline`, `prediction`, `edges`, and `nodes`. Query responses are rejected as input. Output arrays include the baseline, original prediction, whole-correction shrinkage, endpoint correction, partner correction, and partner-only adjustment. The query panel must be specified as a batch; predictions need not be invariant to regrouping the batch. A panel with no partner directions retains the original predictor.

All arrays must share response units, ordering and preprocessing. No per-query readout masks or arbitrary row weights are supported by this interface. The original core and growth model fitting interfaces remain in the parent `software/` directory.

## Outputs

- `results/FOUR_ARM_SCORES.csv`: 55 configurations × four policies; full, endpoint and partner SSE and gains.
- `results/FOUR_ARM_BY_CONDITION.csv`: condition summaries obtained by pooling SSE before taking ratios. The 45 K562 layout configurations share evaluation pairs.
- `results/REPEAT_DECISIONS.csv`, `REPEAT_SCORES.csv`: coefficients from each repeat and every evaluation repeat, alongside pooled calibration.
- `results/QUERY_PANEL_SCORES.csv`: full, identity-halves, identity-thirds and singleton panels, evaluated on the original query graph.
- `results/ALGEBRA_CHECKS.csv`, `INTERFACE_TESTS.json`: numerical identities and executable-interface checks.
- `results/RNA_FOUR_ARM_SCORES.csv`: 265 archived RNA records × four policies, replayed from saved sufficient statistics. These supplementary projection scores retain their original model identities and do not replace the primary rectangle analysis.
- `results/FUNCTIONAL_DIAGNOSTIC_TABLE.csv`, `RNA_DIAGNOSTIC_TABLE.csv`: model, data-repeat measures where available, inner direction/amplitude, coefficient and external outcome.
- `results/INDEPENDENT_CHECK.json`, `PORTABLE_REPLAY.json`: independent SVD checks and six-table exact replay in a new directory.
- `METHOD_DERIVATION.md`: definitions, proofs and scope.

`inputs/ARCHIVED_*` retains historical decisions and regression summaries for comparison. It is not read by the operational prediction interface. Historical paths in source ledgers identify provenance; numerical replay uses the local arrays and no external absolute paths.

## Rights and source relationship

The new interface is author-owned scientific code; the package-level `RIGHTS_AND_RELEASE.md` records the pending author-approved release licence. Existing third-party and adapter licences remain separate. No new licence grant is inferred from a successful replay. Source measurements and biological predictor implementations retain their original attribution.

## Archive coverage

`results/HISTORICAL_COVERAGE.csv` maps every one of the 788 historical ledger rows to the 55 selected functional configurations, 265 matching RNA candidate/fold records, or historical context only. Historical context rows are retained rather than presented as newly recomputed four-arm comparisons. Eligibility follows the saved primary functional selection and the availability of matched inner/outer RNA sufficient statistics; it is not defined by improvement. Candidate rows and repeated archive snapshots are not independent experiments.
