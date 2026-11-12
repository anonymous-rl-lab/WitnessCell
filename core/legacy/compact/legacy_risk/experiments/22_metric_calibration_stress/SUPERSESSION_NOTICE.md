# Supersession notice

Experiment 22 v2 supersedes the external draft named `METRIC_CALIBRATION_STRESS_PROTOCOL.md` with SHA-256:

`ea35a8b05d6eeaa3f8857e1203650616dde6378e935bbb62a6c27967eaaf3002`

That draft is retained only as provenance outside this executable experiment directory. It must not be copied into a formal runner or cited as the frozen protocol.

Material corrections in v2 include:

- source weights are absolute target-vs-rest scores, min–max scaled and squared—not raw t-statistics squared and sum-normalized;
- weighted delta R² replaces the noncanonical weighted-Pearson primary endpoint;
- NIR follows the locked condition-identity retrieval implementation, not a top-K DEG reciprocal-rank formula;
- DRF uses fixed positive and negative controls and has no invented universal cutoff;
- formal metric validity is candidate-blind, but its duplicate controls come from the same formal dataset;
- no two-dataset formal pilot is allowed;
- endpoint-gate substitution is shadow-only and cannot rewrite v14;
- Gate 21 reuses the exact original threshold and accepted set;
- Gate 21 requires gene-level arrays for WMSE and fails closed if only scalar MSE is available;
- exact v14→v13 fallback is kept distinct from abstention and from dataset-mean prediction.
