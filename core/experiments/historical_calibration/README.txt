WitnessCell model-matched partner gain and calibration — 2026-11-08

Start: REPORT_CN.txt
Run: python reproduce.py
Dependencies: requirements.txt (recorded runtime versions; Python version in ENVIRONMENT.json).

Scope: retrospective secondary analysis of recovered experiments, not new preregistration or new biological confirmation.
788 eligible model outcomes / 64 training configurations / 8 panels / 5 studies / 110 candidate pools. Historical input preserves all 1280 prior records. Functional calibration includes 55 configurations.

inputs/INNER_GROWTH_REPLAY.npz: training inner outcomes, M0 and regularized predictions, fold labels, graph edges, selected-zero indicators.
inputs/CALIBRATION_OUTCOMES.npz: outer outcomes, M0, original predictions and query identities.
results/CALIBRATION_DECISIONS.csv: alpha uses training inner reward and cost only.
results/REGRESSION_SENSITIVITY.csv: terminal-checkpoint and panel exclusions retained alongside the main analysis.
results/KEY_RESULT_TABLE.csv: original, whole-correction and partner-only calibration results.
results/VERIFICATION.json, PORTABLE_REPLAY.json: validation and independent sufficient-array replay.
results/FIGURE_NOTES.txt: interpretation and scope of both scientific figures.

Extraction scripts enrich.py and calibrate_growth.py require the original archive paths listed in the source manifests. Portable reproduce.py does not. common.py contains the fixed regression specification. verify.py uses independent SVD scoring.

Known limitations: overlapping historical experiments; only five studies; selected inner scores can be optimistic; heterogeneous repeat definitions and model implementations; small Shen pooled gains are not positive in both repeats. Partner-only shrinkage uses fixed query-batch identities and is not batch-invariant. No claim of surpassing all earlier repairs.

Previous package: https://chatgpt.com/library?file_id=libfile_0089fa8bf37881919f6094146d490ac8&account_id=personal
