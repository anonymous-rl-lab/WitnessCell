# Gate 21 sandbox — frozen selective prediction

This sandbox tests whether training-only estimated Witness Risk can support a
frozen accept/abstain decision on target identities that were absent from the
calibration panel.

The roles are fixed:

- `response`: held-out per-target squared error of the estimated-Witness
  predictor;
- `risk`: `estimated_witness_risk`, computed without the held-out double
  outcome;
- `decision`: accept when risk is no greater than a threshold frozen from the
  five-seed calibration panel; otherwise abstain. Geometry-only fallback is a
  secondary diagnostic, not the primary claim.

The five engineering seeds (0–4) from Gate 07 are the only calibration data.
The once-only test consists of formal seeds 100–129 restricted to perturbation
pairs never appearing in calibration. Test reliability, oracle risk and test
loss are sealed until the protocol and thresholds are hashed.

Execution order:

1. `prepare_assets.py`
2. `run_calibration_smoke.py`
3. `freeze_gate.py`
4. `reveal_formal_test.py`
5. `audit.py`

This is a retrospective frozen reanalysis of an existing Norman development
panel. It is not a prospective biological validation or a model-agnostic
uncertainty guarantee.
