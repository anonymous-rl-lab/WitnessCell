# Model card

## Model details

- Name: WitnessCell
- Package version: 0.1.0
- Algorithm contract: frozen v14
- License: Apache-2.0
- Attribution during double-blind review: WitnessCell Authors
- Output: condition-level expression means

## Intended use

WitnessCell is intended for research on unmeasured two-endpoint genetic
combinations when condition-level moments, measured singles, a training
interaction graph, and an explicit validation role are available. It is suited
to reproducible benchmarking, target-free deployment archives, and selective
accept/abstain analyses.

## Out-of-scope use

Do not use the package as a clinical decision system, a cell generator, a
causal estimator, a guarantee of perturbation safety, or a universal
uncertainty score. Do not silently reuse the retrospective Norman threshold on
another dataset. Version 0.1.0 supports at most two endpoints per condition.

## Inputs

- Unique gene names in a fixed order.
- Per-condition expression means and population variances in that order.
- Positive per-condition cell counts.
- Explicit and disjoint training/validation condition lists.
- A gene-to-GO mapping. Missing annotations are allowed and produce the exact
  mean-program behavior; no GO asset is bundled.

Raw counts, preprocessing, batch correction, cell filtering, and condition
definition remain the user's responsibility and should be documented.

## Outputs

`PredictionBatch` contains WitnessCell means/effects and the factorized
means/effects, allowing the interaction residual to be inspected. Optional risk
and decision fields are separate and never masquerade as target truth.

## Reliability and limitations

Evidence gates reduce unsupported endpoint complexity but cannot establish
external validity. Performance can degrade under unseen cell types, assay
shifts, weak or biased single-perturbation measurements, poor GO coverage,
disconnected interaction graphs, or endpoint naming mismatches. Correlated
genes are not independent replicates; the frozen gates operate over
leave-one-single-out units.

The estimated Witness Risk is target-conditioned self-risk. Low risk is not a
probability of correctness; high risk should trigger abstention or additional
measurement, not automatic substitution with an unvalidated predictor.

## Evaluation expectations

Report split construction, number of training/validation/test conditions,
subgroup definitions, pair-balanced aggregation, gate activation, factorized
and witness performance, selective coverage, and the provenance of any risk
threshold. Keep final target outcomes outside fit and calibration inputs.

## Privacy and security

The package sends no telemetry and performs no network access. Model bundles
are pickle-free and integrity checked. They may still encode aggregate
biological measurements, so access control and governance remain necessary.
