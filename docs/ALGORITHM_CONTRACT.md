# Frozen algorithm contract

## Scope

WitnessCell 0.1.0 predicts expression **condition means** for control, single
endpoints, and two-endpoint combinations. The default configuration is
identified as `witnesscell-v14-frozen`. It is a productized transcription of
the frozen v14 research path; experimental branches are not public options.

## Fit roles and leakage boundary

`train_conditions` may contain control, measured singles, and measured doubles.
`validation_conditions` may contain validation doubles used to select the
saturation coefficient, interaction-kernel noise ratio, and residual amplitude.
The two sets must be disjoint. Final evaluation or deployment-target outcomes
must not appear in either set and are not arguments to `predict`.

All endpoint-head weights and endpoint activation decisions use training
singles only. Their leave-one-single-out evidence gates use condition-level
moments rather than treating cells as independent replicates.

## Endpoint program

For an unseen endpoint, the v13 baseline head constructs:

- the mean effect over known training singles;
- a GO-neighbor program minus that mean (top 10, Jaccard similarity); and
- a one-coordinate self anchor at the target gene.

Nonnegative weights are fit by pooled leave-one-single-out least squares and
clipped to `[0, 2]`.

The dense head is active when its one-sided 95% lower confidence bound on
all-gene MSE gain is positive, or when the joint sparse head is active. The
joint sparse head requires both:

- positive one-sided 95% lower bound on all-gene MSE gain; and
- positive one-sided 95% lower bound on top-100 Pearson-correlation delta.

There is no third requirement that sparse must beat dense. A failed gate
returns the exact simpler head: mean only, or dense without sparse.

## Incremental response-fingerprint amplitude

The v14 correction contains three one-coordinate features: self statistic,
sign-conditioned 10th/90th response-fingerprint quantile, and signed RMS.
Weights use nested leave-one-single-out fitting with ridge `0.2` and clipping to
`[-1.5, 1.5]`.

Activation requires positive one-sided 95% lower bounds for both all-gene MSE
gain and top-100 Pearson-correlation delta over v13. Top-100 MSE is recorded as
a diagnostic only and is not an activation gate. Failure returns v13 exactly.

## Combination backbone and interaction witness

Endpoint effects are added and transformed elementwise as

`effect / (1 + alpha * abs(effect))`.

Training-double residuals relative to the uncorrected v13 factorized backbone
are transferred through the unsigned endpoint-incidence kernel
`X X^T / 2`. Validation doubles select `alpha`, kernel noise ratio, and a
closed-form residual amplitude `gamma` clipped to `[0, 1]`. At inference, the
v14 endpoint head forms the factorized target backbone and the learned witness
residual is added only for two-endpoint targets.

If no training doubles exist, no interaction residual is added. If no
validation doubles exist, the contract uses the frozen uncalibrated defaults
defined in the source; diagnostics explicitly report unavailable validation
MSE as non-finite internally and serialize it safely.

## Witness Risk and selective action

`exact_witness_risk` exposes the universal-kriging/BLUP adequacy-plus-geometry
formula. `WitnessRiskEstimator` implements the frozen Gate 07 training-only
descriptor and nested-CV kernel selection. The resulting score is
**WitnessCell self-risk only**. It is not universal predictive uncertainty and
must not be used as evidence that a different model is correct.

`SelectivePolicy` emits only `accept` and `abstain`. The bundled Norman
threshold is retrospective to that frozen evaluation and has no transport
guarantee; new datasets should calibrate their own threshold under a declared
coverage and weighting protocol.

## Non-claims

The contract does not provide cell-level samples, distributional prediction,
causal identification, arbitrary higher-order combinations, clinical advice,
or a geometry-only automatic fallback/router.
