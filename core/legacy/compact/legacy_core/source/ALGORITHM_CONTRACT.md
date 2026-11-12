# v14 algorithm contract

## Frozen v13 background

The official GEARS split, v13 mean+GO+self endpoint head, saturation,
factorized backbone, validation selection, endpoint-incidence interaction
kernel, and scoring protocol are unchanged.  The incremental head is excluded
from validation selection, so `alpha`, `noise_ratio`, and `gamma` remain the
v13 values.

## Incremental amplitude head

For an unseen endpoint `g`, v14 starts from the v13 endpoint effect and adds a
one-coordinate correction at gene `g`.  Its three training-only scalar
features are:

1. the generic on-target self response among known training singles;
2. the sign-matched tail of gene `g`'s response across known perturbations;
3. the sign-matched RMS of that response fingerprint.

The correction coefficients solve top-100 residual ridge regression against
v13 with relative ridge `0.2` and bounds `[-1.5, 1.5]`.  The value `0.2` was
selected before test reveal from a training-only nested-LOO grid.

## Dual upgrade gate

For every outer training endpoint, that endpoint is removed before both the
v13 refit and inner correction fit.  The correction activates only when the
one-sided 95% lower bounds of both quantities are positive:

- relative all-gene MSE gain over v13;
- top-100 PCC delta over v13.

Top-100 MSE is recorded but is not a third activation gate.  If either gate
rejects, the split returns v13 exactly.

## Leakage boundary

- Only official training singles enter features, coefficients, ridge
  selection, and gates.
- Validation combinations retain their frozen v13 role in interaction-model
  selection; the v14 correction never enters that selection.
- Test means, test variances, subgroups, and official scores never enter fit,
  selection, or activation.
- Deployment archives contain no truth arrays.
