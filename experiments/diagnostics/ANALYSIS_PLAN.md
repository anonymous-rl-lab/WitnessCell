# WitnessCell v37: fixed retrospective attribution analysis

This plan was written on 2026-11-09 before computing the two new controls. Han outcomes and earlier results were already known. This is retrospective method attribution, not preregistration or an independent confirmation. Original models, outcomes, folds and published scores remain unchanged.

## Scope

1. Recompute fixed-direction oracle amplitudes, with both nonnegative and [0,1] constraints, from archived outcomes. These outcome-informed quantities are diagnostics only. Sum fold-specific SSE before percentages; never infer an oracle from pooled gains and amplitude ranges.
2. Tabulate training-side endpoint-eliminated repeat correlations and held-out repeat correlations separately. No correlation is labelled an information-theoretic ceiling. Shared controls can inflate apparent repeatability.
3. Test attribution of Han SF12 gains using existing five seeds and three whole-pair folds, both complete readout panels, and the two nonzero SF12 source functions (HeLa, 293T). A549 remains the original zero-function control. Original V4 results remain reported but are not reinterpreted as twelve-feature models.

## Random-direction reference

Use 999 fixed Gaussian coefficient vectors (seed 20261109), independently for each source. Draw in the source training projected-feature standardized coordinates, giving each vector the same coefficient norm as its fitted source vector. Within each allowed calibration set, rescale its projected output to the fitted source direction's energy, using no outcomes. Fit the original clipped [0,1] amplitude from that calibration set, then evaluate on the original query set. Use the same vector across folds and response definitions. Report the entire distribution, empirical exceedance counts and per-seed results. This chosen random reference is a specificity diagnostic, not a universal biological null or an exchangeability-based population p-value.

## Nonlinear references

Use the same target calibration pairs as the original scalar procedure. Single-effect-derived features have no gene annotations or source outcome weights.

Primary: additive baseline plus a saturating single-sum correction and fitted endpoint effects. Saturation is z/(1+lambda*abs(z)/r), r is the RMS of the permitted single responses. Select lambda in {0,.05,.1,.2,.5,1,2,5} by three inner whole-pair folds, complete squared error, favouring lower lambda for exact ties. Fit the saturation coefficient and endpoint effects by least squares on allowed calibration pairs, separately by repeat.

Sensitivity: additive baseline plus the four symmetric local SF12 response terms formed with target single-response RMS, and fitted endpoint effects. Fit standardized projected features using ridge in {.01,.1,1,10}; select by the same inner complete-response loss, favouring stronger ridge for ties.

Refit each selected nuisance baseline on the outer calibration set. Fit a common clipped partner coefficient for each unchanged frozen source direction against that baseline's projected calibration residual. Evaluate baseline alone and baseline plus the controlled direction on the untouched query fold. Retain original additive results. These are comparator/nuisance fits and scalar fits; no source biological predictor is retrained. Report baseline improvement and additional direction gain with distinct denominators. No choice between these references is made using query scores.

## Reporting and validation

All new controls are retained regardless of sign. Query outcomes enter scoring and oracle diagnostics only, never selection. Record exact source hashes, splits, coefficients, summed losses, random draws and numerical identities. Cross-check selected scores by direct array residualization. Preserve the complete WitnessCell evidence structure: theory, released core, K562/Jurkat, RNA, then control and Han attribution. Avoid causal claims that source noise alone explains target transfer.
