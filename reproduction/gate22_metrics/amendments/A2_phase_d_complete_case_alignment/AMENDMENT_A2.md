# Experiment 22 Amendment A2 — Phase D complete-case geometry alignment

After A1 reset the complete-case DataFrame index, Phase D progressed to the
geometry matched-coverage control and stopped before writing a result. The
frozen script passed the original 213-length `geometry_risk` array together
with 207 complete-case pair weights to `weighted_quantile`. The engineering
smoke had all rows evaluable, so this missingness-only branch was not triggered.

A2 applies both A1's dense-index correction and one additional positional
alignment: if and only if the quantile input has the 213-row Gate21 length and
its weights have the frozen evaluable-row length, select the exact immutable
`scoring_evaluable` mask stored in the pre-freeze gene contract.

No risk value, query identity, loss, threshold, random seed, replicate count,
criterion, prediction, truth or gene weight is changed. No Phase D output file
existed at either crash. Phase M and Phase P remain untouched.
