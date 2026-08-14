# Experiment 22 pre-freeze corrections

These corrections were made while the protocol status was `FREEZE_CANDIDATE` and
before any WitnessCell prediction was scored under Experiment 22.

1. The interpolated-duplicate definition was corrected to use the locked
   implementation's `pvals_adj_df_dict` (second-half Benjamini–Hochberg-adjusted
   target-vs-rest p-values). The upstream unadjusted-p-value assignment is a
   commented-out alternative and is not executed at the locked commit.
2. No endpoint, threshold, model prediction, accepted/rejected identity or
   formal result was changed by this correction.

The correction is included in the eventual frozen manifest.

## Comparator-availability correction

The protocol already required `PRED_LINEAR=NOT_ADJUDICATED` when neither exact
raw predictions nor an exact regeneration path is available. The original C2
wording incorrectly made scientific freeze depend on that asset being
available, rather than on its availability being resolved without imputation.
C2 now passes only when the source and local assets have been audited and the
formal status has been fixed before scoring. In the current package the raw
linearModel condition means are absent and the required R PCA/ridge runtime and
processed input state are unavailable; the method is therefore explicitly
`NOT_ADJUDICATED`, never replaced by leaderboard summaries.

This correction narrows the possible verdict and cannot improve a score.

## Gate 21 random-seed correction

The original Gate 21 executable uses seed `20260811` for the pair-cluster
bootstrap and seed `20260812` for the within-seed random-selection control. The
draft protocol collapsed these into one seed. Experiment 22 now reuses both
original seeds exactly. No result had been computed when this was corrected.

## Gate 21 canonical-weight availability

Before any new weighted loss was computed, the cell-level asset audit found
that scPerturBench's Norman panel contains no ELMSAN1 condition. Consequently,
the six `ELMSAN1+MAP2K3` query instances (one of 33 pair identities) cannot
receive a canonical cell-derived target-vs-rest weight. They are not imputed.
The original 213 identities and accept/reject decisions remain sealed; the full
WMSE verdict is `NOT_ADJUDICATED`, while a separately labelled 207-row
complete-case sensitivity is permitted. This change narrows the claim.
