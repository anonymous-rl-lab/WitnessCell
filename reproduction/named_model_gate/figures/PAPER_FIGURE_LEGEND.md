# Paper figure legend

## Main figure

**Target-conditioned witness geometry improves strict-seen1 combinatorial
prediction.** **A,** Frozen evaluation design. GEARS, CPA and WitnessCell used
the same source-half-A single perturbations, controls, 12 training doublets and
5 validation doublets. Five disjoint destination-half-B doublets calibrated the
Witness amplitude and model-specific scalar fusion weights. Fourteen
destination-half-B strict-seen1 pairs were reserved for final scoring. The
formal matrix comprised three independent cell splits and three model seeds per
named model (18 GPU fits); inference used the 14 biological target pairs.
**B,** Mean interaction-residual MSE for six frozen strategies. Points are pair
means after averaging the nine paired scoring runs; error bars are pair-
bootstrap 95% confidence intervals. **C,** Pair-level residual MSE for
WitnessCell against the saturation construction control. All 14 targets lie
below the identity line (42.49% mean reduction; exact one-sided sign-flip
`p=6.10e-5`). **D,** Relative residual-MSE reductions with pair-bootstrap 95%
confidence intervals. The open point marks the only interval crossing zero.

## Extended diagnostic figure

**Named-model complementarity is small after Witness calibration.** **A,**
Full-effect cosine for the six frozen strategies; points and intervals use the
same pair-level aggregation as the main figure. **B,** Calibration-selected
Witness fusion weights across nine paired runs. GEARS+Witness retains a mean
8.17% GEARS contribution, whereas CPA+Witness retains a mean 2.10% CPA
contribution. Horizontal bars denote means.

