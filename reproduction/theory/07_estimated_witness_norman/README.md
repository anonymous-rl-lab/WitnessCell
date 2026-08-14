# Gate 07 — estimated Witness Risk on Norman

This gate asks whether the oracle quantities from Gate 06 survive estimation
from training double perturbations and whether they improve held-out
prediction.

The order is fixed:

```bash
# known-covariance positive and negative controls
python validate_estimator.py --seeds 12

# five-seed engineering smoke test
python run_norman_estimated_witness.py --seeds 5 --seed-offset 0 --out results/smoke

# untouched formal seeds after the smoke test and protocol freeze
python run_norman_estimated_witness.py --seeds 30 --seed-offset 100 --out results/formal_30split

# frozen-result numerical audit and seed-level confidence intervals
python audit_formal_results.py --results results/formal_30split
python summarize_formal.py --results results/formal_30split --bootstrap 20000
```

`build_results_workbook.mjs` uses the bundled `@oai/artifact-tool` runtime to
rebuild the formatted XLSX summary from the frozen CSVs.  Set
`WITNESSCELL_PREVIEW_DIR` to retain rendered QA previews.

The estimated model is compared with a strong nested-CV geometry-only kernel
ridge baseline and a leakage-allowed empirical oracle.  The oracle is used
only to measure risk alignment and the remaining performance ceiling.

Each outer split also writes `estimated_covariances/seed_XXX.npz`, containing
the exact deployed `K_hat`, `k_t_hat`, `k_tt_hat`, noise estimate, target risk,
indices, pair names, and selected hyperparameters.

Formal gate, frozen before seeds 100–129 are opened:

1. mean within-split Spearman between estimated and oracle risk >= 0.20;
2. mean within-split Spearman between estimated risk and realized MSE >= 0.25;
3. held-out residual MSE improves by >= 3%, with >= 70% split wins and
   one-sided paired Wilcoxon p < 0.05;
4. mean full-effect MSE is no worse than geometry-only.

The shipped pseudobulk is a development panel.  A publication-level claim
still requires replay on independently split cells so measurement noise and
biological signal can be separated.
