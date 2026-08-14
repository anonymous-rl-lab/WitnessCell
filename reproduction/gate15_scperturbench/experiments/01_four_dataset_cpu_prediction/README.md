# Experiment 1 — four-dataset CPU prediction gate

## Question

Can WitnessCell generate combination-perturbation predictions on every official
scPerturBench combination dataset using the exact published split contract,
without GPU training or access to held-out outcomes?

## Frozen design

- Datasets: Norman, Wessels, Schmidt, and Replogle_exp6.
- Split seeds: 1, 2, and 3.
- Total fits: 12.
- Model input: control, official training conditions, and official validation
  conditions only.
- Split reproduction: GEARS 0.1.0 simulation split after its GO-support filter.
- Output used downstream: pickle-free, target-free `deploy_predictions.npz`.
- Zero-witness rule: when no training double perturbation exists, the witness
  correction is exactly zero and the method reduces to the factorized backbone.

All four dataset split audits match the published test conditions exactly.
The experiment produces 654 held-out condition-seed predictions.

## Run

No installer is called.

```bash
bash run.sh /path/to/official_h5ad_directory
```

Optional output and asset locations:

```bash
bash run.sh /path/to/data /path/to/predictions /path/to/gears_assets
```

## Audit

```bash
python audit.py --package-root ../..
```

The audit verifies 12 archives, exact split matches, expected test counts,
pickle-free arrays, and absence of `truth` and `truth_variance` from deployment
files. Frozen per-fit settings are in `prediction_summary.csv`.

## Interpretation

Passing this gate establishes a reproducible, target-isolated prediction
artifact for formal scoring. It is necessary evidence, but it is not itself a
SOTA claim.
