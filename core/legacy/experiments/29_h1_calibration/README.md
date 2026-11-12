# H1 measurement calibration and complete-response prediction

WitnessCell paper v22: Results 6, Figure 7a, Supplementary Note 11 and Table S18.

This module reproduces the four-arm calibration comparison from compact
sufficient statistics and independently refits the calibration operators.
The primary result is **12.92474038% lower complete RNA MSE** for nine-library
calibration relative to raw cell-weighted means. This is a complete-response
gain, not evidence of an additional partner component.

## Run the included experiment replay

Python 3.10+ and NumPy are sufficient. No network, GPU, large raw matrix or
other experiment directory is required. From this directory:

```bash
python run.py --out /path/to/new_external_result_directory
```

The unified entry point verifies `MODULE_MANIFEST.json` and writes
`REPLAY_REPORT.json` to the new external directory. It refuses a nonempty
directory or a path inside this archived module. For a read-only stdout-only
check, `python replay.py` is also available.

The script can also be called by absolute path from another working directory.
It verifies the archived source hashes, recomputes all 168 condition/half
errors from the full-RNA score Gram matrix, repeats all 20 inner candidate
selections, and independently reconstructs all 80 preprocessing operators
with constrained weighted least squares or the specified direct pooling rule.
It checks that held-out combination coefficients are exactly zero. A mismatch
exits with an error. The reference assets remain unchanged.

`qa/REPLAY_VERIFIED.json` records this integration's executed numerical check.
The stored score Gram is sufficient for the reported Euclidean scores but
does not recover individual predictions for all 37,292 genes. Five illustrative
gene predictions are retained separately in `reference/FIVE_GENE_EXAMPLES.npz`.

## Expected results

| Input processing | Full RNA MSE | Gain relative to raw |
| --- | ---: | ---: |
| Raw nine-library cell-weighted mean | 0.001623698041296788 | 0% |
| Direct fixed seven-library composition | 0.001476556934290247 | 9.062097956% |
| Seven-library common-offset calibration | 0.001436203065359193 | 11.547404207% |
| Nine-library common-offset calibration | 0.001413839284906622 | 12.924740380% |

At the same seven-library input support, learned calibration improves on
direct composition standardization by **2.732970737%**. Nine-library calibration
improves 20/21 conditions and 5/5 outer folds; omitting any one condition leaves
an aggregate gain of approximately 11.53–17.97%. The nine-versus-seven-library
comparison also changes measurement support and is not a pure estimator test.

## Experimental identity and boundaries

- Source: Joung et al., *A transcription factor atlas of directed differentiation*,
  Cell 186, 209–229.e26 (2023), DOI 10.1016/j.cell.2022.11.026. The processed
  H1/hESC TF-overexpression cache is the Joung/Zhang 2023 source used in the
  archived experiment. Source study accession: GSE216481; it is not a claim
  that every packaged intermediate is a deposited file with that accession.
- Historical development pool: 25 combination conditions and 10 singles;
  5 outer folds, with 20 training pairs per outer fold. Twenty-one held-out
  conditions meet the common seven-library target coverage criteria.
- Each arm uses identical fixed-reference targets, generated from unmodified
  technical-half measurements and the same seven-library composition.
- Calibration, single-profile origin and complete M0 fitting are rebuilt
  within each permitted training split. Seven saturation candidates (including
  an unanchored endpoint alternative), five noise parameters and analytic
  shrinkage are selected using three inner folds.
- The H1 M0 uses the permitted-single mean as an origin because this cache
  lacks an untreated control. It is the archived native H1 predictor family,
  not the numerical v14 implementation, and the mean is not an untreated sample.
- The two halves are technical measurement references, not independent
  biological replicates. This is historical development evidence.
- The module contains no real-partner success claim from artificial injection.

## Optional rerun from full expression intermediates

`source/run_m0.py`, `source/origin_free_model_frozen.py` and
`source/replay_m0.py` are unchanged archival code. The first two implement the
original prediction experiment. Their original paths are historical; use the
portable wrapper below rather than editing those files:

```bash
python retrain_from_intermediates.py \
  --inputs /path/to/batch_preprocessor_20261020 \
  --output /path/to/new_calibration_run
```

This requires four external files listed with SHA-256 in
`reference/SOURCE_INPUT_HASHES.json`: full-depth condition-by-library expression,
technical-half expression, and their two metadata archives. They are deliberately
not bundled. The wrapper validates exact source hashes before invoking the
unchanged trainer and refuses to overwrite a nonempty output directory. The
original five-minute training bound remains in effect. The compact package has
been executed through `replay.py`; this integration did not rerun full-expression
training or original cell processing.

## Files and provenance

- `reference/`: unchanged sufficient statistics, original protocols, all
  candidates/selections, per-condition errors, summaries and five-gene examples.
- `metadata/design.npz`: exact copies of the count/design arrays needed to
  refit calibration operators; no expression or target values are added.
- `source/`: unchanged historical trainer, model and replay source.
- `SOURCE_PROVENANCE.json`: original experiment archive hash and per-file
  source hashes; documents the small metadata extraction.
- `replay.py`: independent integration checker, with no absolute workspace paths.
- `run.py`: standard external-output entry point with module checksum verification.
- `retrain_from_intermediates.py`: path-resolving wrapper for the optional
  unchanged full-expression trainer.

The original experiment archive also contains artificial injection and
counterfactual diagnostics. These are not needed to replay Figure 7a and are
not duplicated here.
