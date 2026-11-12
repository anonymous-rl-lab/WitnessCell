# Experiment 25: support controls and witness quantity

This is the compact, completed supplementary experiment payload for
`grn_repo_v20`. It contains the original frozen analyses and their finished
results. This release reorganizes existing evidence; it adds no new scientific
experiment and changes no result or model prediction.

Both the payload directory and its standalone ZIP are below **5,000,000 bytes**.
The standalone payload reproduces the reported selection and bootstrap
statistics from retained per-query intermediates. Reconstructing the Wessels
endpoint states and full prediction vectors additionally uses the **canonical
assets in the main code package**, identified by SHA256 in
`CANONICAL_INPUT_MAP.json`. The small ZIP is not an independent copy of those
large assets.

## Reproduce

Run from the extracted `grn_repo_v20` repository root. Restore the main package's
losslessly compressed data to their canonical paths first:

```bash
python tools/restore_data.py --verify
python experiments/25_support_controls_and_witness_quantity/verify_bundle.py
python experiments/25_support_controls_and_witness_quantity/reproduce.py \
  --repo . --output /tmp/witnesscell_followup_replay --mode all --max-seconds 600
```

Choose a new output directory on each run. Inputs, reference outputs and
historical locks are read-only. The output directory contains fresh logs,
reconstructed numerical files, rerun confidence intervals and
`reproduction_status.json`. The wrapper enforces a total subprocess wall limit
and limits BLAS/OpenMP to one thread.

For a standalone statistics replay, with no main repository or raw cells:

```bash
python reproduce.py --output /tmp/witnesscell_followup_statistics \
  --mode statistics --max-seconds 600
```

The other modes are `seen` and `witness`. `witness` requires `--repo` and performs
the full endpoint-to-prediction-to-loss-to-bootstrap replay. `statistics` uses
the saved Wessels query-loss rows and does not reconstruct expression vectors.

The reference runtime was Python 3.12.14 and NumPy 2.3.5. NumPy is the only
dependency for these replay commands; `requirements_reproduction.txt` preserves
the NumPy pin. The supplied endpoint reconstruction code is NumPy-based.
Fresh floating-point computation on another runtime is not automatically
guaranteed to have identical bits: the wrapper verifies the original `.npy`
member bytes and archive hashes and fails on a discrepancy. It never rounds
numbers or relaxes the historical locks to obtain a pass.

## What is retained and what is reconstructed

| Asset | Treatment |
| --- | --- |
| Original analysis scripts, plans, source locks and chronology | Retained unchanged |
| Original OOF risk training rows, frozen features, coefficients and accepted identities | Retained as full-precision CSV/JSON |
| Per-query losses, effect controls, final statistics and figure-ready tables | Retained as full-precision CSV/JSON |
| Wessels split/order/fraction subset identities and frozen prediction inventories | Retained unchanged |
| Wessels allowed condition aggregates, evaluation truth/weights, full replay parameters and GO assets | Reuse one canonical copy in the main package |
| Three frozen numerical states and 21 unique prediction archives | Reconstruct from the canonical inputs using the original source functions |
| Both 5,000-draw bootstrap archives | Reconstruct using the original scripts and original PCG64 seeds |
| Old raw cell matrices, repeated copies of evaluation NPZs and generated plotting binaries | Not included in this supplement |

`RECONSTRUCTED_ASSETS.json` records each omitted NPZ's original file hash and
every original `.npy` member hash. Those member hashes preserve array dtype,
shape and all numerical bits. Compression and deterministic regeneration
replace repeated storage; neither quantization nor a reduced-precision array
was introduced. The 29 omitted numerical archives total 62,251,590 bytes before
ZIP-level packaging. The combined package retains the inputs required to
recreate them.

`reproduce.py` maps old absolute paths to hash-verified canonical locations and
calls the archived numerical functions without editing them. For the original
unmodified verifier, historical metadata files are copied into the fresh replay
workspace **only after** the reconstructed numerical files match their original
hashes. Fresh timing/metadata are separately retained in
`fresh_reconstruction_metadata.json`; this does not create or impersonate a
new preregistration or historical experiment freeze.

The retained `build_review_outputs.py` and the two original numerical scripts
are archival executed sources. Use the portable wrapper for the compact
release; some historical default paths and old inventory entries intentionally
remain in those audit records.

To regenerate the historical descriptive tables and summary figure, copy this
payload to a new working directory, install the extra dependencies in
`requirements_plotting.txt`, and run `python build_review_outputs.py` in that
copy. The script uses the retained reference CSV/JSON inputs and writes its
`tables/` and `figures/` outputs in the working copy. Its old source-path lookup
falls back to the identical retained `original_step1_context/` files.

## Evidence and current interpretation

The primary seen-control comparison fixes the final v14 prediction, excludes
the two perturbed coordinates, and delivers `floor(0.5*n)` queries in each seed
batch. `B` uses predicted effect plus the eight-member condition-subsampling
ensemble variance. `S` adds missing-endpoint and inverse-degree support
features. `R` orders seen2, seen1 and seen0 strata, then uses `B` within a
stratum. `C` jointly fits the same base features with seen0/seen2 indicators;
`C0` omits the degree feature. All fitted risks use the same original weighted
ridge procedure on training-side out-of-fold records.

At 50% coverage, S versus C gives:

| Comparison | Relative reduction in non-target MSE, 95% pair-bootstrap interval |
| --- | --- |
| Norman, all queries | −2.09% [−11.23%, 10.84%] |
| Wessels, all queries | 8.27% [1.74%, 15.30%] |
| Wessels, seen1 queries | 17.57% [5.52%, 26.11%] |
| Wessels, predicted-effect quartile controlled | 6.68% [−0.44%, 10.02%] |

The Norman primary comparison does not establish additional value beyond the
joint seen control. Wessels provides conditional positive evidence, including
seen1, but the effect-stratified interval includes zero. These findings do not
establish superiority to all output-based uncertainty measures, nor do they
rule out every effect-size explanation. The detailed all-coverage and metric
results remain in `seen/evaluation/results.json` and the accompanying CSVs.

The Wessels witness experiment keeps final endpoints, the factorized backbone,
hyperparameters and evaluation queries fixed, and varies only the residual
witness library. It uses three original splits, three predetermined nested
orders and 25%, 50% and 100% of the available combination witnesses. MSE gains
against the fixed factorized predictor rise from 19.10% to 27.42% to 33.25%; the
paired 100% minus 25% gain is 14.15 percentage points [12.42, 16.06]. This is a
component intervention in the fixed system. It is not an end-to-end lower-data
retraining experiment, a prospective experiment-acquisition result, or evidence
that witness density is the sole reason Norman and Wessels differ.

Repeated gene pairs carry all their seed records in the 5,000 shared cluster
bootstrap draws. Subset orders, cells and overlapping seeds are not independent
biological replicates. Both datasets were previously inspected public data.
The witness bootstrap implementation was locked after its point estimates had
been inspected; the original chronology and explicit retrospective status are
preserved. The matched strong-competitor prediction-vector and prospective
biological-validation gaps remain outside this completed supplement.

## Packaging verification

`verification/FULL_REPLAY.json` records the bounded packaging replay. It
reconstructed three endpoint states, all 21 distinct expression-vector files
(6,360,340 coordinates), 1,629 query-loss rows and 543 comparator-loss rows. All
state and prediction archive bytes matched their original hashes. Both
bootstrap arrays also reproduced exactly, including all 482 `.npy` members in
the seen archive and all eight in the Wessels archive.

All 24 new frozen score, decision, model and loss files matched byte-for-byte.
The seen replay compared 10,096 numeric fields in 480 result groups plus two
effect-stratification groups, including 962 confidence intervals. The Wessels
replay compared every reported summary field (173 numeric values). Maximum
numerical differences were zero. The full replay completed in 10.41 seconds
under its 600-second bound in the reference runtime. These checks verify the
compact representation and reproducibility; they add no new independent
scientific evidence.

`verification/STATISTICS_ONLY_REPLAY.json` separately verifies the standalone
statistics mode without canonical repository inputs; it passed in 4.04 seconds.
`verification/VERIFICATION_BINDING.json` binds both reports to the adapter
source hash.

## Files to read first

`KEY_RESULTS.json`, `RESULTS_AND_DECISIONS_CN.txt` and
`LOCKED_REVISION_PLAN_CN.txt` give the completed experiment decisions.
`analysis_plan.json`, `seen/analysis_execution_lock.json`,
`seen/all_new_decisions_frozen.json`, `wessels_witness/protocol_freeze.json` and
`wessels_witness/statistics/analysis_freeze.json` preserve the actual locks and
their scope. `manuscript_source_corrections/` records the source and terminology
corrections, including the shared Wessels/CaRPool origin and Cas13d technology.
