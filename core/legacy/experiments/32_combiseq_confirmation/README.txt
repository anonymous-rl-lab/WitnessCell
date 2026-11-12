WitnessCell frozen program confirmation — complete offline reproduction

RESULT: NOT_CONFIRMED. On 374 Combi-seq drug pairs, adding the frozen G to the
identical endpoint control increases full-program error by 0.0232823% and
endpoint-orthogonal error by 0.0302401%. A correct implementation and repeatable
computation do not make this a positive scientific confirmation.

Read deliverables/WitnessCell_Frozen_Program_Confirmation_20261020_Report.txt
for the complete Chinese report and the relationship to the full WitnessCell
system. run1/ contains the original frozen predictions and all scored results.

To reproduce from the included compressed published count matrix, offline:

    python REPRODUCE.py --out /a/new/writable/output_directory

Required: Python 3.12, numpy 2.3.5, pandas 2.2.3. Matplotlib 3.10.8 is optional
for the included plotting script. ENVIRONMENT.json records the exact original
environment. No network, external dataset, Library access or old project path
is needed. The optional source-exclusion zip is not needed for model execution.

REPRODUCE.py uses unchanged model/run_confirmation.py and verifies 108 named
arrays by dtype/shape/byte hash, exact CSV file bytes, exact inner choices,
conditional ranges, leakage checks and scientific completion fields. Large
CONDITION_MEANS and FIT arrays are regenerated rather than stored redundantly.
EXPECTED_SCIENTIFIC_OUTPUTS.json stores their exact expected hashes.

REPRODUCTION_VERIFICATION.json documents one successful replay. A replay uses
already revealed outcomes: it is not another independent confirmation. The
newly emitted software access flags are logs of that replay; historical freeze
claims refer only to the untouched run1/ records. Original prediction hashes
are preserved in run1/ALL_PREDICTIONS_FROZEN.json. Repacked NPZ container bytes
can differ by ZIP metadata; verification is of every actual named array.

Cross-environment BLAS or NumPy versions may produce last-bit differences.
The verifier reports these rather than silently passing a numerical mismatch.

Files frozen before expression processing: CONFIRMATION_PROTOCOL.txt,
COHORT_AND_SPLITS.json and every path hashed by FREEZE.json. Do not modify them
to obtain another result. The reproduction wrapper, plot, report and package
manifest were prepared after scoring for delivery and do not change the model.

Primary comparison: endpoint control + G versus the same endpoint control.
The original scalar saturation + G remains an additional disclosed row.
Lambda selection remains the frozen scalar-saturation-based procedure; no
post-outcome optimization was added for the complete-system comparator.

Source: Mathur et al. Nature Communications 13,4450 (2022),
DOI 10.1038/s41467-022-32197-0, GEO GSE174695 (parent GSE174696).
This is one new-source, cross-intervention-type RNA confirmation with refitting
on allowed source training pairs. It is not prospective, not a new genetic
double-perturbation RNA replication, and not numerical parameter transfer.

All 3 runs and orientations of a query pair are hidden together. Internal
DMSO_DMSO controls and measured single drugs are allowed. 36 drug endpoints are
all measured; the target task is unseen combinations among these endpoints.
The main program uses 51 genes from the existing fixed mapping; no post-score
program or alias selection occurred.

The conditional endpoint-multiplier ranges are not unconditional biological
confidence intervals. One explicitly labelled post-score measurement diagnostic
in access_audit does not change the negative primary verdict. No run was dropped,
no program switched and no post-score model tuning was performed.

PACKAGE_MANIFEST.json hashes every delivered member other than itself.
The .zip.sha256 sidecar hashes the complete distributed archive.
