WitnessCell_Qualified_Rectangles_Fixed_Test_20261026

Start: python replay_compact.py
Requires Python3.10+ and NumPy; executed with Python3.12.

Self-contained compact replay
- All five tasks, all three methods, all failures retained.
- Float64, no quantization.
- Each case_R coefficient is reproduced from all saved episode normal equations.
- Every partner SSE is recomputed directly from unique measured RNA targets and frozen predicted contrast vectors.
- Full RNA scores use per-query saved losses and quadratic case_R evaluation statistics; all full vectors were independently scored before packaging.
- All5, excludingQ4, and edge-balanced aggregates are reproducible.
- Full RNA vectors and raw-cell inputs are NOT inside this compact archive. Their original hashes remain in ALL_PREDICTIONS_FROZEN.json. It would be incorrect to describe this compact replay as raw-cell/full-feature reconstruction.

Full reconstruction sources
1. WitnessCell_Neural_Partner_Training_Inputs_20261021.zip
   Library libfile_ee62d07de8c881919088359f4ae3b89c
   TRAINING_BATCH_ONLY.npz sha2566901f8aa807ccb42de387fce2ba62ac9f40fba1ba918bab3e040edcb0e62767e
2. WitnessCell_Case_Transfer_Qualified_20261026.zip
   Library libfile_0826dd7b6a1881919d7c70c29410b64c
   Restored archive sha256a8e15017824c11b03a9118a27eb11c13d6480695bc075d454635654ea4e297f2
   Supplies fixed source algorithm, metadata and knowledge; raw technical halves are external to that compact archive.
3. WitnessCell_RNA_Clean_Manual_20261020.zip
   Library libfile_39757f9138b481918cbd85df7a1bd797
   JOUNG_STRICT_QC_CONDITION_MOMENTS.npz sha256cd677a71ff7a58fc9840fc7c2535b9c4582da79fe1db8f1f664a917910262c3d
   Previous package audit/factory/extract_train_halves.py reconstructs TRAINING_RAW_HALVES.npz.

run_five.py contains exact source paths and PROTOCOL.json records all input hashes.
To rerun full prediction, restore the sibling workspace source assets to those paths and use a fresh run directory in a separate working copy. Keep the archived run untouched. Stage order: freeze, smoke, predict, then score_five.py. The Q4 smoke also requires the previous full working prediction archive; it is an archived reproduction test, not training input.
source/legacy_model.py, source/fold_factory.py and source/case_features.py are unmodified copies; run_five.py calls only their required routines.

Evidence status
Historical-panel development extension. Q4 known result; other four pre-defined eligible tasks newly scored in this run. Five tasks share8edges and3 identity contrast directions, not five independent experiments. Case_R raw-half reliability is not calibrated to the7-library target. No new tuning or endpoint repair was performed.

Package integrity
PACKAGE_MANIFEST.json hashes every included file except itself. It is an archival checksum list, not preregistration. The original prediction freeze contains some omitted large-file hashes and is explicitly separate from this compact manifest.
