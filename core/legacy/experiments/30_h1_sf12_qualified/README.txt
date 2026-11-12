WitnessCell H1 SF12: five qualified rectangle tasks (paper v22, Fig. 7 / SI Note 12)

【自包含复算 / Compact replay】
  python run.py --out /path/to/fresh/output
  Equivalent: python run_replay.py --output /path/to/fresh/output
  Python 3.10+ and NumPy. CPU, no GPU, no network, no new training.
  Output: REPLAY_REPORT.json and ALL_TASK_SCORES.csv.

This module retains all five predefined tasks, all three models (M0, old12/SF12,
case_R), every failure, fixed protocols and original source code. It adds a
portable replay entry and an optional external-input reconstruction entry.

Expected headline values (gain = 100 * [1 - model/reference SSE]):
- SF12 complete RNA: edge-balanced gain 2.318114597% vs M0; positive on all 5 tasks.
- SF12 partner contrast: pooled gain 0.855498601%; positive on 3/5 tasks.
- SF12 excluding Q4: partner gain 0.384154477%.
- SF12 excluding Q5: partner gain -0.428933921%.
- case_R complete RNA: edge-balanced gain -94.539049921% (harmful).
- case_R partner: pooled gain 0.074367952%.
- 5 overlapping rectangles use 8 unique pair edges and 3 identity-contrast directions.
  They are NOT 5 independent biological replications.

Scientific interpretation / 结论边界
完整 RNA 的五块同向改善值得保留；伙伴收益仍不稳定，不能写成稳定的新伙伴优势。
case_R 是应保留的失败对照，不进入推荐模型。本实验属于历史面板的开发扩展，
不是新的独立来源确认。抽样/组合筛选规则来自计数和支持元数据，全部合格任务均保留。

What is exactly reproduced without external data?
1. Actual measured RNA vectors for 8 unique targets and all frozen partner
   prediction vectors directly reproduce every partner SSE and gain.
2. Archived per-edge full-response SSE values reproduce complete-response
   aggregates. Full RNA prediction arrays are omitted to keep the module small.
3. Added exact SF12_SHARED_TRAIN_STATS.npz independently re-solves all five
   12-coefficient SF12 fits. These statistics were derived from original verified
   intermediate inputs at frozen settings. All 60 coefficients matched the
   original TRACE exactly; rebuilt partner corrections matched within 3.82e-16.
4. All saved inner-episode normal equations also re-solve the one case_R
   coefficient per task; quadratic statistics check its full/partner losses.
5. Task splits, overlapping contrasts, rank=3, all aggregates and leave-one-out
   sensitivity are verified. Compact replay does not reconstruct RNA features
   from cells or re-estimate preprocessing. SF12 parameter solves are exact
   sufficient-statistic solves, not raw-cell training.

【完整中间数据重建 / Optional reconstruction from external intermediate inputs】
  python reconstruct_from_inputs.py \
    --asset /path/TRAINING_BATCH_ONLY.npz \
    --halves /path/TRAINING_RAW_HALVES.npz \
    --output /path/to/fresh/reconstruction
  Add --check-inputs-only to check identity/hashes without fitting.

The runner verifies input SHA-256, uses packaged KNOWLEDGE_NODE.npz, preserves
all historical fixed settings, refits calibration/M0/SF12 per support graph,
checks Q4 archived contrast reproduction and withheld-count poisoning, freezes
all five predictions, then scores. No hyperparameter search is performed.
Original scripts run_five.py and score_five.py are unmodified historical sources;
use the wrapper for portable paths, not their historical default paths.

The added statistics can themselves be reconstructed without technical halves:
  python derive_sf12_statistics.py --asset /path/TRAINING_BATCH_ONLY.npz \
    --out /path/to/fresh/sf12-statistics
This bounded command was executed in 2.5 seconds and reproduced all archived
SF12 coefficients exactly. It changes no model settings and reads only each
allowed support graph when constructing features.

External inputs (not bundled; provenance/RECONSTRUCTION_INPUTS.json has exact hashes):
- TRAINING_BATCH_ONLY.npz, 32,146,451 bytes:
  WitnessCell_Neural_Partner_Training_Inputs_20261021.zip,
  Library libfile_ee62d07de8c881919088359f4ae3b89c.
- TRAINING_RAW_HALVES.npz, 6,905,980 bytes:
  derived from JOUNG_STRICT_QC_CONDITION_MOMENTS.npz in
  WitnessCell_RNA_Clean_Manual_20261020.zip,
  Library libfile_39757f9138b481918cbd85df7a1bd797.
  Original extraction code is in WitnessCell_Case_Transfer_Qualified_20261026.zip.
This is intermediate-data reconstruction, not extraction from raw cells.

Code and data map
- source/legacy_model.py::old12: actual 12 shared-feature SF12 model and ridge solve.
- source/fold_factory.py: fold-local calibration, preserved M0 endpoint shrinkage,
  GenePT axes and training-side technical-repeat reliability.
- source/case_features.py: empirical-case comparator features.
- data/KNOWLEDGE_NODE.npz: exact 10 x 64 frozen GenePT node inputs.
- provenance/GenePT_provenance/: official GenePT v2 download and fixed projection
  provenance. FROZEN_RUN.json is a historical upstream record; its neural-training
  settings do not describe SF12. SF12 uses frozen external embeddings, no encoder training.
- run/: original protocol, scores, compact vectors, training statistics,
  per-task traces and original prediction-freeze hashes, plus explicitly named
  SF12_SHARED_TRAIN_STATS.npz added during integration. The new statistics were
  not retroactively inserted into the historical prediction freeze.
- audit/: original historical independent full-vector audit; it requires omitted
  original full vectors and inputs. Use run_replay.py for self-contained verification.
- replay_compact.py: unchanged original compact replay.
- provenance/ORIGINAL_PACKAGE_MANIFEST.json: original archive member hashes;
  provenance/ORIGINAL_README.txt: original scope and asset recovery notes.
- MODULE_MANIFEST.json: current module member hashes, excluding itself.

No scientific numerical results or fixed model settings were changed by integration.
The SF12 normal-equation replay is a newly added reproducibility layer, not a new
experiment or additional validation task.
