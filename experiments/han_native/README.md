# Han native refit v1 — completed experiment

Start with REPORT_CN.md and results/KEY_RESULTS.csv. All 60 prespecified batches completed; the manuscript itself has not been overwritten.

## Replay saved fitted parameters (no fitting of predictor weights)

Python 3.12; install requirements.txt. From this directory run:

```bash
python verification/compact_replay.py
```

This reconstructs selected inner and final corrections from saved parameters, refits no predictor weights, recomputes the scalar control from allowed inner predictions, and verifies all 1,080 reference score rows. The expected precision is given in replay_results/REPLAY_CHECKS.json. Numerical differences in another BLAS environment are assessed in original squared-error units (absolute tolerance 1e-8); bitwise equality of every intermediate is not asserted across different environments.

## Refit the frozen scientific procedure

Work in a separate copy to preserve the delivered records:

```bash
python code/engineering_checks.py
python code/prepare.py
python code/launch.py --workers 6
python code/score.py
python verification/independent.py
```

Preparation regenerates outcome-free query inputs and training-only files, plus a new execution freeze record. Workers checkpoint every candidate and completed outer batch. Interrupted workers resume completed artifacts; scientific parameters do not change. The scorer refuses to run until all 60 prediction locks exist. Candidate convergence failures remain distinct from selecting zero correction.

## Contents and lossless reduction

- inputs/: unchanged Han panel inputs, outcomes, and original physical-pair partitions. This reproduction starts from the archived v37 phenotype matrices; it does not reprocess FASTQ. The upstream count/readout reconstruction remains in grn_repo_v28_paper_v36.zip and the v37 analysis supplement.
- code/: frozen fitting, selection and scoring scripts; unchanged native scientific engines under vendor/.
- fitted/: final and selected inner model parameters, exact selected inner predictions, baseline fits and control coefficients. Symmetric spectral matrices store the upper triangle losslessly; rank-one matrices store factors and strength, checked to reconstruct the original matrix exactly.
- results/: all comparisons, per-repeat/per-partition results, all candidate inner selection scores and convergence diagnostics.
- verification/: original-engine equality and independent residualization checks; compact model replay.
- engineering_notes/: the rejected pre-scoring cache-collision run's incident and original freeze. Invalid model outputs are not used as evidence.

The redundant unselected inner matrices and expanded copies of six algebraically derived query predictions are omitted. Full refitting recreates all candidates; saved selected weights recreate every reported outer score. Original LOCK.json files identify the full run before compaction, so some hashed expanded files are intentionally absent here. MANIFEST.json describes this deliverable; replay_results verifies its actual sufficient contents.

## Statistical and usage scope

The primary analysis is reconstructed readout, native baseline, rank-one model. Spectral models, saturation baselines, published readout and component control retain their prespecified secondary roles. Gains are relative to the named baseline. Repeated partitions do not count as independent studies. This is a frozen-procedure reanalysis of previously accessed Han data, not an untouched prospective validation. Published aggregate outcomes are one readout, not biological duplicates.

Scientific-source rights follow the author-provided WitnessCell code archive and its rights statement; this experiment does not assign a new blanket licence. NumPy, SciPy, pandas and optional plotting dependencies retain their own licences.
