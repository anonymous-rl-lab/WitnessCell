# v1.1 formal SOTA audit release

- Completed 12 CPU fits across four official genetic-combination datasets and
  three split seeds.
- Added exact official six-metric scoring, resumable seed checkpoints, and an
  algebraically equivalent multithreaded E-distance implementation.
- Reproduced every published aggregate rank exactly before inserting
  WitnessCell.
- Added the official top-100 leaderboard, all-gene robustness leaderboard,
  dataset-equal sensitivity analysis, perturbation-cluster bootstrap, audit
  report, and publication figure.
- Corrected the comparison scope: 27 methods exist benchmark-wide, but 15
  published methods have eligible genetic-combination results. WitnessCell is
  the 16th entrant in that track.
- Bundled the two verified read-only GEARS GO-support assets needed to recreate
  the published splits. No GEARS, CUDA, PyTorch, or CPA installation is used.
- Fixed the official-repository adapter so the GEARS asset path is passed
  explicitly.
- Formal scoring now consumes target-free `deploy_predictions.npz` archives.

Primary result: WitnessCell position 1, mean per-operation rank 5.1437, under
the official top-100 six-metric aggregation over 654 operation-seed units.

## v1.2 evidence organization update

The formal evidence is now packaged as two standalone experiment units:

- `01_four_dataset_cpu_prediction` freezes the 12 CPU predictions, split
  identity, target isolation, selected settings, and deployment audit.
- `02_official_full_metric_sota` freezes all six official metrics, exact
  published-rank reconstruction, formal SOTA verdict, robustness panels, and
  perturbation-cluster uncertainty.

Both units have independent one-command runners and audit scripts. No
scientific number or model output was changed in this organizational update.
