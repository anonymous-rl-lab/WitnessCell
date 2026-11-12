# Reproduction levels and resources

| Level | Entry | Data and compute scope |
|---|---|---|
| Integrity | `tools/verify.py` | Exact shipped-file checks; no fitting |
| Default paper replay | `reproduce.py` | Sixteen tasks, saved scores/arrays/weights, new scalar calibration and designated small reference fits, source checks and all 12 figures |
| Matched-product refit | `--tasks global_refit` | Seven configurations; frozen learner/selection rules |
| Saturation-reference refit | `--tasks saturation_refit` | Four configurations; frozen learner/selection rules |
| Native Han training | `experiments/han_native/code/fit_task.py` and its retained task plan | Full per-task baseline and matrix fits; opt-in, potentially lengthy |
| Historical Torch inference | `--tasks han_frozen_inference_optional` | Optional dependency environment; source predictors remain fixed |
| Manuscript rendering | `--tasks paper_pdf` | Pandoc, XeLaTeX, TeX packages and fonts listed by paper/render and paper/requirements.txt |
| Full raw-data regeneration | Original acquisition/processing indexes under core | May require external FASTQ and large neural assets; not included in the default claim |

Use Python 3.12 and the exact versions in requirements-tested.txt. The recorded integration environment is in verification/ENVIRONMENT.json. Default replay uses no GPU and requires no Torch. Installation can need Internet access; default scientific replay uses local data after installation. Keep at least 2 GB of free disk space for restored views and outputs; peak memory and runtime depend on task concurrency and BLAS implementation.

For an interrupted run use the same output with `--resume`. A completed task is skipped only when its fingerprint matches the current release/task specification. Unknown tasks, repeated task names, and unsafe output locations are rejected.
