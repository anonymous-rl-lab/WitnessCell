# WitnessCell: paper reproduction

Code release **31.0.0**, for **WitnessCell v41, language revision R1**.

This repository integrates the complete base archive and the v41 companion. No separate historical ZIP is required for the default paper replay. It covers the released response predictor, target-information theory, RNA and functional partner prediction, component control, native Han analysis and matched global-reference analyses.

## Quick start

Tested on Linux with Python 3.12.14. Use a fresh environment and install the pinned dependencies:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -B tools/verify.py
python -B reproduce.py --out ../witnesscell-replay --jobs 2
python -B tools/verify.py --restored
```

The default run has 16 tasks and does not retrain biological predictors. It reconstructs selected predictions from saved weights, recomputes scores, ranks, intervals, calibration coefficients and diagnostics, regenerates the 12 figures, and checks current manuscript/table identities. Each task writes a log, a status and its outputs under the selected output directory. Successful task status means the named task's assertions passed, not that raw sequencing and all original training were rerun.

Start with a short subset or continue an interrupted run:

```bash
python -B reproduce.py --list
python -B reproduce.py --out ../witnesscell-smoke --tasks paper_scores theory_growth interfaces coverage
python -B reproduce.py --out ../witnesscell-replay --jobs 2 --resume
```

Outputs must be outside the repository. Resume checks the release/task fingerprint before skipping a completed task. Failed runs retain their previous outputs. The preparation step restores byte-identical aliases and regenerates 120 Han prediction caches from saved weights; it does not fit predictors or read query outcomes.

## Reproduction levels

| Command/task | What is executed |
|---|---|
| Default `reproduce.py` | Saved-array/weight replay, specified scalar/reference fits, scoring and figure regeneration |
| `--tasks global_refit` | Seven original matched-product-reference refits with fixed folds and grids |
| `--tasks saturation_refit` | Four saturation-reference refits with fixed folds and grids |
| `--tasks paper_pdf` | Current Markdown to main/Methods and SI PDFs; requires Pandoc, XeLaTeX and fonts |
| `--tasks han_frozen_inference_optional` | Historical SF12/V4 Torch inference; separate optional environment |

The native Han full fitting scripts and frozen plans remain in `experiments/han_native/`; the default native-Han task uses its saved weights. See [reproduction levels](docs/REPRODUCTION_LEVELS.md) for exact scope and tolerances. Large raw FASTQ and some historical neural-baseline weights are represented by original source/acquisition indexes; the default replay does not require them.

## Repository map

| Path | Contents |
|---|---|
| `reproduce.py`, `REPRODUCTION_INDEX.json` | One current task entry and explicit scope |
| `paper/manuscript_sources/` | Authoritative v41 R1 Markdown, legends and vectors |
| `paper/reference_figures/` | Current figure comparison references |
| `experiments/reference_sensitivity/` | Integrated v41 companion: common denominators, saturation references and Han direction decomposition |
| `experiments/global_reference/` | Seven matched-product fits, mode diagnostics and theoretical checks |
| `experiments/han_native/` | Fixed inputs, splits, weights, native/saturation analyses and fitting code |
| `experiments/diagnostics/` | Repeatability, fixed-direction reference optima and historical attribution checks |
| `core/` | Established theory, benchmark, RNA, growth, control, interfaces and their source assets |
| `docs/PAPER_COVERAGE.csv` | Every SI table and figure mapped to its evidence and task |
| `provenance/` | Source hashes, integration map, pruning and lossless deduplication records |
| `verification/` | Actual checks for this code release |

Some historical paths retain version strings because frozen programs or provenance refer to them. They are supporting dependencies, not alternative current entry points. Current instructions are this README and the task index.

## Numerical fidelity

Input arrays, weights, full-precision scores, folds, primary/secondary designations and scientific algorithms are preserved. Exact file checks are separate from numeric tolerance checks. `--compatible` permits deterministic cache bytes to differ across numerical environments; it does **not** relax downstream scoring tolerances. Use the pinned environment for the strict route. See [tolerances](docs/NUMERICAL_FIDELITY.md).

## GitHub use and licensing

The package contains a CI workflow, smoke/regression tests, pinned dependencies, an environment diagnostic, a citation file and a minimal container recipe. No repository has been published and no hosted CI run is claimed. Large read-only inputs are losslessly packed and restored locally; ordinary output/cache paths are ignored by Git.

Licence boundaries are inherited. The existing MIT notice covers only its five named interface/test files. No repository-wide open-source licence has been granted for the scientific modules; rights-holder confirmation remains necessary before public distribution. See [licence scope](LICENSE.md) and [third-party notices](THIRD_PARTY_NOTICES.md).

The paper is a manuscript; no publication DOI or public repository URL is invented in `CITATION.cff`.
