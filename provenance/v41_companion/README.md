# WitnessCell v41 revision supplement — language revision R1

This companion to `grn_repo_v30_paper_v39.zip` supplies the v41 manuscript and figures, with the unchanged v40 numerical analyses. It supersedes the v40 companion for this paper version; the base archive is still required. v41 is an editorial correction, with no new model fitting or scientific results. Language revision R1 updates the manuscript and captions only; LANGUAGE_AUDIT_CN.md records the prose changes.

## Run

Unpack both archives. The base argument is the extracted `grn_repo_v30` directory. Use Python 3.12 and `requirements.txt`. Choose a new output directory outside this supplement.

```bash
python reproduce.py --base /path/to/grn_repo_v30 --out /path/to/v41_replay
python reproduce.py --base /path/to/grn_repo_v30 --out /path/to/v41_refit --refit
```

Default replay uses saved four-configuration fitted predictions and recomputes the nine analysis tables, Han attribution and figures. `--refit` also fits the saturation reference, rank-one and spectral learners afresh for the three K562 high-support folds and Jurkat. The fixed model grids and query identities are unchanged. Both routes automatically run the base archive's deterministic alias/cache restoration; this writes restored files into the extracted base archive and involves no predictor fitting. Keep output paths separate.

Add `--render` to rebuild the two manuscript PDFs. Rendering additionally requires pandoc, XeLaTeX, Liberation Serif/Sans and DejaVu Sans Mono. Individual and combined figure PDFs are generated without `--render`. The saved manuscript PDF can also be read in the light package.

## Contents and exact scope

- `ANALYSIS_PLAN.json`: original v40 plan recorded before those fits, input hashes and stopping rule; retrospective analysis of previously accessed data.
- `inputs/`: compact original permitted data for four configurations and the v39 common-reference source table.
- `fitted/`: selected models, inner predictions, locked query predictions and coefficients.
- `code/refit_saturation.py`: same Han saturation-reference family, with unchanged rank-one/spectral fitting and zero options.
- `code/analyze_v40.py`: common-reference comparison, saturation scoring and whole-correction Han decomposition.
- `code/validate_v40.py`: independent direct-design least-squares check, comparison with original Han gains and the nine expected tables.
- `results/`: nine numerical analysis tables and preserved v39 figure-source tables.
- `paper_sources/`, `render/`: current Markdown, vectors and rendering tools.
- `validation/`: v40 refit records under `historical_v40/` and v41 saved-prediction replay/editorial checks. The two execution scopes are separate.

The nine analysis tables are COMMON_REFERENCE, SATURATION_SCORES, SATURATION_POOLED, SATURATION_BY_CONFIGURATION, SATURATION_BY_REPEAT, HAN_DIRECTION_ATTRIBUTION, HAN_ATTRIBUTION_POOLED, HAN_ATTRIBUTION_BY_PARTITION and HAN_ATTRIBUTION_BY_REPEAT (CSV).

Han attribution uses all 60 saved batches and both learner families from the base archive: 120 fitted outputs and 360 repeat/template records. It neither fits Han predictors nor chooses coefficients using query outcomes. The decomposition of observed error gains itself uses query outcomes and is descriptive. All signed contributions are retained.

K562 saturation rank-one partner gain is 11.3241%. Jurkat selects zero rank-one and yields 1.9872% spectral gain over saturation. Han saturation spectral gain orthogonal to the specified single-product template is 4.7865 and 7.0300 percentage points in the reconstructed and published readouts. These statements use their own baseline denominators and do not exclude every possible global response model.

## Validation and tolerances

The v40 clean `--refit` run reproduced 44 numerical arrays and nine result tables with maximum absolute difference zero in the recorded environment. Sixteen independent direct-design baseline solves differed by at most 1.5543e-14. Floating-point comparison in the supplied independent checker uses explicit absolute/relative tolerances; byte identity of ZIP/PDF files is not the model reproducibility criterion. MANIFEST.json provides file hashes.

The v41 delivery executes saved-prediction replay and document checks; it does not repeat that v40 refit. Scientific filenames `analyze_v40.py` and `validate_v40.py` identify unchanged analysis code, while rendering uses v41 sources.

The historical benchmark, RNA models and Han predictors are not retrained by this supplement. Their replay/refitting routes remain in the base archive. `RIGHTS_AND_RELEASE.md` states the existing licence scope; no new blanket software grant is implied.
