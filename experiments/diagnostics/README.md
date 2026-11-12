# WitnessCell v37 analysis supplement

This companion contains the new retrospective diagnostics and the two targeted Han attribution controls. It is separate from the light manuscript. Existing biological results and source weights are unchanged. The analysis fits target comparator baselines and scalar amplitudes; it does not retrain the source SF12, V4 or rank-one predictors.

## Reproduce the new analyses

Use Python 3.12 and the versions in ENVIRONMENT.json. From this directory:

```bash
python -m pip install -r requirements.txt
python code/han_attribution.py --han inputs/han --out replay
python code/functional_diagnostics.py --control inputs/control --out replay
python code/summarize.py --results replay
python code/check_analysis.py
python code/make_figures.py --results replay --out figure_replay
```

The default run uses 999 fixed random directions per source, all 30 outer folds across the two Han readouts, and the original calibration rule. No download or Torch installation is required for this prepared-feature replay. `replay/` is a new output directory; preserve `results/` as the supplied reference. `check_analysis.py` compares replay with references and checks query-outcome exclusion, independent least squares, SVD projection and the SSE identity. Its new validation report is written under `replay/`.

`inputs/han/FROZEN_SF12_FEATURES.npz` stores exact source-derived features, not fitted target labels. Replaying them preserves the complete fixed source direction and all original target splits. The native source-weight route was also run in this revision. To repeat that route, extract `grn_repo_v28_paper_v36.zip` and run the same script with `--han /absolute/path/grn_repo_v28/experiments/han2017/snapshot --out native_replay`, using the base archive's PyTorch dependencies. Its pretrained inputs remain frozen. The two routes produced identical reported numerical tables in the recorded environment.

## Sources and outputs

- ANALYSIS_PLAN.md: fixed retrospective plan written before these new calculations; Han outcomes had already been accessed.
- inputs/control: exact archived arrays and internal scalar decisions for 55 configurations; the ten primary functional fits also have repeat diagnostics.
- inputs/han: original full-screen targets, singles, predictions, pair splits and scalar parameters, plus frozen SF12 features and source coefficient information.
- results: fold scores, all 999-direction summaries, nonlinear comparator parameters, exact oracles, repeat correlations and validation records. The native RUN_RECORD preserves its original script hash; the current wrapper additionally supports prepared input. PREPARED_RUN_RECORD records that replay and its hash. Neither record is retroactively rewritten.
- inputs/figures: compact, unchanged original sources for retained Fig.4 panels; code/make_figures.py regenerates updated Fig.4 and new supplementary Figs.6–7. It does not replace the unchanged Fig.1–3 or Fig.5.
- CLAIM_SOURCE_MAP_V37.csv: new claim/table/figure to source mapping.

Full-precision input and output bytes are verified by MANIFEST.json. Numerical replay uses explicit tolerances, not a promise of cross-environment byte identity. Random tail fractions describe the selected coefficient reference; repeat correlations are not theoretical noise ceilings; outcome-informed oracles are never deployed.

## Render the manuscript

With Pandoc, XeLaTeX, Liberation fonts and DejaVu Sans Mono installed, point the renderer at the extracted light paper directory:

```bash
python render/render_documents.py --paper /absolute/path/WitnessCell_v37_Light
```

The output PDFs replace the existing PDFs there. Render logs are local build outputs. Publication figure data and scripts stay in this analysis supplement.

## Base release and permissions

Existing theory, benchmark, RNA, functional and application results are reproduced from `grn_repo_v28_paper_v36.zip` (SHA256: 0cf6df0ae3db1b1ba0bdedbdf7723c42678ba4e6bedfb8254b7ff861071a2adb). This is an increment, not a replacement for that complete release. The base file remains available with the previous manuscript delivery.

This archive does not invent a new blanket software licence. The base archive's five-file MIT scope and third-party notices remain intact. Author authorization for additional scientific code and redistribution must be resolved before public release; providing a review copy does not itself supply that authorization.
