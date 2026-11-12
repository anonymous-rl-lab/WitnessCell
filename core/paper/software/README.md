# WitnessCell: new-input prediction example

This delivery provides runnable interfaces to the **released transcriptomic core** and the **functional residual model**. The transcriptomic adapter uses the original `multihead_v14` identity head. Seven source modules are copied without modification from `grn_repo_v27/compact/legacy_core/source/v14/`; their hashes are recorded in `INTERFACE_TESTS.json`. It predicts new query identities from user-supplied condition summaries rather than replaying stored predictions.

## Run

Python 3.11 or later, NumPy and SciPy are required. From the extracted package directory:

```bash
python -m pip install -r software/requirements.txt
python software/witnesscell_predict.py --input software/example_input.json --output my_predictions
python software/test_interface.py
```

The supplied example has 24 readouts, seven observed training singles, 18 training doubles, five validation doubles and five query doubles. G7 is a query endpoint without an observed training single. Its annotation groups and measurements are synthetic and demonstrate the interface only. They are not added to any biological result. `make_example.py` regenerates this input.

## Prepare your own input

Copy the JSON schema in `example_input.json`. `genes` is the ordered readout axis. All `mean` and sample `variance` vectors must use exactly that order and the same normalization. `n` is the cell count for the condition. Provide control, permitted training observations, validation observations, and query **identities only**. The core does not normalize raw counts. Apply the preprocessing appropriate to your assay before computing these summaries; do not jointly preprocess using query outcomes for model selection. `gene2go` maps gene symbols to lists of GO annotation identifiers from a documented annotation release.

Single identities are gene symbols; doubles use `GENE1+GENE2`. Query identities may be unseen singles or doubles. The adapter requires at least one training single and training/validation doubles for its supported fitting path. Every endpoint must occur on the readout gene axis. Specify the condition partition before fitting. Canonicalized duplicate pairs, reversed-pair leakage, query measurements, mismatched vectors and negative variances are rejected. This interface assumes one consistently defined intervention mode for each gene; mixed activation/repression inputs require the separately evaluated mode-aware implementation.

`predictions.csv` contains the predicted normalized condition mean and its nonlinear baseline for every query/readout. `fit_metadata.json` identifies the selected saturation, residual regularization and amplitude, observed single identities, fitting partitions and input hash. Differences between complete and baseline predictions describe endpoint refinement, not a new partner contrast. The outputs do not include calibrated uncertainty.

## Choose a configuration before examining query outcomes

| Input and scientific objective | Configuration | Output and scope |
| --- | --- | --- |
| Normalized RNA summaries, controls, measured or missing singles, training/validation doubles; complete-response prediction | Released core, demonstrated here | Response means; Fig. 3 benchmark belongs to this configuration |
| Matched RNA screen with measured endpoint responses, training cell halves, enough pair overlap to estimate residual structure | R partner extension | Additional readout-weighted partner correction; use the recorded experiment implementation in the full archive; no claim of universal transfer |
| A fixed anchor and a specified transcriptional objective, with mode-aware candidate responses | PairPick with R-G | Utility and partner choice; do not substitute the released-core demo for R-G results |
| Repeat-specific functional measurements and single/double support | Functional residual model | Predicted growth and partner-sensitive residuals; refit in the new context |

Choose by input type and intended output, then use training/validation data for fitting. There is no validated automatic selector among these implementations. Inadequate partner-support geometry gives an undefined or uninformative partner diagnostic; it does not authorize a positive partner claim. A conservative user can retain the complete-response configuration and report its scope.

## Validation and reproduction scope

`test_interface.py` checks exact parity against a direct call to the unchanged native core, finite output, unseen-single handling, absence of query truths, rejection of leakage and malformed inputs, and equivariance to a common shift of all observed responses. The stored predictions are an example output; the command recomputes them. This verifies the interface, not performance on a new biological source. Original-array score replay, neural-weight replay and full retraining remain distinct operations described in Supplementary Note 12.

## Functional screen completion: real processed input

```bash
python software/witnesscell_growth.py --input software/growth_example_input.json --output my_growth_predictions
python software/test_growth_interface.py
```

The input contains the original Jurkat training measurements: 79 endpoints, 2,464 training pairs and two source repeat channels. The 617 query records contain identities only. The fitter reconstructs the native three training folds, chooses penalties using training outcomes, and refits the repeat-specific baseline and shared rank-one/spectral corrections. It reads neither query outcomes nor the stored expected predictions. The supplied example reproduces all three archived prediction arrays exactly. This is a new input interface to an existing evaluation, not a new biological experiment.

For new data, use the schema in `growth_example_input.json`: unique endpoint names in `genes`; finite single-perturbation responses of shape endpoints × 2 in `singles`; training records with `pair` and two repeat values in `response`; and identity-only `queries`. Use the same control-centred growth units throughout. Define the training/query split before fitting. The adapter requires two repeat channels and enough training identities for three nonempty native inner folds. It rejects self-pairs, unknown endpoints, repeated/reversed pairs, training–query overlap, nonfinite measurements and query outcomes. It does not ingest raw sequencing counts or claim cross-assay preprocessing support.

The output CSV contains one row per pair and repeat: `pair`, `repeat`, `baseline`, `rank1_prediction`, `rank1_increment` and `spectral_prediction`. The fitting output also records selection and provenance. Select rank one as the primary configuration before looking at query outcomes; the spectral extension is supplementary. The adapter always applies the archived high-support fitting procedure, refitting parameters to the supplied context. The archived low-support condition is an evaluation contrast, not an automatically selected operating mode. `growth_example_expected.npz` is used only by the parity test.

## Complete scientific source

The original **unmodified** `grn_repo_v27.zip` is included at `reproduction/grn_repo_v27.zip`. Its internal README and replay-scope documents distinguish score reconstruction, checkpoint replay and full retraining. Extract it separately and follow the asset-specific instructions; the archive is not claimed to contain every external raw sequencing file or every large neural checkpoint. The current paper adds no new neural training. See `reproduction/README.md` and `RIGHTS_AND_RELEASE.md`.

## Licensing and public release

`LICENSE-ADAPTER.txt` grants MIT terms for the new transcriptomic and functional adapters, generator and interface tests. The seven transcriptomic and three functional scientific source modules retain their provenance and existing rights; no upstream licence was present in the supplied archive. A unified licence approved by the authors/rightsholders and a persistent public repository must therefore be supplied before public software release. Third-party dependencies retain their own licences. This reviewer delivery does not claim that public release has already occurred.

## Fresh v28 execution check

From the package root, run `python code/check_interfaces.py`. This fits both examples in a new temporary directory and runs the interface suites against those newly generated outputs. It preserves the distributed example outputs. Detailed logs and measured runtimes are saved under `revision/interface_checks/`. The measured environment is recorded in `FRESH_INTERFACE_RUNS.json`; installation time was not measured because dependencies were preinstalled. These examples run on CPU; GPU hardware is not required. Full neural-model retraining is a separate workflow.

The functional example measures 2,464 / 3,081 pairs (approximately 80%) and predicts 617 / 3,081 pairs (approximately 20%). This is the evaluated setting, not an empirically established minimum support threshold. Define an appropriate validation split for a new screen before examining query outcomes.

## Component-specific adjustment

After fitting a response model, `../partner_control/software/witnesscell_partner_control.py` fits a scalar from supplied internal held-combination predictions and adjusts the partner part on a specified query panel. See `../partner_control/README.md` for the runnable real-data example and strict input schema. It does not replace the biological-model fitting commands above.
