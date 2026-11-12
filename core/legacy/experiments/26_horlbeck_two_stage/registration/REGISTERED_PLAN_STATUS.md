# Registered plan status and relationship to the two-stage results

Prepared on 18 October 2026 (UTC), after the analyses described below. This is a retrospective status report accompanying the results deposit; it is not a new preregistration and does not change the original registration.

**Research question retained:** Can measured combinations improve prediction of unmeasured gene pairs beyond single-perturbation information, endpoint correction and a training-fitted nonlinear response baseline, and how does that benefit depend on the available combination measurements?

The research programme began with the OSF-registered K562 analysis in `OSF_PREREGISTRATION.txt`, version 1.1, distributed with `WitnessCell_Horlbeck_Freeze_v1.zip`. Subsequent method development and the held-environment Jurkat analysis address this same question with revised methods and budgets. They are reported separately from the original registered tests.

## 1. Original registered analysis

The original analysis used 80 K562 genes and 3,160 eligible unordered gene pairs. One fixed test set contained 632 pairs. Nested training budgets were 158, 316, 632, 1,264 and 2,528 pairs (5%, 10%, 20%, 40% and 80% of all eligible pairs). **The 20% training budget was the sole primary budget.** All guide combinations, orientations and both experimental replicates of a pair shared its partition.

Each replicate was fitted separately. Training-only three-fold validation selected among ranks 0–3 and the registered penalty grid. The principal comparison was the selected procedure against its corresponding rank-zero baseline. The registered contrast analysis used up to 1,000 fully held-out rectangles.

| Registered objective | Execution and result | Status |
|---|---|---|
| H1: lower held-out growth-prediction error at the 20% training budget | Both replicates selected rank zero. Selected and baseline predictions were identical: gain 0; conditional 95% interval [0, 0]; both replicate gains 0. | Positive-evidence criterion not met. Substantially limited by failure of positive-rank candidate optimization. |
| H2: lower error for held-out rectangle contrasts at the same budget | Gain 0; conditional 95% interval [0, 0]; both replicate gains 0. | Positive-evidence criterion not met, with the same computational limitation. |
| H3: cross-replicate agreement after separately fitted baselines | On 1,000 fixed overlapping test rectangles, adjusted Pearson correlation 0.276160 and Spearman correlation 0.280207. Raw correlations were 0.347450 and 0.338632. | Completed as a descriptive analysis. Positive correlation does not establish low rank, extrapolation, or rescue H1/H2. |
| H4: gains and support across five measurement budgets | All five budgets and both replicates completed; the selected procedure had zero incremental gain at every budget. Baseline MSE improved with increasing training information. | Completed as a descriptive analysis, but positive-rank numerical failures limit inference about the absence of usable partner information. |
| Cross-replicate prediction and descriptive controls | Saved predictions were evaluated in the other replicate without destination refitting; selected-versus-baseline gains were zero. Additive, endpoint and shared-quadratic controls were scored. | Completed, secondary/descriptive. |
| Original progression criterion for Jurkat | The original criterion required positive H1 and H2 evidence. It was not satisfied. | Later Jurkat work followed a separately frozen, revised workflow motivated by subsequent K562 development; it is not a successful passage through this original gate. |

At the primary budget, 266 of 270 positive-rank inner-fold fits failed all three prescribed starts. Neither replicate had a positive-rank candidate with three eligible validation folds. Of 810 initialization runs, 806 exhausted the 2,000-iteration limit. Final rank-zero selection therefore occurred during inner selection; it was not a positive-rank outer-refit fallback. The exact zero gains and intervals describe identical predictors, **not a biological effect estimated to be exactly zero**.

All 10 budget-by-replicate jobs completed in approximately 39 minutes of accumulated fitting time, below the six-hour cap. One interrupted session was resumed from the same checkpoints without resetting the budget or changing the registered rules. The original test outcomes were opened on **2026-10-18 at 12:43:30 UTC**, and scoring completed at **12:43:31 UTC**. Original outputs were retained.

For completeness, the original selected-model MSEs (identical to M0 at each budget) were:

| Training pairs | Budget | MSE | Incremental gain |
|---:|---:|---:|---:|
| 158 | 5% | 0.005960586219 | 0% |
| 316 | 10% | 0.004503994795 | 0% |
| 632 | **20%, primary** | **0.003964133649** | **0%** |
| 1,264 | 40% | 0.003696804444 | 0% |
| 2,528 | 80% | 0.003560078858 | 0% |

These results remain the results of the original registered tests. Later positive findings do not replace them.

## 2. Subsequent K562 development

After the original results were inspected, numerical diagnosis, regularization and repeat-sharing experiments informed method revisions. The completed two-stage package contains the final K562 comparison: fixed-rule M0, a shared signed rank-one correction and a shared symmetric spectral-shrinkage correction. The penalties act on the correction matrix, and complete predictions are selected by training-only validation.

Three outer development folds evaluated **632 distinct pairs drawn from the original non-test pool**. This is a different set of 632 pairs from the original formal test set. Low-density training contained 421, 421 and 422 pairs; high-density training contained 2,317, 2,317 and 2,318 pairs. Each fold excluded its own evaluation pairs at both budgets. The original formal test pairs did not enter these later fits or scores.

| Development setting | Rank-one gain versus same-budget M0 | Spectral gain versus same-budget M0 |
|---|---:|---:|
| Low density | +0.2652% | −0.2879% |
| High density | +11.0292% | +11.3865% |
| High density, after eliminating free endpoint-additive effects | +11.6838% | +10.6710% |

These are development results: the first fold was inspected during the smoke test, and K562 had already informed method choices. High-density gains are not a rerun or a positive replacement of registered H1 at 20%. The later endpoint-eliminated projection metric also differs from the original fixed-rectangle H2 endpoint; its improvement is reported as complementary evidence, not a retrospective H2 success. The two-budget comparison establishes neither a minimum sample size nor an asymptotic sampling law.

## 3. Jurkat held-environment confirmation

The original registration reserved Jurkat for subsequent work. The final K562 model code and selection rules were frozen before Jurkat fitting and held-out scoring. Rank one was designated the primary method, with the spectral method secondary. Parameters were re-estimated from Jurkat training data using those fixed rules.

The first availability check found no usable single-perturbation measurement for GINS4 in either Jurkat replicate. A dated, pre-fit availability amendment removed that gene without replacing it, leaving 79 genes and 3,081 eligible pairs. The fixed split used 2,464 training pairs and 617 test pairs, with both replicates and all equivalent measurements grouped by pair.

Local records show the amended freeze at **15:32:20 UTC**, saved-prediction lock at **15:32:52 UTC**, and test scoring at **15:33:24 UTC**, all on 18 October 2026. These are execution records, not earlier OSF upload timestamps. The protocol was locally frozen before evaluation; this results deposit does not claim a separate OSF registration of that revised protocol before evaluation.

The primary rank-one method reduced held-out MSE by **6.4308%**, with **4.4542%** and **8.8845%** improvement in the two replicates. Its endpoint-eliminated error improved by **6.1750%**. The spectral method's corresponding improvements were **8.1606%**, **8.0730%**, **8.2692%** and **5.8975%**. The primary direction-replication condition was met. Its nominal 95% gene-composition reweighting interval was **1.7335%–12.1290%**; this is conditional sensitivity evidence, not a calibrated guarantee over new training experiments or biological batches.

This is confirmation of a frozen workflow in a held-out cell environment within the same original study. It is not an independent-study replication, transfer of K562-fitted parameter values, RNA-prediction validation, or newly collected prospective wet-laboratory evidence. High training coverage remains part of the demonstrated scope.

## 4. What has and has not been completed

- **Completed:** the original registered K562 execution and reporting, including null and numerically constrained results; the revised K562 three-fold development comparison; the locally frozen Jurkat held-environment analysis; reconstruction and full computational reproduction of the final two stages.
- **Not established:** positive original H1/H2 results at the registered 20% budget; high predictive gains at sparse budgets; exact biological rank one; a minimum measurement count; an independent-study or RNA validation; a 75% biological hit rate.
- **Not claimed:** that all K562 outcomes remained unseen. The original 632-pair formal test was scored before subsequent method development. The final K562 development evaluation used a separate pair set. Jurkat held-out scores were unavailable to fitting and model selection until the recorded scoring step; computational data aggregation and sealing preceded that step.

## 5. Files and evidence trail

The companion `ANALYSIS_AMENDMENT.md` describes the scientific and computational changes. The compact archive is:

`WitnessCell_Horlbeck_TwoStage_Repro_20261018.zip`

SHA-256: `b60e01230f3f43f156271697fe1fdc52b6c4ce5f68b522c3664c6e4aa6ff4082`

Inside it, `stages/K562/` and `stages/Jurkat/` contain the final model code, processed inputs, exact partitions, saved predictions, results, freeze/provenance records and reproduction entry points. The compact archive is a reproduction package for these two final stages, not an archive of every earlier development attempt or the full source GEO gzip.

The original execution is separately retained in `WitnessCell_Horlbeck_Registered_Results_20261018.zip` and its accompanying report. Relevant source records are `results/TEST_OUTCOMES_OPENED.json`, `results/results.json` and `CONVERGENCE_SUMMARY.json`; original specification sections 4, 8, 19–26 and Appendices A–C define the registered objectives and decisions. This status document summarizes those original outcomes so the current deposit can be interpreted without relabeling them.
