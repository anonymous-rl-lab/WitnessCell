# Analysis amendments and data-access chronology

Prepared on 18 October 2026 (UTC), after completion of the analyses below. This document accompanies the results deposit for the OSF-initiated WitnessCell Horlbeck study. It records changes and their timing; it does not retrospectively preregister those changes or replace the original registration.

## 1. Continuing research objective

The objective remained to determine whether measured combinations provide usable predictive information beyond single-perturbation effects, endpoint correction and a training-fitted nonlinear response baseline, and how that benefit depends on combination measurement coverage.

The study began with an OSF preregistration. Following numerical failures and development analyses in K562, the modeling workflow was revised, frozen before inspection of Jurkat evaluation outcomes, and assessed in that held cell environment. The original registered execution and subsequent analyses are reported separately. `REGISTERED_PLAN_STATUS.md` reports the original hypotheses, their outcomes and the later results.

## 2. What changed, and why

| Component | Original registered procedure | Revised two-stage procedure and rationale |
|---|---|---|
| Positive-rank fitting | Separately fitted replicate models; training-only selection among ranks 0–3 | Numerical diagnosis was followed by development of a shared correction. The original primary run had 266 failures among 270 positive-rank inner fits. Later converged independent-factor fits still did not improve development predictions, motivating changes to representation rather than simply longer optimization. |
| Baseline | Training-fitted endpoint/nonlinear baseline under the registered selection procedure | Fixed M0 method with alpha = 10, refitted separately for each replicate and each permitted training subset. Holding this rule fixed makes incremental correction gains interpretable. The same fitted predictions are used by competing correction methods within a budget. |
| Replicate information | Partner factors fitted separately | One symmetric correction matrix G is shared by the two replicates in original growth-response units. Replicate baselines remain separate. Averaging residual information does not by itself prove that the recovered signal is biological or partner-specific. |
| Model family | Parametric low-rank factors with ranks 0–3 | Shared signed rank-one G versus a symmetric spectral-shrinkage G, with G = 0 available. Multiple positive and negative modes are allowed in the spectral model. This is a substantive model change. |
| Regularization | Registered factor-parameter penalties | Penalties on the resulting G: nuclear norm plus squared Frobenius norm. The loss is mean-edge squared error against the mean replicate training residual, plus tau/n times the nuclear norm and nu/n² times the squared Frobenius norm. Tau fractions are 0.3, 0.1, 0.03 and 0.01 of the subset-specific zero-solution threshold; nu is 0.1 or 1. |
| Residual construction | Original joint procedure | Fit M0 on each inner-training subset, fit G to that subset's training residuals, then validate the complete M0 + G prediction on hidden inner pairs. Refit the same sequence on all permitted outer training data after selection. Training residuals are not treated as held-out performance. |
| K562 evaluation and budget | One 632-pair formal test; five training budgets; 20% the sole primary budget | Three development folds over a different 632-pair set from the original non-test pool. Compare low training counts 421/421/422 with high counts 2,317/2,317/2,318 on the same fold-specific evaluation pairs. This tests an information-budget hypothesis and does not replace the original primary endpoint. |
| Structural diagnostic | Fixed fully held-out rectangle error | In addition to raw prediction error, report error after eliminating free endpoint-additive effects using the test graph incidence matrix. This projection diagnostic is not the original registered rectangle endpoint, and eliminating additive effects alone does not exclude all nonlinear-baseline explanations. |
| Follow-up environment | Jurkat reserved for later work and a dated follow-up protocol | After K562 development, freeze the complete procedure locally and fit it on Jurkat training data. Rank one is primary; spectral shrinkage is secondary. This is a subsequent held-environment assessment under revised rules. |

The revised method-selection rule uses three training-only folds and pooled original-unit SSE. Candidates must pass numerical eligibility checks; G = 0 remains available. The zero baseline is replaced only when the selected correction improves inner SSE beyond the fixed numerical tolerance. Evaluation truth does not choose the correction family, penalties or initialization after the freeze. Exact implementation, split identities, loss normalization, convergence criteria and tie-breaking rules are retained in each stage's code and plan files.

## 3. Data use and chronological record

All listed times are UTC on 18 October 2026; China Standard Time is UTC + 8 hours. These are local execution records, not externally witnessed OSF deposit timestamps.

| Event | Recorded time | Implication |
|---|---|---|
| Original K562 formal test outcomes opened | 12:43:30 UTC (20:43:30 CST) | The original 632-pair formal test had been scored before subsequent development. It is inaccurate to describe all K562 answers as never inspected. |
| Subsequent K562 diagnosis and model development | After the original scoring and before the Jurkat freeze | The revised workflow was developed using K562 results. These analyses remain development analyses even though their new outer folds hide each evaluation pair from fitting. |
| Initial Jurkat freeze | 15:31:12 UTC (23:31:12 CST) | The initial protocol preceded numeric preparation. The first asset check stopped because GINS4 lacked usable single-perturbation measurements in both replicates. |
| Availability-only Jurkat amendment and new freeze | 15:32:20 UTC (23:32:20 CST) | Remove GINS4 without replacement, before combination outcome summaries, fitting or scoring. Preserve the initial version and amendment. The remaining panel has 79 genes. |
| Jurkat prediction lock | 15:32:52 UTC (23:32:52 CST) | Predictions and selection records were saved and hashed before the scorer opened the held-out truth. |
| Jurkat scoring | 15:33:24 UTC (23:33:24 CST) | Evaluate the frozen procedure; no outcome-driven model or grid changes followed. |

Two different K562 sets each contain 632 pairs. The original formal test set is disjoint from the later 632-pair development evaluation set. The latter is drawn from the original 2,528-pair non-test pool. Original formal pairs did not enter later fitting or scoring, but their earlier aggregate results were already known. A later script's statement that it did not access formal outcomes refers to that script's data boundary, not to the entire study's history.

Public source files were available before these analyses. Computational preparation necessarily read and aggregated source values and wrote sealed evaluation arrays. The Jurkat trainer read training outcomes and label-free query inputs; it did not read those sealed evaluation outcomes. Thus the evidence concerns a recorded separation of fitting, prediction lock and evaluation, not independent third-party custody or newly collected data.

The amended Jurkat protocol uses 3,081 eligible unordered pairs: 2,464 training and 617 test pairs. Both replicates, guide combinations and equivalent directions of a pair share its partition. Its three inner folds are for training-side selection, not three independent confirmation experiments.

## 4. How the results should be described

The original registered positive-evidence criteria were not met. Positive-rank optimization failures limit what the original zero gains imply about biology. Those outcomes remain reported, including the original progression gate that was not passed. Subsequent K562 development supplied a new rationale for the revised Jurkat assessment; it does not convert the original run into a successful registered test.

The two-stage archive supports a K562 development comparison and a locally frozen Jurkat held-environment confirmation. Jurkat is another cell environment within the same original study. Parameters were re-estimated from Jurkat training observations. Accordingly, the result concerns generalization of the fitting workflow, not direct transfer of K562-fitted parameters, independent-study replication, RNA validation or prospective collection of new biological measurements.

Jurkat rank one remained the primary method; the secondary spectral result does not replace it. Nominal gene-reweighting intervals are conditional sensitivity summaries for the observed data and fitted models, not coverage-calibrated uncertainty over new experiments. High training coverage is part of the demonstrated scope. Neither this analysis nor its deposit establishes a 75% biological hit rate, an exact biological rank, or an asymptotic measurement-complexity law.

This upload records completed results and the associated amendments. It is not evidence that the revised models or Jurkat protocol were uploaded to OSF before their execution. The earlier local freeze and prediction-lock records are retained so their actual sequence can be inspected. The original registration is preserved unchanged.

## 5. Exact accompanying archive and reproduction

- Archive: `WitnessCell_Horlbeck_TwoStage_Repro_20261018.zip`
- Size: 2,491,533 bytes.
- SHA-256: `b60e01230f3f43f156271697fe1fdc52b6c4ce5f68b522c3664c6e4aa6ff4082`
- Contents: processed numeric inputs, pair partitions, final-stage source code, plans, local freeze/history records, saved predictions, results, checks and reproduction commands for both stages.
- Entry point after extraction: `python3 reproduce.py --out ./reproduced`; verification instructions are in `README_zh.md` and `verify.py`.

The compact archive omits the large original source gzip and regenerable intermediate fit files. It is not the complete history of every earlier model-development run. Original registered-run outcomes and limitations are summarized in `REGISTERED_PLAN_STATUS.md`; their separate original result archive remains retained. The current deposit comprises this amendment, that status document, and the unchanged two-stage reproduction archive.
