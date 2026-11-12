# Han native functional refitting — fixed analysis v1

Written 2026-11-09 before fitting the new native Han predictors. Han outcomes and earlier transfer/attribution results have already been examined. This is a fixed-procedure reanalysis of an independent study, not untouched prospective confirmation or external preregistration. Source manuscript scores remain unchanged.

## Question and hierarchy

Can the frozen Jurkat functional learning procedure learn useful partner corrections from allowed Han training pairs, beyond its own fitted functional baseline? Does the increment persist over a reasonable target-fitted nonlinear reference? Can the already defined component control use this correction without shrinking its endpoint contribution?

Primary: reconstructed two-repeat Han panel, native Horlbeck/Jurkat baseline, selected signed rank-one correction, endpoint-eliminated squared prediction error relative to that same baseline. The primary correction is evaluated before any new control. Secondary: selected spectral correction, full error, original-vs-controlled and partner-off comparisons. Attribution sensitivity: replace only the baseline with v37's saturating-single-sum plus endpoint reference. Readout sensitivity: published aggregate panel, fitted as one response, never counted as two independent repeats. No analysis is promoted after inspecting results.

## Inputs and splits

Use unchanged v37 eligible identities and outcomes: 20,359 reconstructed pairs and 20,791 published pairs; all available single-response inputs. Reconstructed input has 207 single identities, of which 206 occur in eligible pairs. Keep that original identity list. Published input has 205 identities. Five archived identity partitions, each with three outer whole-unordered-pair folds. All readings of a pair follow its identity. Source weights, previous Han transfer predictions, and external gene annotations are not fitted or used in the new native predictors.

Preparation may read the existing outcome arrays to write disjoint training and scoring files. Fit workers accept training responses plus outcome-free query identities only. Query outcomes are opened by a separate scorer after all predictions are locked. Record pair counts, degrees, repeat correlations and exact input hashes. Repeated partitions and real repeats do not create independent target studies.

## Frozen learning rule

Native baseline: unchanged `frozen_models.fit_model`, rank=0, alpha=10, single-effect standard deviation scaling per response, constrained endpoint/slope/curvature effects, as in the released Jurkat interface. Generalize the adapter's response loop to its supplied column count; this changes no two-repeat arithmetic.

Residual matrix: unchanged `spectral_models.py`, mean residual across the supplied response columns; signed rank-one primary and symmetric spectral secondary. Nuclear fractions {0.3,0.1,0.03,0.01}, Frobenius penalties {0.1,1}; original six starts for rank-one, original spectral warm starts across decreasing nuclear penalties inside each inner fold. Optimizer cap 3,000 and relative stationarity tolerance 1e-6. No new architecture, starts, grid expansion or response-dependent fallback.

Three inner whole-pair folds use the released interface hash `wc_spectral_smoke_v1|` plus the canonical `gene1|gene2`. Refit each baseline on each inner training subset. Select candidates by summed complete-response SSE across all inner held-out pairs and response columns; exact ties prefer larger nuclear fraction then larger Frobenius penalty. Include zero correction using the original improvement tolerance. A candidate must meet the original convergence requirements in all three inner folds. Retain all diagnostics, including unselected candidates. An engineering timeout or failed fit is an explicit computational failure, never a zero biological gain. Retrying an interrupted identical fit is allowed and logged; changing the scientific grid is not.

Saturation baseline: v37 formula B=z+Xb+c[z/(1+lambda|z|/r)-z], z=sum of target single responses, r=RMS of target singles. Lambda {0,.05,.1,.2,.5,1,2,5}; three internal whole-pair folds and complete SSE select lambda, ties prefer lower lambda. Least-squares endpoint and saturation coefficients fitted separately by response column. This baseline must be selected and refitted inside every relevant training subset, including partner inner folds. The v37 identity ordering `v37-inner|outer_seed|canonical_pair` is retained for these nested baseline folds. It is a prespecified alternative baseline, not an unchanged copy of the original full Jurkat pipeline.

## Control and comparisons

For the selected predictor, pool inner-held-out residuals r_k=P_k(Y_k-B_k) and corrections d_k=P_kD_k, and set alpha=clip[0,1](sum<r_k,d_k>/sum||d_k||^2), with alpha=0 when energy<=1e-14. Zero-selected models retain D=0. This estimates one common strength across response columns. The selected hyperparameter and alpha use the same inner data, as in the existing control development; the external fold evaluates that entire rule.

On each fixed outer query graph, compute the endpoint-orthogonal projection P_Q from identities only. Save B; original B+D; whole attenuation B+alpha D; partner-off B+(I-P_Q)D; component control B+(I-P_Q)D+alpha P_QD; and a secondary independently full-loss-optimized whole attenuation using gamma=clip[0,1](sum<Y_k-B_k,D_k>/sum||D_k||^2). No query outcomes enter these operations. Query-panel dependence is retained and disclosed.

## Lock, scoring and stopping

Save each completed inner fit and outer batch atomically with parameters, predictions, convergence, timing, hashes and resume state. Lock all 60 outer baseline/readout batches (2 baselines x 2 readouts x 5 partitions x 3 folds) before outer scoring. Preserve final and selected inner models and all candidate inner predictions/selection scores for exact replay. Numerical engineering equivalence tests use synthetic inputs or existing Jurkat examples, not Han outer accuracy.

Score full and endpoint-eliminated SSE, both repetitions separately where available, totals by partition, and pooled SSE over all five partitions. Percent improvement is 100*(SSE_reference-SSE_method)/SSE_reference, with the denominator named. Baseline-to-addition and correction-to-own-baseline are different comparisons. Report all models, zeros and signs. Do not treat repeated partitions as independent biological replicates or claim population significance from them.

Validate source/input hashes, pair disjointness, training-only selection, native-adapter equivalence, predictions from saved weights, independent QR/lstsq residualization, and the full-equals-partner absolute SSE change for component control vs original. Complete a bounded run and report the result even if it is negative. Do not retune after scoring.
