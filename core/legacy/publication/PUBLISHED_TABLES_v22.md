## Supplementary tables



### Supplementary Table S1 | Released-split endpoint decomposition

The zero-contrast reference is one. All-gene scoring is shown. E is the number of combination–seed evaluation records; d is the endpoint-eliminated dimension within that split. Joint backbone denotes the backbone paired with final v14, whereas saturation-only was independently selected. A dash denotes an unavailable matching prediction, not an inferred score.

| Dataset / seed | E / d | Final v14 | Joint backbone | Saturation-only |
| --- | --- | --- | --- | --- |
| Norman / 1 | 73 / 24 | 0.995140 | 0.995140 | 0.995140 |
| Norman / 2 | 66 / 16 | 1.000736 | 1.000736 | 1.003976 |
| Norman / 3 | 69 / 22 | 1.018454 | 1.018454 | 1.009911 |
| Wessels / 1 | 63 / 39 | 1.000000 | 1.000000 | 0.999904 |
| Wessels / 2 | 54 / 33 | 0.993565 | 0.993565 | 0.993565 |
| Wessels / 3 | 64 / 45 | 1.000000 | 1.000000 | 0.997968 |
| Schmidt / 1 | 22 / 8 | 0.999722 | 0.999722 | — |
| Schmidt / 2 | 19 / 4 | 0.997796 | 0.997796 | — |
| Schmidt / 3 | 19 / 5 | 0.996738 | 0.996738 | — |
| Replogle_exp6 / 1 | 35 / 10 | 1.000000 | 1.000000 | — |
| Replogle_exp6 / 2 | 29 / 4 | 0.998739 | 0.998739 | — |
| Replogle_exp6 / 3 | 27 / 2 | 0.998246 | 0.998246 | — |

The six available Norman/Wessels additive references equal one. The complete component table and additional scoring variants remain in the compact archive; the separate matched RNA comparison is Table S2.

**Dataset-level complete-response and projected effects**

Within each dataset, squared errors and truth energies are summed across the three seeds before forming ratios. Gains compare final v14 with its own jointly selected backbone; projected ratios use projected observed-response energy as denominator.

| Dataset | Records / d | Readouts | Full-response gain % | Projected ratio |
| --- | --- | --- | --- | --- |
| Norman | 208 / 62 | 5025 | 5.7433 | 1.005527 |
| Wessels | 181 / 117 | 5020 | 31.9396 | 0.998160 |
| Schmidt | 60 / 17 | 5005 | 0.0658 | 0.998415 |
| Replogle_exp6 | 91 / 16 | 5019 | 1.2342 | 0.999490 |

Projected gains against the paired backbone are zero within saved-array quantization. Wessels uses 5,020 raw readouts here, versus 4,995 non-target readouts in Table S2.

### Supplementary Table S2 | Matched RNA comparison, seeds and fitting durations

All methods are scored on the same held-out combinations, 11 same-fit rectangles and 4,995 readouts. Neural values are means of the three seed-specific scores; fixed predictions have one realization. MSE columns are multiplied by 1,000 for display. Positive gain means lower rectangle error than RNA M0.

| Model | Rectangle MSE × 1,000 | Gain vs M0 % | Response MSE × 1,000 | Projected ratio |
| --- | --- | --- | --- | --- |
| Measured-single additive | 4.082949 | 3.5945 | 3.717907 | 1.000000 |
| Saturation-only | 4.344639 | -2.5844 | 3.254583 | 1.038787 |
| RNA M0 | 4.235183 | 0.0000 | 1.281651 | 1.018972 |
| RNA M1 | 3.835457 | 9.4382 | 1.214510 | 0.918185 |
| RNA R | 3.559154 | 15.9622 | 1.162690 | 0.865165 |
| GEARS | 3.990340 | 5.7812 | 5.154565 | 0.976397 |
| scouter | 4.085221 | 3.5409 | 18.063703 | 1.000400 |

M1/R full-response gains against M0 are 5.2387%/9.2819%. Note 2.3 explains the neural complete-response errors; Note 2.4 defines the distinct rectangle and projected outcomes.

**Pairwise comparisons and individual optimization seeds**

Neural references are the mean score or the lowest-error seed; best-seed contrasts are descriptive.

| Reference | R gain vs mean/fixed % | R gain vs best seed % |
| --- | --- | --- |
| RNA M0 | 15.9622 | 15.9622 |
| RNA M1 | 7.2039 | 7.2039 |
| Measured-single additive | 12.8289 | 12.8289 |
| GEARS | 10.8058 | 8.3950 |
| scouter | 12.8773 | 12.5895 |

**Individual optimization seeds.** The same 55 biological combinations are reused; seed variation is not biological replication.

| Model / seed | Rectangle MSE × 1,000 | Response MSE × 1,000 | Projected ratio |
| --- | --- | --- | --- |
| GEARS / 0 | 4.027845 | 6.843357 | 0.966473 |
| GEARS / 1 | 3.885328 | 4.261756 | 0.969947 |
| GEARS / 2 | 4.057849 | 4.358583 | 0.992771 |
| scouter / 0 | 4.071770 | 1.717753 | 0.997512 |
| scouter / 1 | 4.090416 | 27.309634 | 1.001689 |
| scouter / 2 | 4.093476 | 25.163722 | 1.001998 |

**Selected final-fit epochs (folds 0–4).**

| Seed | GEARS | scouter |
| --- | --- | --- |
| 0 | 8, 7, 5, 7, 16 | 12, 4, 17, 5, 37 |
| 1 | 14, 14, 5, 17, 9 | 1, 5, 6, 3, 12 |
| 2 | 14, 9, 3, 5, 23 | 1, 5, 7, 11, 5 |

All GEARS fits selected schedule step 5. Scouter upstream commit: `fe2a73bdda02c9ec639bd55dabe758da2f7ff335`; see Note 2.3 for selection scope.

### Supplementary Table S3 | Reliability controls and readout-scale diagnostics

Positive percentages indicate lower error than M0. All models, folds and readouts are fixed. Randomized rows summarize conditional perturbations of the retained predictions, not additional biological experiments.

| Predictor/control | Rectangle-error reduction % | Complete-response reduction % |
| --- | ---: | ---: |
| M1 | 9.4382 | 5.2387 |
| Training-inner uniform amplitude | 9.6988 | 5.3343 |
| Reliability-weighted R | 15.9622 | 9.2819 |
| R relative to uniform amplitude | 6.9362 | 4.1700 |
| Weight permutation, median of 200 | 3.2895 | See retained draw table |
| Weight permutation, best of 200 | 8.0423 | See retained draw table |
| Same weights on M0 minus measured additivity | 0.2368 | −134.596 |
| Same-norm random increment, weighted, median of 200 | −0.704 | See retained draw table |

R improved on the uniform-amplitude control in all 11 rectangles. None of the 200 permuted-weight predictions matched R; none of the 200 weighted random increments beat M0. These counts are descriptive controls. The scalar coefficients were 0.838645, 1, 1, 1 and 1 for folds 0–4.

**Across-condition variance standardization and fixed-ten-readout deletion.**

All entries below are rectangle-error reductions relative to M0 under the same transformation. Negative values mean deterioration. Standardization uses only outer-training combination variances; deletion uses the previously fixed ten-gene list in Note 3.3.

| Model | Raw gain % | Standardized gain % | Ten-gene deletion gain % |
| --- | --- | --- | --- |
| Measured-single additive | 3.5945 | 0.0244 | -0.3639 |
| Saturation-only | -2.5844 | -0.0999 | -0.6658 |
| RNA M0 | 0.0000 | 0.0000 | 0.0000 |
| RNA M1 | 9.4382 | -8.6708 | -10.0314 |
| RNA R | 15.9622 | -0.1180 | -0.3798 |
| GEARS | 5.7812 | -21.2418 | -9.9918 |
| scouter | 3.5409 | -7.7030 | -0.4520 |

**Concentration of positive gains.** Denominators are defined in Note 3.3; the top-ten share is not a share of net improvement.

| Model | Improved readouts / 4,995 | Ten-gene share of positive gains % | Negative cancellation % |
| --- | --- | --- | --- |
| RNA M1 | 1596 | 82.2127 | 54.6826 |
| RNA R | 1929 | 92.1253 | 9.5238 |

R also gave +2.0407% under within-condition cell-variance standardization. These variances measure a different scale from across-condition response variance. Cross-readout rectangle Spearman/direction agreement were −0.003834/49.35% (M0), 0.027902/50.91% (M1) and 0.026450/50.19% (R).

### Supplementary Table S4 | K562 training budgets, pooled and fold-specific outcomes

All gains are percentage reductions relative to M0 fitted on the same allowed training pairs. K562 combines SSE from three outer development folds; Jurkat uses one held-out evaluation. Each pair has two source replicate measurements.

| Environment / budget | Training pairs | Evaluation pairs | M0 pooled MSE | Rank-one / spectral gain % |
| --- | --- | --- | --- | --- |
| K562 low | 421 / 421 / 422 | 632 total | 0.00425513 | 0.2652 / -0.2879 |
| K562 high | 2,317 / 2,317 / 2,318 | 632 total | 0.00382441 | 11.0292 / 11.3865 |
| Jurkat high | 2,464 | 617 | 0.00092919 | 6.4308 / 8.1606 |

**Endpoint-eliminated error and replicate direction.** Projected gains use projected M0 error as denominator and sum within-fold errors.

| Environment / budget | Contrast dimension | Rank-one gain % | Spectral gain % | Positive raw-error replicate–fold counts |
| --- | --- | --- | --- | --- |
| K562 low | 132 + 131 + 130 | -0.0326 | -0.8591 | 2/6; 3/6 |
| K562 high | 132 + 131 + 130 | 11.6838 | 10.6710 | 6/6; 6/6 |
| Jurkat high | 538 | 6.1750 | 5.8975 | 2/2; 2/2 |

Positive counts are rank one; spectral. High-budget M0 improved 10.1225% over low-budget M0; the tabulated partner gains use their own same-budget baseline.

**K562 fold-specific gains and selected settings**

Fold labels 0–2 match the stored source identifiers. Parentheses give replicate 1 and replicate 2 gains; all entries are percentages versus the same-budget M0. Evaluation counts are 211, 211 and 210.

| Fold | Budget | Rank-one gain (replicates) % | Spectral gain (replicates) % | Projected gains: rank one / spectral % |
| --- | --- | --- | --- | --- |
| 0 | low | 0.0000 (0.0000, 0.0000) | -1.1441 (-0.9392, -1.5961) | 0.0000 / -2.6890 |
| 0 | high | 5.7110 (8.2111, 0.6430) | 5.4163 (6.2223, 3.7825) | 6.5305 / 2.9572 |
| 1 | low | 0.0000 (0.0000, 0.0000) | 0.0159 (0.1854, -0.2272) | 0.0000 / -0.3035 |
| 1 | high | 16.8664 (12.1849, 23.6444) | 14.1850 (11.1525, 18.5756) | 18.6244 / 15.8522 |
| 2 | low | 0.7412 (0.0487, 1.6475) | 0.1415 (0.1715, 0.1021) | -0.0948 / 0.3397 |
| 2 | high | 10.0215 (7.4512, 13.3659) | 13.5532 (10.3977, 17.6591) | 9.8256 / 13.0592 |

**Training-only choices.** Entries are \((\tau/\tau_0,\nu)\); zero means the selected predictor retained M0. The best penalized fit was saved separately when zero was chosen.

| Fold | Low: rank one | Low: spectral | High: rank one | High: spectral |
| --- | --- | --- | --- | --- |
| 0 | Zero | (0.01, 0.1) | (0.03, 0.1) | (0.3, 0.1) |
| 1 | Zero | (0.03, 1) | (0.03, 0.1) | (0.3, 0.1) |
| 2 | (0.3, 1) | (0.01, 1) | (0.03, 0.1) | (0.3, 0.1) |

All 288 inner fits, 12 outer correction fits and 144 objective-order checks met the recorded criteria.

### Supplementary Table S5 | Jurkat primary and secondary results

The panel contains 79 genes, 2,464 training pairs and 617 test pairs. Rank one is primary and spectral correction is secondary.

| Quantity | Rank one: primary | Spectral: secondary | M0 replicate-mean reference |
| --- | --- | --- | --- |
| Pooled MSE | 0.00086943 | 0.00085336 | 0.00112689 |
| Pooled gain % | 6.4308 | 8.1606 | -21.2774 |
| Endpoint-eliminated gain % | 6.1750 | 5.8975 | 0.9177 |
| Replicate 1 gain % | 4.4542 | 8.0730 | -19.4882 |
| Replicate 2 gain % | 8.8845 | 8.2692 | -23.4982 |
| Conditional reweighting interval % | [1.7335, 12.1290] | [3.0017, 13.1873] | [-40.5760, -8.9395] |
| Minimum gain after one-gene deletion % | 4.9492 | 7.4284 | -22.9565 |
| Positive one-gene deletions | 79/79 | 79/79 | 0/79 |

Both partner families selected \((\tau/\tau_0,\nu)=(0.3,0.1)\). All 48 inner fits, 24 objective-order probes and two final correction fits passed. The 10,000 paired gene-reweighting draws condition on the observed measurements and fitted predictions; their intervals exclude refitting and are sensitivity summaries rather than calibrated coverage guarantees. The 79 overlapping deletions are likewise descriptive.

### Supplementary Table S6 | Original plans, amendments and result identities

Times below are UTC on 18 October 2026; China Standard Time is UTC + 8 hours. Local model freezes are distinct from external registration and later results deposit.

| Study stage | Fixed comparison or event | Outcome and interpretation |
| --- | --- | --- |
| RNA M1 development | Original partner-kernel contrast criterion: at least 10% improvement | 9.4382%; original criterion not met |
| RNA R development and subsequent diagnostic freeze | Same 55-pair panel; readout reliability added, then retained for diagnostics | 15.9622%; developmental result, not independent RNA confirmation |
| Original registered K562 | Separate-repeat ranks 0–3, five budgets, original 632-pair test; 20% primary | Both primary gains zero; 266/270 positive-rank inner fits failed at primary budget |
| Original K562 scoring, 12:43:30–12:43:31 | Formal outcomes evaluated | Known before revised K562 development |
| Revised K562 development | New matrix penalties/shared correction; distinct 632-pair evaluation from non-test pool | Three-fold budget comparison; original formal-test pairs excluded |
| K562 completion plan, 15:09:48 | Finalized smoke procedure applied to remaining folds | Fold 0 previously examined; three-fold results remain development evidence |
| Initial Jurkat freeze, 15:31:12 | Fixed procedure and original 80 gene identities | Availability check stopped for missing GINS4 singles |
| Availability-only amendment, 15:32:20 | Remove GINS4 without replacement and refreeze | 79 genes; before fitting or combination outcome summaries |
| Jurkat prediction lock, 15:32:52 | Predictions and selection records saved and hashed | Trainer had not read sealed evaluation outcomes |
| Jurkat scoring, 15:33:24 | Evaluate 617 held-out pairs | Primary rank-one criterion met; spectral remains secondary |
| Subsequent OSF results deposit | Two-stage archive and disclosure documents | Preserves follow-up, without retroactively preregistering revised methods |

The original five K562 budgets used 158, 316, 632, 1,264 and 2,528 training pairs. All selected the baseline, with incremental gains zero and pooled MSEs 0.00596059, 0.00450399, 0.00396413, 0.00369680 and 0.00356008. Zero denotes identical selected predictions, not a biological effect estimated to be exactly zero. Original outcomes and computational records remain in the registered-plan status documents; subsequent development does not replace them.


### Supplementary Table S7 | Complementary breadth and selection results

These comparisons use earlier endpoint or analytic-risk implementations and are not additional partner-learning replications.

**Released benchmark and combination-only weighted analysis**

| Outcome | Value and scope |
| --- | --- |
| Released aggregation | Mean top-100 operation rank 3.803; first among 15 published entrants plus WitnessCell |
| Benchmark population | 654 condition–seed units, including 540 combinations and 114 unseen singles |
| Combination-only WMSE ratio vs released additive training mean | 0.7513 [0.7031, 0.8019] |
| Combination-only weighted direction-\(R^2\) increment | +0.1908 [+0.1648, +0.2171] |
| All 654-unit WMSE ratio | 0.6829 [0.6255, 0.7434]; includes unseen singles |

**All-query non-target loss ratios**

| Dataset | Comparator | MSE ratio [95% pair interval] | WMSE ratio [95% pair interval] |
| --- | --- | --- | --- |
| Norman | Released training mean | 0.8701 [0.8304, 0.9070] | 0.9439 [0.8972, 1.0001] |
| Norman | Shared-endpoint additive | 0.9049 [0.8563, 0.9475] | 0.9996 [0.9464, 1.0622] |
| Norman | Independent saturation | 0.9355 [0.9023, 0.9669] | 0.9612 [0.9256, 0.9981] |
| Norman | Paired factorized backbone | 0.9392 [0.9083, 0.9685] | 0.9677 [0.9377, 0.9997] |
| Wessels | Released training mean | 0.6459 [0.6036, 0.6860] | 0.6976 [0.5650, 0.8610] |
| Wessels | Shared-endpoint additive | 0.6459 [0.6036, 0.6860] | 0.6976 [0.5650, 0.8610] |
| Wessels | Independent saturation | 0.6764 [0.6318, 0.7180] | 0.6185 [0.5184, 0.7330] |
| Wessels | Paired factorized backbone | 0.6675 [0.6238, 0.7084] | 0.6234 [0.5228, 0.7390] |

**Fixed-system Wessels residual-library contribution against the same backbone**

| Loss | Quarter library gain % | Half library gain % | Full library gain % | Full-minus-quarter percentage points |
| --- | --- | --- | --- | --- |
| MSE | 19.10 [16.17, 22.32] | 27.42 [23.81, 31.22] | 33.25 [29.28, 37.49] | 14.15 [12.42, 16.06] |
| WMSE | 26.81 [17.07, 34.52] | 33.91 [22.99, 42.74] | 37.66 [26.18, 47.74] | 10.85 [7.52, 14.36] |

**Separate analytic-risk threshold experiment**

| Outcome | Original MSE endpoint | Subsequent WMSE complete case |
| --- | --- | --- |
| Queries / canonical pairs | 213 / 33 | 207 / 32 |
| Accepted records | 110 | 106 |
| Pair-balanced coverage | 49.48% | 48.95% |
| Accepted/all loss ratio | 0.5343 | 0.4775 |
| Pair-cluster percentile interval | [0.3646, 0.7591] | [0.3175, 0.7071] |

The six weighted-unevaluable queries are from one pair; no replacement weights or outcomes were introduced. Pair-balanced coverage and losses use the same record-weight definition in Note 5.3.

**Released-predictor empirical ranking at 50% nominal coverage**

| Dataset / pool | Support vs categorical-exposure MSE improvement % | 95% conditional pair-bootstrap interval |
| --- | ---: | --- |
| Norman, all queries | −2.09 | [−11.23, 10.84] |
| Norman, seen1 reselection | 1.83 | [−18.42, 16.10] |
| Wessels, all queries | 8.27 | [1.75, 15.30] |
| Wessels, seen1 reselection | 17.57 | [5.52, 26.11] |

Exposure-aware controls narrow the interpretation of earlier gains against effect/variance-only ranking. Coarse effect-adjusted comparisons did not establish a general effect-independent advantage. Full retained coverage/loss grids are in machine-readable sources.


### Supplementary Table S8 | Source and reproduction map

Bare CSV filenames and the `RNA_controls/` directory are relative to this paper package's `source_data/`. Stage paths refer to the compact growth archive, and `official12/` refers to the compact second-step archive. Large original assets are identified separately rather than represented as included; their identities and extraction checks are retained in the source manifests and the compact archives. Original and amended protocol records are in `registration/`.

| Scientific object | Retained source | Reproduction level |
| --- | --- | --- |
| Released-split predictions and decomposition | `Step2_Fig5_Official12.csv`; compact archive `official12/` | Per-split sufficient-statistic replay; original-vector scorer uses original archives |
| Matched RNA comparison | `Step2_Fig5_Matched55.csv`; `Step2_Fig5_SeedMetrics.csv`; `Step2_ED5_Folds.csv`; `Step2_ED5_Rectangles.csv` | Aggregate, seed, fold and rectangle scoring |
| Readout contributions | `Step2_ED5_GeneGains.csv`; `RNA_displayed_signed_examples.csv`; retained scale outputs | Positive/negative gain decomposition and fixed-readout deletion |
| Uniform amplitude and weight randomization | `RNA_controls/`; original `confirmation_20261018/controls/` | Fixed-model outputs, selection metadata, stratified draw results and code |
| Growth core and selection | `stages/K562/core/{frozen_models,spectral_models,run_smoke}.py` | Processed-input refitting, matrix objectives and complete-pipeline selection |
| K562 partitions and predictions | `stages/K562/cases/outer*/inputs`; `stages/K562/SUMMARY.json`; `Horlbeck_predictions.csv` | Exact pair membership, two budgets and three outer evaluations |
| Jurkat protocol and results | `stages/Jurkat/{PLAN.json,FREEZE.json,PREDICTIONS_LOCKED.json,SUMMARY.json}` | Frozen rule, availability amendment, fitted parameters and held-out scoring |
| Jurkat sensitivities | `Horlbeck_Jurkat_gene_deletion.csv`; `Horlbeck_Jurkat_node_reweighting.csv`; `Horlbeck_Jurkat_intervals.csv` | Conditional reweighting and fixed-prediction deletion checks |
| RNA four-condition examples | `biological_examples/rna/`; all 54,945 candidate coordinates and selected four-corner values | Signed-response reconstruction, training-weight provenance and post-hoc display selection |
| Jurkat named combinations | `biological_examples/jurkat/`; all 617 pair scores and 9,493 shared-endpoint comparisons | Two-repeat response display and complete-pool context from frozen primary predictions |
| Earlier prediction and ranking | Retained prediction-ratio, risk-comparison and witness-library source CSVs | Complementary full-response and selection analyses |

\clearpage

### Supplementary Table S9 | RNA four-corner responses and signed contrasts

Responses are relative to the common control mean in source-processed expression units. Rectangle and fold identifiers here follow the one-based display identifiers in the example tables. Signs specify the fixed contrast orientation; exact precision, control means, cell-half summaries and uniform-amplitude predictions are supplied in the source CSVs. The rows describe selected explanatory cases (Note 3.4), not additional experimental replicates.

| Readout / rectangle / fold | Combination | Sign | Cells | Observed | M0 | M1 | R |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| S100A9 / 2 / 2 | EP300 + FLI1 | + | 104 | -0.495284 | -1.852908 | -1.468034 | -1.497866 |
| S100A9 / 2 / 2 | EP300 + INTS1 | − | 160 | -0.207039 | -0.137643 | -0.047588 | -0.054568 |
| S100A9 / 2 / 2 | FLI1 + GFI1 | − | 89 | -0.327134 | 1.168718 | 0.336637 | 0.401131 |
| S100A9 / 2 / 2 | GFI1 + INTS1 | + | 256 | 2.994962 | 2.347484 | 2.923488 | 2.878842 |
| LYZ / 6 / 1 | GSE1 + HDAC3 | + | 135 | 1.819147 | 1.263848 | 1.366401 | 1.354516 |
| LYZ / 6 / 1 | GSE1 + MED24 | − | 158 | 0.235604 | 0.418220 | 0.248548 | 0.268211 |
| LYZ / 6 / 1 | HDAC3 + SMARCD1 | − | 267 | 0.185446 | 0.180293 | 0.018331 | 0.037100 |
| LYZ / 6 / 1 | MED24 + SMARCD1 | + | 192 | -0.305902 | -0.799841 | -0.493738 | -0.529212 |
| LYZ / 5 / 5 | GFI1 + RING1 | + | 151 | 1.778080 | 2.009728 | 2.224961 | 2.199126 |
| LYZ / 5 / 5 | GFI1 + SMARCD1 | − | 143 | 1.443725 | 1.555983 | 1.358412 | 1.382127 |
| LYZ / 5 / 5 | RCOR1 + RING1 | − | 74 | 0.346673 | 0.515970 | 0.464533 | 0.470707 |
| LYZ / 5 / 5 | RCOR1 + SMARCD1 | + | 95 | -0.102989 | -0.041853 | -0.040526 | -0.040686 |

**Signed contrasts reconstructed from the four rows for each example.** The same-fit ordering is \((A,B),(A,D),(C,B),(C,D)\), with signs \(+,-,-,+\). Training weights are shared across that fitted model's queries for the displayed readout.

| Readout / rectangle | Observed contrast | M0 | M1 | R | Training weight |
| --- | ---: | ---: | ---: | ---: | ---: |
| S100A9 / 2 | 3.033851 | -0.536498 | 1.166406 | 1.034414 | 0.922490 |
| LYZ / 6 | 1.092195 | -0.134506 | 0.605784 | 0.519993 | 0.884111 |
| LYZ / 5 | -0.115307 | -0.104078 | 0.361490 | 0.305606 | 0.879965 |

The adverse LYZ contrast is close to the M0 prediction before correction; its relative percentage deterioration is consequently unstable as a summary. The absolute contrast values and the complete candidate table are retained instead.

\clearpage

### Supplementary Table S10 | Named Jurkat combinations in both experimental repeats

Values are original growth-phenotype \(\gamma\) units. The primary rank-one predictor is M0 plus the frozen shared correction. Pair identities and both repeats were held out together. Cases were selected after evaluation as described in Note 4.8; NIP7–POLN is the largest absolute pooled-error deterioration.

| Combination | Repeat | Observed | M0 | M0 + rank one | Shared correction |
| --- | ---: | ---: | ---: | ---: | ---: |
| PGK1–TTI2 | 1 | -0.296624 | -0.255375 | -0.267979 | -0.012603 |
| PGK1–TTI2 | 2 | -0.345788 | -0.292557 | -0.305160 | -0.012603 |
| OGFOD1–TTI2 | 1 | -0.030502 | -0.060076 | -0.049662 | 0.010415 |
| OGFOD1–TTI2 | 2 | -0.041278 | -0.077034 | -0.066620 | 0.010415 |
| MED1–RPL10 | 1 | -0.336827 | -0.310029 | -0.324868 | -0.014838 |
| MED1–RPL10 | 2 | -0.330016 | -0.300927 | -0.315765 | -0.014838 |
| NIP7–POLN | 1 | -0.508615 | -0.408077 | -0.398085 | 0.009993 |
| NIP7–POLN | 2 | -0.508551 | -0.453696 | -0.443703 | 0.009993 |

For TTI2–OGFOD1 minus TTI2–PGK1, observed differences were 0.266123 and 0.304509 in repeats 1 and 2. M0 gave 0.195299 and 0.215523, and the primary partner model gave 0.218317 and 0.238540. The common pair correction is shared by construction; improved agreement with the separately withheld observations, rather than repeat correlation between the predictions, is the relevant display result. This difference retains the PGK1 and OGFOD1 endpoint terms.


### Supplementary Table S11 | Original-system integration, fixed-G controls

Positive values indicate percentage error reduction relative to the unchanged original v14. Each entry is the geometric aggregation of three original seed-specific loss ratios.

| Control | Norman complete | Norman projected | Wessels complete | Wessels projected |
|---|---:|---:|---:|---:|
| Original v14 / zero | +0.0000% | +0.0000% | +0.0000% | +0.0000% |
| Unit correction | -3.2257% | +2.6758% | +0.4921% | +0.1680% |
| Uniform shrinkage | +0.3410% | +0.9191% | +0.2535% | +0.1214% |
| Amplitude-only | -1.1021% | +2.5682% | +0.6993% | +0.2004% |
| Half-agreement / amplitude-noise | -1.1653% | +2.3473% | +0.6335% | +0.1614% |
| Training-selected policy | +0.3410% | +0.9191% | +0.6283% | +0.1716% |

Norman: 208 pair–seed records, 118 distinct pairs. Wessels: 181 records, 106 distinct pairs. The complete and projected outcomes use the distinct non-target masks defined in Note 7.

### Supplementary Table S12 | Six original split selections and realized increments

| Dataset / seed | Test pairs | Selected policy | Complete gain | Projected gain | Test-graph contrast dimension |
|---|---:|---|---:|---:|---:|
| Norman / 1 | 73 | Uniform shrinkage | +0.3228% | +0.6241% | 24 |
| Norman / 2 | 66 | Original v14 / zero | +0.0000% | +0.0000% | 16 |
| Norman / 3 | 69 | Uniform shrinkage | +0.6991% | +2.1213% | 22 |
| Wessels / 1 | 63 | Uniform shrinkage | +0.0514% | -0.0011% | 39 |
| Wessels / 2 | 54 | Half-agreement / amplitude-noise | +1.5144% | +0.4513% | 33 |
| Wessels / 3 | 64 | Amplitude-only | +0.3130% | +0.0641% | 45 |

Selection used training predictions only. A uniform coefficient of zero in Norman seed 2 equals the unextended core. Every inner selection log, endpoint-availability stratum and output array is retained in the original-system run package.

### Supplementary Table S13 | Independent-source RNA, complete frozen comparison

All entries are percentage reductions relative to M0; positive values indicate lower error. A negative value for a comparator indicates that it has larger error than M0, not the reciprocal improvement of M0.

| Method | Full response | Projected primary | Full G2M | Projected G2M |
|---|---:|---:|---:|---:|
| M0 / zero | +0.00000% | +0.00000% | +0.00000% | +0.00000% |
| Unit | +0.22204% | -0.12901% | +3.02328% | +6.60838% |
| Uniform | +0.26980% | -0.03260% | +2.56969% | +4.03391% |
| Amplitude-only | +0.25987% | -0.03322% | +2.45074% | +3.60235% |
| Amplitude-noise | +0.07927% | -0.00355% | +0.56885% | +0.39142% |
| Selected | +0.26980% | -0.03260% | +2.56969% | +4.03391% |
| Measured additive | -30.60428% | +0.00473% | +19.14509% | +0.01230% |
| Saturation-only | -30.26808% | +0.00000% | +19.21281% | +0.00000% |

All 45 pairs were withheld once in the three fixed outer folds. The gene program was specified before this evaluation and uses 52 of 54 exact unique mappings. The primary partner result is negative; the secondary application does not outperform measured additivity in complete error.

### Supplementary Table S14 | Independent RNA fold identity, selections and conditional sensitivity

| Fold | Train / test pairs | Selected c | Core gamma | Full gain | Projected gain |
|---|---|---:|---:|---:|---:|
| 0 | 30 / 15 | 0.50 | 0.77691 | +0.20136% | +0.02930% |
| 1 | 30 / 15 | 0.75 | 0.78138 | +0.41490% | -0.04779% |
| 2 | 30 / 15 | 0.75 | 0.76559 | +0.19528% | -0.07727% |

All folds selected uniform shrinkage, lambda=10, alpha=0.5 and eta=1. Each test graph has rank ten and five contrast dimensions. The complete core is reselected within every allowed nested fit, not selected on these three outer outcomes.

| Sensitivity | Complete gain | Projected gain | Interpretation |
|---|---:|---:|---|
| Positive-node multiplier, central 95% | +0.10351% to +0.44214% | −0.10916% to +0.07366% | Fixed predictions; projections refitted to weights |
| Delete one endpoint, minimum to maximum | +0.24285% to +0.35761% | −0.08717% to +0.05049% | Ten deletions; test projection rebuilt |

These ranges are conditional descriptions rather than biological-replicate confidence intervals. Original protocol, implementation correction and freeze receipts are retained in the source record.


### Supplementary Table S15 | Dixit program development and fixed-reference delivery

All rows evaluate the same 45 outer queries and the same program units. Positive gains mean lower error than the independently fitted scalar saturation baseline. The selected full-roster reference is the primary development result; all other rows retain their specified comparator identities.

| Construction | Total SSE | Projected SSE | Total gain % | Projected gain % | Query batching invariant |
|---|---:|---:|---:|---:|---|
| Scalar saturation | 6.538922596 | 2.307093784 | 0.0000 | 0.0000 | Yes |
| Unprotected, matched ridge | 6.695915797 | 1.977192992 | -2.4009 | 14.2994 | Yes |
| Unprotected, own inner selection | 6.550751732 | 1.961337847 | -0.1809 | 14.9866 | Yes |
| Current-query projection | 6.209021805 | 1.977192992 | 5.0452 | 14.2994 | No |
| Fixed full-roster reference | 6.390677159 | 1.977192992 | 2.2671 | 14.2994 | Yes |
| Fixed unmeasured-pool reference | 6.278276738 | 2.046447926 | 3.9861 | 11.2976 | Yes |
| Reference kernel, full residual loss | 6.710468769 | 2.234667612 | -2.6235 | 3.1393 | Yes |
| Reference kernel, training-projected loss | 6.390677159 | 1.977192992 | 2.2671 | 14.2994 | Yes |

The two fixed-reference implementations with 2.2671% total gain are mathematically equivalent and are counted once. Keeping projected gain unchanged at fixed coefficients follows from endpoint cancellation. Conditional sensitivity and adverse endpoint deletion are in Note 9.4.


\clearpage

### Supplementary Table S16 | Frozen Combi-seq complete comparison

All methods predict the same 374 held-out pairs. Gains in the table are percentages relative to the shared endpoint baseline B+E. Positive values mean lower error. Program units and outer-fold scaling are identical across methods; projected losses pool the 266 fixed test-graph contrast dimensions.

| Method | Complete SSE | Projected SSE | Complete gain % | Projected gain % |
|---|---:|---:|---:|---:|
| Measured-single additive | 177.237788737 | 125.915191495 | -3.006695 | +0.475131 |
| Scalar saturation B | 172.818152532 | 126.516309700 | -0.438100 | +0.000000 |
| Endpoint baseline B+E | 172.064339459 | 126.516309700 | +0.000000 | +0.000000 |
| Previous \(B+g_F\) | 172.876243066 | 126.554568339 | -0.471860 | -0.030240 |
| Principal \(B+E+g_F\) | 172.104399943 | 126.554568339 | -0.023282 | -0.030240 |

The principal system reduces complete error by 0.413008% relative to scalar saturation, while the endpoint baseline alone reduces it by 0.436189%. The preceding saturation-plus-program assembly gives -0.033614% versus saturation. The principal partner increment is negative, and no alternative comparator changes that decision.


\clearpage

### Supplementary Table S17 | Combi-seq folds, model selection and run-target outcomes

Fold outcomes compare the principal system \(B+E+g_F\) with the same endpoint baseline \(B+E\). All three folds selected alpha=0.5. Lambda concerns the unchanged partner operator; noise and gamma concern the shared endpoint baseline.

| Fold | Train / test | Lambda | Noise / gamma | Complete % | Projected % |
|---|---|---:|---|---:|---:|
| 0 | 249 / 125 | 10 | 0.3 / 0.997148 | +0.035808 | -0.000757 |
| 1 | 249 / 125 | 0.1 | 0.01 / 0.445541 | -1.798455 | -0.767833 |
| 2 | 250 / 124 | 0.1 | 0.01 / 0.327347 | +0.433845 | -0.500088 |

All three retained experimental-run targets were scored against the same predictions fitted to pooled training means. They were not used to select different models or independent success thresholds. Each row includes the same 374 pair identities.

| Run target | Endpoint SSE | Principal system SSE | Complete gain % |
|---|---:|---:|---:|
| rep1 | 82.344587094 | 82.541384556 | -0.238993 |
| rep2 | 1372.160073814 | 1372.332694951 | -0.012580 |
| rep3 | 57.254705045 | 57.005467898 | +0.435313 |

Only one of three run targets improves. The principal frozen decision is not confirmed. The conditional complete-gain range is -1.098310% to +0.688105%; this is not a biological-replicate confidence interval.


\clearpage

### Supplementary Table S18 | Calibration of the measurement reference improves native M0 response prediction

| Input estimator | Source libraries | Original-scale MSE | Error reduction versus raw (%) |
|---|---:|---:|---:|
| Raw cell-weighted mean | 9 | 0.001623698041 | 0.0000 |
| Direct fixed-composition standardization | 7 | 0.001476556934 | 9.0621 |
| Shared-offset calibration | 7 | 0.001436203065 | 11.5474 |
| Shared-offset calibration | 9 | 0.001413839285 | 12.9247 |

All arms predict the same seven-library equal-composition target. Mean squared error averages 37,292 non-target RNA coordinates, two held-out technical-half references and 21 eligible conditions, with equal condition weight. Five outer folds were drawn from the historical development pool of 25 combinations; each used 20 combinations and ten singles for fitting, with preprocessing and M0 selection repeated inside three inner folds. Seven-library calibration reduced error by 2.7330% relative to direct standardization under the same available-library information. The corresponding nine-library difference was 4.2476% and additionally involved two source libraries. Nine-library calibration improved 20/21 conditions and 5/5 folds versus raw inputs; neither conditions with shared endpoints nor overlapping folds or technical halves are independent biological replications. This table concerns complete-response prediction without a partner head.

\clearpage

### Supplementary Table S19 | Complete five-task shared-feature and case-conditioned development results

SF12 is archived as `old12`; CS is archived as `case_R`. All gains are percentages of squared-error reduction relative to the same task-specific M0; negative values mean increased error. Each task holds out all four combinations and trains on the other 17. Results are reported to six decimal places, with full precision retained in `source_data/v22_qualified_rectangles/` and the archived `run/SCORES.json`.

**A. Task identities.** Each row denotes the four combinations in the Cartesian product. Q4 is the previously inspected positive block.

| Task | Withheld rectangle |
|---|---|
| Q1 | {FERD3L, FLI1} × {MSGN1, NR5A2} |
| Q2 | {FERD3L, FLI1} × {MSGN1, PTF1A} |
| Q3 | {FERD3L, FLI1} × {NR5A2, PTF1A} |
| Q4 | {FERD3L, MSGN1} × {NR5A2, PTF1A} |
| Q5 | {FLI1, MSGN1} × {NR5A2, PTF1A} |

**B. All task-specific comparisons.** Complete gains sum the four query losses within each task. Partner gains use the unit four-corner contrast.

| Task | SF12 complete % | SF12 partner % | CS complete % | CS partner % |
|---|---:|---:|---:|---:|
| Q1 | +0.833319 | +1.193447 | -30.785255 | +0.905671 |
| Q2 | +3.726540 | -5.785755 | -47.928293 | -5.642622 |
| Q3 | +1.730775 | -1.784271 | -434.531746 | -4.528013 |
| Q4 | +4.665223 | +6.398984 | +9.291157 | +6.333680 |
| Q5 | +0.195160 | +4.091145 | +3.349258 | +3.948164 |

**C. All five tasks, with denominators.** Full-response loss is balanced over eight unique query edges; partner SSE pools the five overlapping contrast evaluations. These are different aggregation units.

| Method | Edge-balanced full loss | Complete gain % | Pooled partner SSE | Partner gain % |
|---|---:|---:|---:|---:|
| M0 | 11.953162408 | +0.000000 | 31.398220201 | +0.000000 |
| SF12 | 11.676074405 | +2.318115 | 31.129608867 | +0.855499 |
| CS | 23.253568584 | -94.539050 | 31.374869988 | +0.074368 |

**D. Excluding the previously inspected Q4 block.** Edge-balanced full losses are recomputed for the retained task occurrences before forming ratios.

| Method | Edge-balanced full loss | Complete gain % | Pooled partner SSE | Partner gain % |
|---|---:|---:|---:|---:|
| M0 | 11.919499737 | +0.000000 | 28.937740449 | +0.000000 |
| SF12 | 11.662816946 | +2.153469 | 28.826574823 | +0.384154 |
| CS | 25.641260816 | -115.120277 | 29.070229140 | -0.457840 |

These five overlapping tasks share eight edges and span three linearly independent contrast directions; they are not independent replications. Leave-one-task SF12 partner gains range from −0.428934% to +1.761012%. Complete-loss balancing, rejected CS results and additional denominator checks are detailed in Supplementary Note 12.
