# Additional archived development tables

These tables preserve v22 numbering and wording as source documentation, not v24 SI numbering. They are not additional independent confirmations. Original note references refer to the v22 supplement in the code archive.

## Archived Table S11 | Original-system integration, fixed-G controls

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

## Archived Table S12 | Six original split selections and realized increments

| Dataset / seed | Test pairs | Selected policy | Complete gain | Projected gain | Test-graph contrast dimension |
|---|---:|---|---:|---:|---:|
| Norman / 1 | 73 | Uniform shrinkage | +0.3228% | +0.6241% | 24 |
| Norman / 2 | 66 | Original v14 / zero | +0.0000% | +0.0000% | 16 |
| Norman / 3 | 69 | Uniform shrinkage | +0.6991% | +2.1213% | 22 |
| Wessels / 1 | 63 | Uniform shrinkage | +0.0514% | -0.0011% | 39 |
| Wessels / 2 | 54 | Half-agreement / amplitude-noise | +1.5144% | +0.4513% | 33 |
| Wessels / 3 | 64 | Amplitude-only | +0.3130% | +0.0641% | 45 |

Selection used training predictions only. A uniform coefficient of zero in Norman seed 2 equals the unextended core. Every inner selection log, endpoint-availability stratum and output array is retained in the original-system run package.

## Archived Table S18 | Calibration of the measurement reference improves native M0 response prediction

| Input estimator | Source libraries | Original-scale MSE | Error reduction versus raw (%) |
|---|---:|---:|---:|
| Raw cell-weighted mean | 9 | 0.001623698041 | 0.0000 |
| Direct fixed-composition standardization | 7 | 0.001476556934 | 9.0621 |
| Shared-offset calibration | 7 | 0.001436203065 | 11.5474 |
| Shared-offset calibration | 9 | 0.001413839285 | 12.9247 |

All arms predict the same seven-library equal-composition target. Mean squared error averages 37,292 non-target RNA coordinates, two held-out technical-half references and 21 eligible conditions, with equal condition weight. Five outer folds were drawn from the historical development pool of 25 combinations; each used 20 combinations and ten singles for fitting, with preprocessing and M0 selection repeated inside three inner folds. Seven-library calibration reduced error by 2.7330% relative to direct standardization under the same available-library information. The corresponding nine-library difference was 4.2476% and additionally involved two source libraries. Nine-library calibration improved 20/21 conditions and 5/5 folds versus raw inputs; neither conditions with shared endpoints nor overlapping folds or technical halves are independent biological replications. This table concerns complete-response prediction without a partner head.



## Archived Table S19 | Complete five-task shared-feature and case-conditioned development results

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