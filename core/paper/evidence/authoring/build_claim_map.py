from pathlib import Path
import csv,json,hashlib
root=Path(__file__).resolve().parents[1]
items=[
('M1','Method and workflow','Implementation identities','Note 1; Table S1','source_data/REPRODUCTION_INDEX.json','Separate configurations and reproduction scopes'),
('M2','Joint-intervention information','Equal-world risk lower bound 0.0625','Note 2.3','source_data/formal_bridge/theory/TARGET_RISK_LOWER_BOUND.txt','Specified finite mechanism subproblem'),
('M3','Joint-intervention information','Matched panels and enlarged-family boundary','Notes 2.2, 2.6; Table S4','source_data/formal_bridge/theory/THREE_GATE_DESIGN_ENUMERATION.json','Exact design geometry; not a biological-network claim'),
('M4','Joint-intervention information','92.8409% target-risk decrease','Note 2.4; Table S2','source_data/formal_bridge/scores/summary.csv','Equal-world target MSE; paired noise simulation'),
('M5','Joint-intervention information','Paired interval and GLS comparator','Note 2.4; Table S2','source_data/formal_bridge/scores/paired_comparisons.csv','Monte Carlo interval over 1000 seeds'),
('M6','Joint-intervention information','Candidate-family target coverage','Note 2.5; Table S3','source_data/formal_bridge/scores/world_summary.csv','Known-noise, declared mechanism family'),
('M7','Formal prediction performance','First of 16 eligible entries','Note 3.2; Table S5','source_data/publication/released_benchmark/formal_top100_leaderboard.csv','Released core; 654 operation-seed records'),
('M8','Formal prediction performance','Six-metric strengths and weaknesses','Note 3.2; Table S6','source_data/publication/released_benchmark/official_top100_per_metric_leaderboard.csv','Complete metric-specific positions'),
('M9','Formal prediction performance','Paired rank differences','Note 3.2; Table S7','source_data/publication/released_benchmark/paired_cluster_bootstrap.csv','384 dataset-perturbation clusters'),
('M10','Formal prediction performance','Combination WMSE decrease 24.87%','Note 3.3; Table S8','source_data/publication/Fig2_prediction_ratios.csv','540 combination-seed records; specified trainMean comparator'),
('M11','Formal prediction performance','Non-target error controls','Note 3.3; Table S8','source_data/publication/ED1_coverage_comparisons.csv','Full-cell means; non-target masks; geometric seed aggregation'),
('M12','Partner response','Released endpoint branch cancels under projection','Note 4; Table S9','source_data/publication/Step2_Fig5_Official12.csv','Structural identity and retained split scores'),
('M13','Partner response','RNA R rectangle gain 15.9622%','Note 5; Table S10','source_data/publication/Step2_Fig5_Matched55.csv','Developmental matched panel; mean seed losses for neural models'),
('M14','Partner response','All RNA folds','Note 5.2; Table S11','source_data/publication/Step2_ED5_Folds.csv','Fixed-model seed copies counted once'),
('M15','Partner response','All eleven same-fit rectangles','Note 5.2; Table S12','source_data/publication/RNA_rectangle_definition.csv','One-based source folds mapped to zero-based display'),
('M16','Partner response','Reliability and uniform-amplitude controls','Note 6; Table S13','source_data/publication/RNA_supplement_controls.json','Conditional controls; no new independent replication'),
('M17','Partner response','Scale dependence and ten-readout concentration','Note 6.3; Table S13','source_data/publication/RNA_control_scale.json','Raw versus standardized errors retained'),
('M18','Functional growth','K562 budget comparison and Jurkat results','Note 7; Tables S14-S15','source_data/publication/Horlbeck_figure_metrics.csv','Revised K562 development; Jurkat procedure transfer'),
('M19','Functional growth','Training/test populations','Note 7.1; Table S14','source_data/publication/Horlbeck_design.csv','Protected original test pairs and current folds'),
('M20','Functional growth','Jurkat conditional intervals','Note 7.3; Table S15','source_data/publication/Horlbeck_Jurkat_intervals.csv','Conditional reweighting of fixed predictions'),
('M21','Functional growth','Jurkat all-pair evidence','Note 7.3','source_data/publication/biological_examples/jurkat/Jurkat_all_617_pairs_two_repeats.csv','617 pairs and two source repeats'),
('M22','Risk selection','Analytic acceptance and complete-case follow-up','Note 8.3; Table S16','source_data/legacy_tables_v22.md','Original Table S7; 213 MSE / 207 WMSE records'),
('M23','Risk selection','Support versus exposure ranking','Note 8.4; Table S16','source_data/publication/Fig4_risk_comparisons.csv','Different ranking implementation from analytic risk'),
('M24','PairPick','Complete 15-row result population','Note 9.4; Table S17','source_data/gse220974/SUMMARY.csv','Primary five-model set distinguished from historical outputs'),
('M25','PairPick','All 420 choices and empirical-best hits','Note 9.4','source_data/gse220974/CHOICES.csv','15 models by 28 tasks'),
('M26','PairPick','Paired utility uncertainty','Note 9.4; Table S18','source_data/gse220974/PAIRED_COMPARISONS.csv','Eight physical-anchor clusters'),
('M27','PairPick','All candidate pools','Note 9.3; Table S19','source_data/gse220974/GROUPS.json','28 tasks; five pools with at least four candidates'),
('M28','PairPick','Thirty neural jobs completed','Note 9.5; Table S20','source_data/gse220974/COMPLETION.json','Selection/refit JSON for every job in source_data/gse220974/runs'),
('M29','Source-transfer limits','Frozen Dixit complete comparison','Note 10.1; Table S21A','source_data/legacy_tables_v22.md','Original Tables S13-S14; later program development separately labelled'),
('M30','Source-transfer limits','Later Dixit program development','Note 10.1; Table S21B','source_data/publication/Dixit_fixed_reference/WitnessCell_Dixit_FixedReference_20261020_Key_Results.csv','Different endpoint and comparator from frozen RNA test'),
('M31','Source-transfer limits','Frozen Combi-seq result','Note 10.2; Table S21C-D','source_data/publication/CombiSeq_confirmation/WitnessCell_Frozen_Program_Confirmation_20261020_KeyResults.csv','Principal assembly versus its endpoint baseline'),
('M32','Acquisition boundary','Development target-aligned acquisition','Note 10.3; Table S22','source_data/acquisition/development/summary.json','Expected uniform acquisition comparator; 5% gate unmet'),
('M33','Acquisition boundary','Outer acquisition fails to retain gain','Note 10.3; Table S22','source_data/acquisition/outer/summary.json','All six outer targets retained'),
('M34','Evidence history','Original outcomes and subsequent amendments','Note 11.2; Table S23','source_data/legacy_tables_v22.md','Original Table S6; local locks distinct from external deposit'),
('M35','Reproduction','Exact declared task scopes','Note 12; Table S24','source_data/REPRODUCTION_INDEX.json','Code release 27; legacy publication renderer v22')]
for row in items:assert (root/row[4]).is_file(),row
with (root/'CLAIM_SOURCE_MAP.csv').open('w',newline='') as h:
 w=csv.writer(h);w.writerow(['claim_id','main_section','claim_or_endpoint','SI_location','source_path','scope']);w.writerows(items)
print(json.dumps({'claim_rows':len(items),'all_source_paths_exist':True}))
