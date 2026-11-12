# GSE220974: new-data comparison, 2026-11-02

The scientific target is prediction and partner selection among measured, held-out double-guide identities in the same K562 experiment. A choice means best within its measured candidate pool for a specified single-intervention RNA program, not globally optimal treatment or a causal claim about biological fitness.

## Cohort and limitation

Primary operational cohort: ALL 49 detected distinct-regular-gene, mode-specific double identities in the supplied PRISM GSE220974 file. No expression-dependent condition filtering or minimum-count threshold chosen to reproduce the publication's count. These comprise 28 base-gene pairs. Enhancer and same-gene double identities are outside this comparison. True NTC is determined by guide identity; 860 enhancer-target cells formerly labelled control are excluded from NTC.

**Unresolved provenance limitation:** the original study reports 22 designed /21 detected double constructs. The published supplement was located, but its Excel data could not yet be downloaded through the available endpoints. Therefore the 49 detected labels are NOT asserted to be 49 verified intended constructs, and results are NOT a replication of the publication's designed-construct cohort. Incidental co-delivery, mixed guide assignment or doublets cannot be excluded by this comparison. The completed run will retain this limitation. Any later design whitelist is a separately labelled analysis and will not silently replace these results.

## Common data and boundary

- Raw counts: log1p(10000 × counts / cell library size).
- Subtract capture-matched true-NTC mean, add the equal D1/D2/D3 NTC reference; retain continuous, potentially negative corrected values.
- Condition means weight observed cells equally. A rare one-cell condition/capture stratum is retained once and belongs to only one half; it is never duplicated to manufacture reproducibility.
- Select 5000 variance-ranked genes using ONLY true NTC and single cells; include all physical target genes. Result: 5006 readouts, 4987 non-target evaluation genes, 19 physical target RNA ports and 35 distinct gene/mode intervention identities.
- Existing single-derived TargetRank objectives: signed single effect × nonnegative half agreement, L1 normalized; physical targets excluded. No double outcome enters objective or readout selection.
- Within condition/capture, deterministic complementary cell halves (seed 20261102); each half subtracts its own control-half mean. These are sampling checks, not additional biological replicates. Capture-specific results will be reported without calling captures independent biological replicates.
- Five outer folds: hash-ordered balanced base-gene-pair blocks. All activation/repression variants of a base pair stay together. Three inner folds and three nuisance cross-fit blocks use the same identity-only group rule. No outer query response enters fitting, epoch choice, graph inference or feature selection.

## Models

Original roster retained: R-G, PairPick-Stable, PairPick-TargetRank, SingleB-additive, standard ridge selected by RNA MSE or choice utility; GEARS seeds0/1/2 and their mean; Scouter seeds0/1/2 and their mean, Scouter fixed40 seed24; PRESAGE Combo seed24. TargetRank is the archived round9 parallel-head algorithm, with its exact grids, not a newly selected variant. R-G preserves its archived RNA/partner guard and fallback.

Mode adapters are explicit, not labelled untouched official mixed-mode benchmarks:

- R-G/ridge/TargetRank use separate gene-mode incidence columns.
- Graph relations are fitted over 19 physical RNA genes, with 35 mode-specific single interventions. Any intervention of either mode masks that physical gene's equation. Mode-specific single responses supply clamps. Expansion to ranking coordinates occurs after the graph solve.
- Scouter uses the original model/loss, with direct-sum GenePT vectors: activation [embedding,0], repression [0,embedding]. Thus swapping activation and repression between two genes changes the input.
- PRESAGE uses the archived CellSimBench Combo class and all 40 published priors plus two training-only PCA sources; input feature dimension is doubled into disjoint activation/repression blocks. Architecture widths, loss, optimizer and original priors remain as archived.
- GEARS uses the original model/loss and two disjoint mode-labelled copies of its GO perturbation graph. The RNA coexpression graph has physical genes only and is computed from controls/singles.

## Training budget and selection

Classical grids are unchanged and use all three inner folds. Neural epoch selection uses the same inner fold0 conditions for all architectures. GEARS/Scouter: mean original four-edge contrast MSE, non-target RNA MSE for exact ties, then earlier epoch. This new cohort's five inner0 sets each contain at least one complete rectangle. No missing-rectangle fallback is required.

GEARS: archived two schedules, LR halving every1 epoch (20 cap) or every5 (60 cap), patience10; batch32, Adam lr0.001/weight_decay0.0005; fresh refit using chosen schedule/epoch. Scouter: archived 40 epochs, optional one-time 40-epoch training-side extension when the best is in the last5 and last5 validation improves >1% over preceding epochs; batch256, Adam0.001, exponential LR0.9; fresh refit. Fixed40 runs seed24 for exactly40 epochs. PRESAGE: Adam0.0024, batch32, up to1000 epochs with inner MSE patience10; refit selected epoch. Each architecture retains its archived loss; no best-seed reporting.

Smoke tests use only training inputs and validate gradients, mode distinctions, output shape, finite values, and memory/time. They do not change hyperparameters based on outer outcomes.

## Scoring after prediction freeze

Same fixed-A, same outer-fold candidate pools for every method. Primary: mean achieved objective utility and regret to the best measured candidate in each pool; top-choice accuracy and paired per-task differences. Report pools of size≥2 and the predefined size≥4 subgroup separately. Describe maxima as best observed in this pool, never global optima.

Secondary: non-target RNA MSE, endpoint-eliminated residual error within each SAME fitted fold (with rank reported; zero-rank folds are undefined, not zero error), half/capture sensitivity. TargetRank/Stable ranking-only heads do not acquire fictional RNA predictions. No success threshold at 5% or10%; report actual magnitude, consistency and uncertainty, including null or worse results. No post-score model tuning or favorable subgroup replacement.
