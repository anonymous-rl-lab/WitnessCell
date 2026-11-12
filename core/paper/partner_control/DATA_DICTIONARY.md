# Numerical fields and units

`full_sse`, `endpoint_sse`, `partner_sse` use the same unweighted response axis and sum across the stated pairs and readouts/repeats. Their matching `base_*` columns are M0 losses. Gain is 100 × (1 − loss/base); positive means lower error. A difference between two gains is in percentage points.

`alpha` is the bounded coefficient calculated from internal folds. `rho` is internal projected residual/correction alignment. `k_original` is correction energy divided by the cross-product when the cross-product is positive; it is undefined for zero or adverse alignment. `q_raw` and `q_oof` are historical raw and baseline-residual training repeat agreements; the latter depends on the baseline. Missing RNA residual agreement is retained as missing.

The functional diagnostic table joins these quantities to the actual selected zero/rank-one predictor and all four outputs. The RNA table records `model` as the candidate-workflow label and `actual_family` / `chosen_eta` when the source selection record supplies them. Archive aliases and later snapshots retain their identities.

No absolute SSE is pooled between Horlbeck and Shen. Within a condition, SSE is summed before forming the ratio. K562 layouts and RNA candidate records share queries. Summary win/tie/loss counts describe records or configurations, not independent biological trials.
