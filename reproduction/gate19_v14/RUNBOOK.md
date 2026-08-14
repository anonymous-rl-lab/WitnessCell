# Gate 19 regeneration runbook

Run from the repository root. Gate 19 reuses Gate 18 condition-moment caches and
Gate 15 GEARS assets; no raw h5ad is required for mean-response regeneration.

```bash
for dataset in Norman Replogle_exp6 Schmidt Wessels; do
  for seed in 1 2 3; do
    PYTHONPATH=experiments/19_v14_incremental_amplitude_gate/src \
    python experiments/19_v14_incremental_amplitude_gate/src/run_v14_from_cache.py \
      --cache experiments/18_dual_head_evidence_gate/cache/${dataset}.condition_moments.npz \
      --dataset "$dataset" --seed "$seed" \
      --gears-assets experiments/15_scperturbench_sota/module/data/gears_assets \
      --out /tmp/v14/${dataset}/seed${seed}
  done
done
```

Use `src/compare_v14_v13.py` against Gate 18 deploy archives for the 654-unit
paired audit. Use Gate 15 `score_directional.py` for the exact-mean directional
comparison. Raw `predictions.npz` contains held-out means for internal scoring
and must not be shipped; publish only `deploy_predictions.npz`.

Run `python experiments/19_v14_incremental_amplitude_gate/audit.py` for the
release audit. Full distributional scoring additionally requires the four raw
h5ad files documented by Gate 15 and is outside the packaged Gate 19 claim.
