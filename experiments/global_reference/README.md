# Matched global-reference analysis (paper v39)

Integrated at the root of code release v31. The scientific scripts, permitted inputs, original and newly fitted predictions, parameter grids, and results are preserved from the v39 supplement.

From the release root:

```bash
python -B reproduce.py --out ../v39_global_replay --tasks global_reference
python -B reproduce.py --out ../v39_global_refit --tasks global_refit
```

The first command uses saved matrix predictions and recomputes scalar/reference diagnostics. The second refits all seven configurations before scoring. Both compare the resulting full-precision tables and displayed S54–S57 rows. Default replay is not predictor retraining. Canonical manuscript sources and vector references are under paper/; source identities are recorded in provenance/V39_INTEGRATION_MAP.csv and the current normalization guide.
