# WitnessCell code release 28.0.0 — paper v36

This release separates the fixed v36 manuscript from its scientific reproduction materials. See `README_CN.md` for the full guide, `REPRODUCTION_INDEX.json` for nine offline tasks and `docs/PAPER_COVERAGE.csv` for all five main figures, five supplementary figures and 52 SI tables.

```bash
python -B tools/verify.py
python -B tools/reproduce.py --out ../v36_reproduction_new
```

Use a fresh output directory outside this release. Tested dependency versions are in `requirements-tested.txt`; no network access is needed after installing them. The default workflow does not retrain biological predictors. It replays archived scores/predictions, performs Han frozen inference and target scalar calibration, checks interfaces and redraws figures. `qa/VALIDATION_SUMMARY.json` records actual execution scopes.

Exact stored bytes, exact array reconstruction and precision-aware frozen-model replay are separate guarantees. Read `docs/REPRODUCTION_LEVELS.md`. Large raw sequencing and external competitor assets remain referenced, rather than silently claimed to be bundled. Original rights restrictions are unchanged.
