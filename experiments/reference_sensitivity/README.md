# Reference sensitivity and direction decomposition

Integrated from WitnessCell_v41_Revision_Supplement.zip, language revision R1. Scientific numerical code retains its v40 source names because v41 changed text and figure presentation, not estimation rules.

From the repository root run `python -B reproduce.py --out ../reference-replay --tasks reference_sensitivity`. This replays saved predictor outputs, checks four saturation-reference configurations, recomputes Han direction decomposition, and checks nine full-precision tables. Use `--tasks saturation_refit` to actually refit the four configurations first. Default replay restores necessary Han caches from saved weights without query labels.

Canonical paper text is under paper/manuscript_sources. All original companion files are mapped in provenance/V41_INTEGRATION_MAP.csv.
