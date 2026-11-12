# Experiment 28 — RNA / Jurkat response examples

The original standalone response-display archive is embedded byte for byte, once.
It contains four-corner RNA observations and saved predictions, all 54,945 signed
rectangle/readout entries, both measured repeats and predictions for 617 Jurkat
pairs, selected-case rules, figures, prose, and self-contained checking/drawing code.

From the repository root, using existing Python >=3.9, NumPy and Matplotlib >=3.5. Archive verification uses only the Python standard library:

```bash
python tools/run_response_examples_v23.py --mode verify
python tools/run_response_examples_v23.py --mode reproduce --out /tmp/witnesscell_v23_examples_run1
```

Use a new external output directory. Reproduction checks the included numerical
tables and redraws all four example figures; it does not refit models, install
software, download data, or rerun the archived original-data extraction scripts.
The source archive is left unchanged. Original-cell reconstruction and neural
retraining require their separately retained source assets.

Cases are descriptive, selected after evaluation. RNA retains original v17
predictions and contrast orientations; named Jurkat examples retain the frozen
primary rank-one predictor. They support explanation of the completed results,
not an additional prospective experiment. The display source was integrated into
paper v18; its original README's paper-insertion wording documents its creation.

Integration compatibility: the original figure script expects row/column endpoint
labels that are absent from its distributed four-corner CSV. The v23 drawing
adapter derives these labels from the existing AB/AD/CB/CD pair identities,
checks the rectangle and contrast signs, and supplies display labels in memory.
It does not alter the source ZIP, observations, predictions or contrasts.
