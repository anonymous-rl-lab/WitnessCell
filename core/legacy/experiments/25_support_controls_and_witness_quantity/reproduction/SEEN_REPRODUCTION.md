# Portable seen-analysis reproduction

From a relocated copy of the delivered `step1_followup_20261017` directory, run:

```bash
python reproduction/reproduce_seen.py --output /path/to/a/new/seen_reproduction
```

An explicit package location is also supported:

```bash
python reproduce_seen.py --package-root /path/to/step1_followup_20261017 --output /path/to/a/new/seen_reproduction
```

Requirements are Python 3.9+ and NumPy. Use the same NumPy environment as the
original run when available. No raw cell matrices or original sandbox directory
are required.

The wrapper restores the seven retained `original_*` inputs for each dataset and
seed into a temporary source directory. It copies the frozen `run_seen.py` and
`analysis_plan.json` without editing either, then runs `prepare` and `evaluate`
sequentially with an explicit `--source` path. These commands reproduce the same
C/C0 risk-head fitting, frozen S/B decisions and bootstrap statistics already
specified in the delivered analysis. They do not fit a new mean predictor, add
methods or change the analysis choices.

The output directory must be new. It contains the copied runner and plan, the
reproduced `seen/` outputs, separate `prepare.log` and `evaluate.log` files, and
`reproduction_status.json`. The temporary source directory is automatically
removed. The delivered frozen directory remains untouched.

The machine-readable comparison requires the same 480 dataset/scope/coverage/
loss/comparison groups and both effect-quartile controls. Every statistical field,
all 962 confidence intervals, seed ratios, candidate counts, quotas and valid/
invalid bootstrap counts are compared. Strings and integer counts must match
exactly; floating-point values use absolute and relative tolerances of `1e-12`.
Completion times, elapsed time and the newly generated freeze hash are excluded
from comparison. A matching reproduction prints JSON with `"status": "PASS"`
and exits with code 0; a failure prints `"status": "FAIL"` and exits nonzero.

This checks reproduction of the delivered retrospective analysis. It does not
provide a new experiment or independent confirmation of its scientific claims.
