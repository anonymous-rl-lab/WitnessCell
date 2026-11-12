# Experiment 22 — metric-calibration stress test

This directory contains the formally frozen v2 implementation of Experiment 22.
Its protocol, code, candidate-blind metric contracts, comparator scope and
Gate21 identity contract are sealed by `FROZEN_MANIFEST.sha256` before any new
formal WMSE/weighted-delta-R² result is computed.

## Why v2 exists

The superseded draft mixed metric sources and contained four material errors: it used a noncanonical weight formula, replaced weighted delta R² with weighted Pearson, misdefined NIR, and invented a DRF cutoff. It also allowed formal-data pilots and proposed refitting the already frozen Gate 21 threshold.

v2 fixes those defects and separates the experiment into:

1. candidate-independent metric validity;
2. prediction robustness under WMSE and weighted delta R², with NIR reported separately;
3. decision robustness for the immutable Gate 19 and Gate 21 decisions.

## Files

- `PROTOCOL_v2.md` — complete scientific and execution protocol;
- `SOURCE_LOCK.json` — immutable repositories, commits, files, hashes and semantic boundaries;
- `PRE_FREEZE_CHECKLIST.md` — fail-closed blockers that must pass before hashing a scientific freeze;
- `SUPERSESSION_NOTICE.md` — identity and correction ledger for the non-executable v1 draft;
- `VERSION` — protocol-package version;
- `FROZEN_MANIFEST.sha256` — scientific pre-execution freeze;
- `FREEZE_RECEIPT.json` — hash and file count of that manifest (excluded from
  the manifest to avoid self-reference).

## Formal execution

The once-only runner is:

```bash
bash run_formal_experiment22.sh
```

It verifies the frozen manifest, runs candidate-blind Phase M under a runtime
read firewall, seals those outputs, then opens immutable v14/v13 predictions
for Phase P and fixed Gate21 decisions for Phase D. It refuses to overwrite an
existing formal result directory.

The comparator scope is frozen as `MANDATORY_BASELINES_ONLY`; exact raw
`linearModel` condition means are unavailable, so `PRED_LINEAR` is
`NOT_ADJUDICATED`. The cell-level Norman asset lacks ELMSAN1, leaving six
Gate21 query instances without canonical weights. The full 213-query WMSE
verdict is therefore also `NOT_ADJUDICATED`; a 207-query complete-case analysis
is explicitly labelled sensitivity evidence and cannot replace it.
