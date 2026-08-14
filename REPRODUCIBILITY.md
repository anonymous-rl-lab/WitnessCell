# Reproducibility guide

## Reproduction tiers

| Tier | Command | Inputs | Typical scope |
|---|---|---|---|
| 0 — API tests | `python -m pytest -q` | none | 67 package tests |
| 1 — repository smoke | `bash scripts/run_smoke.sh` | included compact assets | contracts and frozen-result audits |
| 2 — v14 mean response | `bash scripts/reproduce_v14.sh` | included condition moments | 4 datasets × 3 seeds, CPU |
| 3 — cell-level scoring | `bash scripts/download_full_data.sh` then Gate 15/22 runbooks | public h5ad files | official distributional metrics |

Tier 0 and Tier 1 are the continuous-integration gates. Tier 2 is the primary
reviewer reproduction path. Tier 3 is intentionally not run in CI because it
downloads about 0.9 GiB compressed data and invokes optional metric packages.

## Environment A: v14 prediction

Reference platform:

- Python 3.12;
- CPU only; no CUDA or PyTorch is required;
- pinned packages in `requirements-repro-v14.txt`.

Create the environment:

```bash
python3.12 -m venv .venv-v14
source .venv-v14/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-repro-v14.txt
python -m pip install -e . --no-deps
```

Run one split before the full matrix:

```bash
DATASETS=Norman SEEDS=1 OUTPUT_DIR=/tmp/witnesscell-v14-smoke \
  bash scripts/reproduce_v14.sh
```

Then run the frozen matrix:

```bash
OUTPUT_DIR=/tmp/witnesscell-v14-formal bash scripts/reproduce_v14.sh
```

The formal matrix is four datasets (`Norman`, `Replogle_exp6`, `Schmidt`,
`Wessels`) by seeds `1, 2, 3`. The reference CPU execution completed in about
17 minutes; runtime varies with processor and BLAS configuration.

### Equality rule

NPZ containers embed ZIP metadata, so byte hashes can vary despite identical
arrays. Formal regeneration is therefore checked with a canonical semantic
digest over sorted array names, dtypes, shapes and C-order bytes. The original
container hashes remain provenance records but are not used as numerical
equality tests.

## Environment B: metric calibration

Experiment 22 was frozen under Python 3.12 with the versions recorded in
`requirements-metrics-v2.txt`. Keep this environment separate from the v14
environment because it uses a different NumPy/anndata/scanpy stack.

```bash
python3.12 -m venv .venv-metrics
source .venv-metrics/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-metrics-v2.txt
```

Metric primitive tests:

```bash
PYTHONPATH=reproduction/gate22_metrics/src \
python -m unittest discover \
  -s reproduction/gate22_metrics/tests \
  -p 'test_metric_core.py'
```

The complete once-only protocol and phase boundaries are documented in
`reproduction/gate22_metrics/PROTOCOL_v2.md`. The four large Phase-M contracts
are excluded from Git because they can be regenerated from the verified public
h5ad files; their original hashes are retained in
`reproduction/assets/checksums/EXCLUDED_LARGE_ASSETS.sha256`.

## Frozen decisions and boundaries

- v14 correction is active in 8/12 splits; four splits exactly fall back to
  v13.
- Gate 21 is a frozen accept/abstain result, not a general uncertainty
  guarantee.
- Gate 22 keeps `PRED_LINEAR=NOT_ADJUDICATED` and
  `GATE21_WMSE=NOT_ADJUDICATED_FULL_213`.
- Deployment archives never contain held-out target means or variances.

## Failure policy

Stop if any checksum, split identity, target-free archive check, expected row
count or frozen verdict changes. Do not silently recompute a threshold, replace
a missing comparator or relabel an inherited result.
