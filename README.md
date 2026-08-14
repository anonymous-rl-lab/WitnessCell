# WitnessCell

Anonymous, reviewer-ready source and reproducibility repository for WitnessCell.

This repository combines two release layers:

1. **`witnesscell 0.1.0`** — the typed, tested Python API and CLI distributed
   under Apache-2.0;
2. **v18 reproducibility assets** — the frozen v14 mean-response algorithm,
   official four-dataset benchmark tables, selective-prediction analysis and
   metric-calibration stress test.

The repository is deliberately smaller than the original 430 MiB evidence
archive. It keeps executable code, compact frozen results, condition-moment
inputs and integrity records. Raw public single-cell matrices, historical
duplicate predictions, nested archives, caches and machine-specific logs are
excluded.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-repro-v14.txt
python -m pip install -e ".[dev]" --no-deps
bash scripts/run_smoke.sh
```

The smoke suite tests the public API, frozen numerical contract, official
benchmark tables, selective-prediction result and metric primitives. It does
not download raw `.h5ad` files.

## Reproduce the v14 mean-response predictions

The four condition-moment caches and required GEARS gene assets are included.
Run one split:

```bash
DATASETS=Norman SEEDS=1 bash scripts/reproduce_v14.sh
```

Run all four datasets and three seeds:

```bash
bash scripts/reproduce_v14.sh
```

The runner produces only the target-free deployment archives by default and
checks their array-level semantic digests against the frozen v18 references.

## Full cell-level benchmark

The official public matrices total about 0.9 GiB compressed and are not stored
in Git. Download and verify them with:

```bash
bash scripts/download_full_data.sh
```

Full distributional scoring and Experiment 22 use separate, pinned metric
dependencies. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) and
[DATA.md](DATA.md).

## Repository map

| Path | Purpose |
|---|---|
| `src/witnesscell/` | Hardened v14 package implementation |
| `tests/` | Unit, contract, serialization and security tests |
| `reproduction/gate19_v14/` | Frozen v14 research implementation and result tables |
| `reproduction/gate15_scperturbench/` | Four-dataset benchmark code and official frozen scores |
| `reproduction/gate21_selective/` | Frozen accept/abstain analysis |
| `reproduction/gate22_metrics/` | Metric-calibration protocol, code, tests and compact results |
| `reproduction/theory/` | Geometry, exact-risk and estimated-risk verification |
| `reproduction/named_model_gate/` | Compact named-model audit evidence |
| `reproduction/assets/` | Included derived inputs and reference checksums |

## Scope

WitnessCell predicts condition means for unmeasured two-endpoint genetic
combinations. It does not claim cell-level generative simulation, universal
uncertainty calibration, causal identification, or automatic routing to a
different model. Selective prediction has two actions: `accept` or `abstain`.

The repository preserves the frozen `NOT_ADJUDICATED` outcomes. Missing raw
linear-model predictions and incomplete weighted comparisons are not converted
into positive claims.

## Double-blind release

Author, institution, account, DOI and contact metadata for this project are
withheld. Third-party dataset and software links remain because they identify
external sources, not the submission authors. See
[docs/GITHUB_UPLOAD.md](docs/GITHUB_UPLOAD.md) before publishing.

## License

The package code is Apache-2.0. External datasets, dependencies and upstream
benchmark implementations retain their own terms. See `LICENSE` and
`THIRD_PARTY_NOTICES.md`.
