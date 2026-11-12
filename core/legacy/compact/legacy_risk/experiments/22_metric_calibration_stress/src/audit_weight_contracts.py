#!/usr/bin/env python3
"""Audit formal weight contracts against the source-parity implementation path."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


DATASETS = ("Norman", "Wessels", "Schmidt", "Replogle_exp6")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment
    generator = (root / "src/phase_m_metric_validity.py").read_text()
    rows = []
    for dataset in DATASETS:
        path = root / "assets/formal_contracts" / f"{dataset}.phase_m_contract.npz"
        with np.load(path, allow_pickle=False) as data:
            weights = data["weights"].astype(np.float64)
            evaluable = data["weight_evaluable"].astype(bool)
            valid = weights[evaluable]
            checks = {
                "finite": bool(np.all(np.isfinite(valid))),
                "nonnegative": bool(np.all(valid >= 0)),
                "positive_sum": bool(np.all(valid.sum(axis=1) > 0)),
                "source_minmax_max_one": bool(np.allclose(valid.max(axis=1), 1.0, rtol=0, atol=1e-14)),
                "has_zero_after_alignment": bool(np.all(np.any(valid == 0, axis=1))),
            }
            rows.append(
                {
                    "dataset": dataset,
                    "conditions": int(weights.shape[0]),
                    "evaluable": int(evaluable.sum()),
                    "contract_sha256": sha256(path),
                    "checks": checks,
                    "pass": all(checks.values()),
                }
            )
    global_checks = {
        "generator_calls_source_parity_wrapper": "source_weight_transform(" in generator,
        "generator_uses_scanpy_locked_test": "method=\"t-test_overestim_var\"" in generator,
        "generator_uses_adjusted_second_half_pvalues": "result[\"pvals_adj\"]" in generator,
    }
    report = {
        "status": "PASS_FORMAL_WEIGHT_CONTRACT_AUDIT"
        if all(row["pass"] for row in rows) and all(global_checks.values())
        else "FAIL_FORMAL_WEIGHT_CONTRACT_AUDIT",
        "datasets": rows,
        "global_checks": global_checks,
        "source_parity_unit_test": "tests/test_source_parity.py::test_weight_transform",
        "pass": all(row["pass"] for row in rows) and all(global_checks.values()),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
