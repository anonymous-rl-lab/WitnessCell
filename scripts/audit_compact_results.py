#!/usr/bin/env python3
"""Fail-closed structural audit for the compact frozen result release."""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_json(relative: str) -> dict:
    path = ROOT / relative
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AssertionError(f"expected JSON object: {relative}")
    return value


def csv_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        next(reader)
        return sum(1 for _ in reader)


def main() -> int:
    gate19 = load_json("reproduction/gate19_v14/FORMAL_VERDICT.json")
    assert gate19["status"] == "PASS_GATE19_FORMAL_MEAN_RESPONSE_RELEASE"
    assert gate19["official_split_count"] == 12
    assert gate19["active_splits"] == 8
    assert gate19["inactive_exact_fallback_splits"] == 4
    assert gate19["condition_seed_units"] == 654
    assert gate19["official_six_metric_status"] == "INHERITED_V13_NOT_RELABELED"

    gate15_root = ROOT / "reproduction/gate15_scperturbench/results/formal_score/full"
    expected = {"Norman": 2680, "Wessels": 1870, "Schmidt": 870, "Replogle_exp6": 1120}
    rows = {
        dataset: csv_rows(gate15_root / dataset / "official_six_metric_raw.csv")
        for dataset in expected
    }
    assert rows == expected and sum(rows.values()) == 6540
    gate15 = load_json(
        "reproduction/gate15_scperturbench/results/formal_score/aggregate/formal_verdict.json"
    )
    assert gate15["status"] == "PASS_FORMAL_AGGREGATION"
    assert gate15["official_primary_top100_position"] == 1
    assert gate15["dataset_equal_weight_top100_position"] == 2

    gate21 = load_json("reproduction/gate21_selective/results/formal_reveal/FORMAL_VERDICT.json")
    assert gate21["status"] == "PASS_FROZEN_SELECTIVE_PREDICTION_GATE"
    assert gate21["test_rows"] == 213 and gate21["test_pairs"] == 33
    assert gate21["pair_identity_overlap_with_calibration"] == 0

    gate22 = load_json("reproduction/gate22_metrics/results/formal_e22/FORMAL_VERDICT.json")
    assert gate22["status"] == "COMPLETE_EXPERIMENT22_FORMAL_EXECUTION"
    assert gate22["METRIC_VALIDITY"] == "PASS_PHASE_M_EXECUTION_AND_SEAL"
    assert gate22["PRED_LINEAR"] == "NOT_ADJUDICATED"
    assert gate22["GATE21_WMSE"] == "NOT_ADJUDICATED_FULL_213"
    phase_p = ROOT / "reproduction/gate22_metrics/results/formal_e22/P/per_operation_metrics.csv"
    assert csv_rows(phase_p) == 2616

    print(
        json.dumps(
            {
                "status": "PASS_COMPACT_FROZEN_RESULT_AUDIT",
                "gate15_metric_rows": 6540,
                "gate19_units": 654,
                "gate21_rows": 213,
                "gate22_operation_metric_rows": 2616,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
