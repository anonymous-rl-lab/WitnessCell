#!/usr/bin/env python3
"""Independent audit of official six-metric scoring and the SOTA verdict."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_UNITS = {
    "Norman": 268,
    "Wessels": 187,
    "Schmidt": 87,
    "Replogle_exp6": 112,
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.package_root.resolve()
    total_rows = 0
    total_waterstein = 0
    for dataset, units in EXPECTED_UNITS.items():
        score = root / "results/formal_score/full" / dataset
        table = pd.read_csv(score / "official_six_metric_raw.csv")
        assert len(table) == units * 10
        assert np.isfinite(table.performance).all()
        assert set(table[table.DEG == 100].metric) == {
            "pearson_distance", "mse", "edistance", "sym_kldiv",
            "wasserstein", "common_degs",
        }
        assert set(table[table.DEG == 5000].metric) == {
            "pearson_distance", "mse", "edistance", "sym_kldiv",
        }
        total_rows += len(table)
        for seed in (1, 2, 3):
            audit = json.loads((score / f"wasserstein_audit_seed{seed}.json").read_text())
            assert all(row["converged"] for row in audit)
            total_waterstein += len(audit)

    aggregate = root / "results/formal_score/aggregate"
    verdict = json.loads((aggregate / "formal_verdict.json").read_text())
    assert verdict["status"] == "PASS_FORMAL_AGGREGATION"
    assert verdict["official_primary_top100_position"] == 1
    assert verdict["sota_official_primary"] is True
    assert all(
        item["status"] == "PASS" and item["exact_published_ranks"]
        for item in verdict["official_aggregation_reconstruction"]
    )
    board = pd.read_csv(aggregate / "formal_top100_leaderboard.csv")
    assert board.iloc[0].method == "WitnessCell"
    assert board.iloc[1].method == "scouter"
    assert board.iloc[2].method == "linearModel"
    bootstrap = json.loads((aggregate / "paired_cluster_bootstrap.json").read_text())
    assert bootstrap["status"] == "PASS_PAIRED_CLUSTER_BOOTSTRAP"
    assert bootstrap["replicates"] == 20000
    report = {
        "status": "PASS_OFFICIAL_FULL_METRIC_SOTA_AUDIT",
        "raw_metric_rows": total_rows,
        "wasserstein_converged": total_waterstein,
        "official_primary_position": verdict["official_primary_top100_position"],
        "official_primary_mean_rank": verdict["official_primary_top100_mean_rank"],
        "both_panels_robustness_position": verdict["both_panels_robustness_position"],
        "dataset_equal_weight_top100_position": verdict["dataset_equal_weight_top100_position"],
    }
    assert total_rows == 6540 and total_waterstein == 654
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
