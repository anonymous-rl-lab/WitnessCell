#!/usr/bin/env python3
"""Independent audit of the 4-dataset x 3-seed deployment predictions."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


EXPECTED = {
    "Norman": {1: 93, 2: 86, 3: 89},
    "Wessels": {1: 65, 2: 55, 3: 67},
    "Schmidt": {1: 31, 2: 28, 3: 28},
    "Replogle_exp6": {1: 42, 2: 36, 3: 34},
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.package_root.resolve()
    rows = []
    for dataset, seeds in EXPECTED.items():
        split_audit = json.loads(
            (root / "results/formal_combo" / dataset / "audit/split_audit.json").read_text()
        )
        assert split_audit["status"] == "PASS"
        assert all(row["exact_match"] for row in split_audit["rows"])
        for seed, expected_conditions in seeds.items():
            run = root / "results/formal_combo" / dataset / f"seed{seed}"
            manifest = json.loads((run / "manifest.json").read_text())
            archive = np.load(run / "deploy_predictions.npz", allow_pickle=False)
            assert manifest["status"] == "PASS_WITNESSCELL_MEAN_PREDICTION"
            assert manifest["split_counts"]["test"] == expected_conditions
            assert str(archive["dataset"]) == dataset
            assert int(archive["seed"]) == seed
            assert len(archive["conditions"]) == expected_conditions
            assert archive["prediction"].shape[0] == expected_conditions
            assert "truth" not in archive.files and "truth_variance" not in archive.files
            assert all(archive[key].dtype != object for key in archive.files)
            rows.append({
                "dataset": dataset,
                "seed": seed,
                "test_conditions": expected_conditions,
                "training_doubles": manifest["training_doubles"],
                "gamma": manifest["selected"]["gamma"],
            })
    report = {
        "status": "PASS_FOUR_DATASET_CPU_PREDICTION_AUDIT",
        "fits": len(rows),
        "condition_seed_units": sum(row["test_conditions"] for row in rows),
        "target_free_deployment_archives": len(rows),
        "rows": rows,
    }
    assert report["fits"] == 12
    assert report["condition_seed_units"] == 654
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
