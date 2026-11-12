#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


EXPECTED_FALLBACKS = {
    ("Schmidt", 3),
    ("Wessels", 1),
    ("Wessels", 2),
    ("Wessels", 3),
}


def close(actual: float, expected: float, tolerance: float = 1e-12) -> None:
    if not np.isfinite(actual) or abs(actual - expected) > tolerance:
        raise AssertionError((actual, expected, tolerance))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--regenerated", type=Path, required=True)
    parser.add_argument("--paired", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    gate19 = args.repo / "experiments/19_v14_incremental_amplitude_gate"
    frozen_root = gate19 / "predictions/formal_combo"
    generated_paths = sorted(args.regenerated.glob("*/seed*/deploy_predictions.npz"))
    assert len(generated_paths) == 12, len(generated_paths)

    array_comparisons = 0
    global_max_abs_delta = 0.0
    active = set()
    inactive = set()

    for generated_path in generated_paths:
        dataset = generated_path.parent.parent.name
        seed = int(generated_path.parent.name.removeprefix("seed"))
        key = (dataset, seed)
        frozen_path = frozen_root / dataset / f"seed{seed}" / "deploy_predictions.npz"
        assert frozen_path.exists(), frozen_path

        with np.load(generated_path, allow_pickle=False) as generated, np.load(
            frozen_path, allow_pickle=False
        ) as frozen:
            assert len(generated.files) == 13
            assert set(generated.files) == set(frozen.files)
            assert "truth" not in generated.files and "truth_variance" not in generated.files
            for name in generated.files:
                assert generated[name].shape == frozen[name].shape
                assert generated[name].dtype == frozen[name].dtype
                assert np.array_equal(generated[name], frozen[name]), (key, name)
                if generated[name].dtype.kind not in "OUS" and generated[name].size:
                    delta = float(np.max(np.abs(generated[name] - frozen[name])))
                    global_max_abs_delta = max(global_max_abs_delta, delta)
                array_comparisons += 1

        manifest = json.loads((generated_path.parent / "manifest.json").read_text())
        if bool(manifest["identity_head"]["correction_active"]):
            active.add(key)
        else:
            inactive.add(key)

    assert array_comparisons == 156
    assert global_max_abs_delta == 0.0
    assert len(active) == 8
    assert inactive == EXPECTED_FALLBACKS

    v13_root = args.repo / "experiments/18_dual_head_evidence_gate/predictions/formal_combo"
    for dataset, seed in sorted(inactive):
        generated_path = args.regenerated / dataset / f"seed{seed}" / "deploy_predictions.npz"
        v13_path = v13_root / dataset / f"seed{seed}" / "deploy_predictions.npz"
        with np.load(generated_path, allow_pickle=False) as generated, np.load(
            v13_path, allow_pickle=False
        ) as v13:
            assert np.array_equal(generated["conditions"], v13["conditions"])
            assert np.array_equal(generated["prediction"], v13["prediction"])

    paired = json.loads(args.paired.read_text())
    assert paired["status"] == "PASS_V14_V13_PAIRED_SMALL_PANEL_AUDIT"
    assert paired["split_count"] == 12
    assert paired["active_splits"] == 8
    assert paired["inactive_exact_fallback"] is True
    overall = paired["overall"]
    assert overall["condition_seed_units"] == 654
    close(overall["relative_all_gene_mse_gain"], 0.001079217249737252)
    close(overall["relative_top100_mse_gain"], 0.0030453806218846413)
    close(overall["top100_pcc_delta"], 0.0013006984799760163)

    aggregate = args.repo / "experiments/18_dual_head_evidence_gate/results/joint_k10_official/aggregate"
    formal = json.loads((aggregate / "formal_verdict.json").read_text())
    assert formal["status"] == "PASS_FORMAL_AGGREGATION"
    assert formal["official_primary_operation_seed_units"] == 654
    close(formal["official_primary_top100_mean_rank"], 3.871559633027523)
    assert formal["official_primary_top100_position"] == 1
    assert formal["both_panels_robustness_position"] == 1
    assert formal["dataset_equal_weight_top100_position"] == 1

    with (aggregate / "formal_top100_leaderboard.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert rows[0]["method"] == "WitnessCell"
    close(float(rows[0]["mean_operation_rank"]), 3.871559633027523)
    assert float(rows[0]["official_position"]) == 1.0

    bootstrap = json.loads((aggregate / "paired_cluster_bootstrap.json").read_text())
    scouter = bootstrap["comparisons"]["scouter"]
    close(scouter["rank_difference_witness_minus_competitor"], -1.5412844036697249)
    close(scouter["cluster_bootstrap_95_ci"][0], -1.9540636702850922)
    close(scouter["cluster_bootstrap_95_ci"][1], -1.1267079136021945)

    result = {
        "status": "PASS_PREREGISTERED_COMPUTATIONAL_REPLICATION",
        "protocol_version": "2.0.0",
        "prior_success_disclosed": True,
        "C1_execution_schema": "PASS",
        "C2_numerical_identity": "PASS",
        "C3_gate_fallback_identity": "PASS",
        "C4_paired_mean_response": "PASS",
        "C5_inherited_gate18_aggregation_integrity": "PASS",
        "prediction_archives": 12,
        "arrays_compared": array_comparisons,
        "global_max_abs_delta": global_max_abs_delta,
        "correction_active_splits": 8,
        "exact_fallback_splits": 4,
        "paired_units": 654,
        "all_gene_mse_gain_vs_v13": overall["relative_all_gene_mse_gain"],
        "top100_mse_gain_vs_v13": overall["relative_top100_mse_gain"],
        "top100_pcc_delta_vs_v13": overall["top100_pcc_delta"],
        "packaged_official_top100_mean_rank": formal["official_primary_top100_mean_rank"],
        "packaged_official_top100_position": formal["official_primary_top100_position"],
        "packaged_aggregation_origin": "Gate18 inherited evidence; not Gate19/v14 raw-h5ad rescoring",
        "raw_h5ad_six_metric_rescoring_performed": False,
        "new_biological_measurements": False,
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
