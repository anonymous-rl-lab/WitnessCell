#!/usr/bin/env python3
"""Compare package output with 12 frozen Gate 19 deployment archives.

This development-only audit reads trusted research pickles solely to reproduce
the historical GEARS filter. It never writes to the research repository and
does not include those assets in any distribution artifact.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT / "src"))

from witnesscell import ConditionMoments, SplitSpec, WitnessCell  # noqa: E402

DATASETS = ("Norman", "Replogle_exp6", "Schmidt", "Wessels")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--research-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--atol", type=float, default=2e-5)
    arguments = parser.parse_args()
    root = arguments.research_root.resolve()
    source = root / "experiments/19_v14_incremental_amplitude_gate/src"
    sys.path.insert(0, str(source))
    import witnesscell_v4_core  # type: ignore[import-not-found]  # noqa: E402, PLC0415
    from official_split import (  # type: ignore[import-not-found]  # noqa: E402, PLC0415
        filter_gears_supported_conditions,
        load_gears_supported_genes,
        make_official_split,
    )

    assets = root / "experiments/15_scperturbench_sota/module/data/gears_assets"
    supported = load_gears_supported_genes(assets)
    with (assets / "gene2go_all.pkl").open("rb") as stream:
        gene2go = pickle.load(stream)  # noqa: S301 - trusted frozen research asset

    rows = []
    passed = True
    for dataset in DATASETS:
        cache_path = root / f"experiments/18_dual_head_evidence_gate/cache/{dataset}.condition_moments.npz"
        with np.load(cache_path, allow_pickle=False) as cache:
            conditions = cache["conditions"].astype(str)
            genes = cache["genes"].astype(str)
            means_matrix = cache["means"].astype(float)
            variance_matrix = cache["variances"].astype(float)
            count_array = cache["counts"].astype(np.int64)
        means = dict(zip(conditions, means_matrix, strict=True))
        variances = dict(zip(conditions, variance_matrix, strict=True))
        counts = dict(zip(conditions, count_array, strict=True))
        moments = ConditionMoments.from_mappings(
            genes=genes, means=means, variances=variances, counts=counts
        )
        official = filter_gears_supported_conditions(conditions.tolist(), supported)
        for seed in (1, 2, 3):
            split = make_official_split(official, seed)
            model = WitnessCell().fit(
                moments,
                SplitSpec.create(split.train, split.validation),
                gene2go=gene2go,
            )
            frozen_path = root / f"experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo/{dataset}/seed{seed}/deploy_predictions.npz"
            manifest_path = frozen_path.with_name("manifest.json")
            with np.load(frozen_path, allow_pickle=False) as frozen:
                targets = tuple(frozen["conditions"].astype(str))
                expected = frozen["prediction"].astype(float)
                expected_factorized = frozen["factorized_prediction"].astype(float)
            reference_fit = witnesscell_v4_core.fit_predict(
                means,
                variances,
                counts,
                genes,
                list(split.train),
                list(split.validation),
                list(targets),
                identity_mode="multihead_v14",
                gene2go=gene2go,
            )
            prediction = model.predict(targets)
            witness_delta = float(np.max(np.abs(prediction.means - expected)))
            factorized_delta = float(
                np.max(np.abs(prediction.factorized_means - expected_factorized))
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            expected_selected = manifest["selected_v2"]
            expected_head = manifest["identity_head"]
            reference_head = reference_fit.identity_head
            diagnostics = model.diagnostics()
            settings_match = bool(
                np.isclose(diagnostics["alpha"], expected_selected["alpha"])
                and np.isclose(diagnostics["noise_ratio"], expected_selected["noise_ratio"])
                and np.isclose(diagnostics["gamma"], expected_selected["gamma"])
            )
            gates_match = bool(
                diagnostics["dense_active"] == reference_head.baseline_head.dense_active
                and diagnostics["sparse_active"] == reference_head.baseline_head.sparse_active
                and diagnostics["amplitude_active"] == reference_head.correction_active
                and diagnostics["dense_active"] == expected_head["dense_active"]
                and diagnostics["amplitude_active"] == expected_head["correction_active"]
            )
            local_pass = witness_delta <= arguments.atol and factorized_delta <= arguments.atol and settings_match and gates_match
            passed = passed and local_pass
            rows.append(
                {
                    "dataset": dataset,
                    "seed": seed,
                    "max_abs_witness_delta": witness_delta,
                    "max_abs_factorized_delta": factorized_delta,
                    "settings_match": settings_match,
                    "gates_match": gates_match,
                    "passed": local_pass,
                }
            )
    report = {
        "contract": "frozen-gate19-v14-deployment-parity",
        "splits": len(rows),
        "tolerance": arguments.atol,
        "passed": passed,
        "max_abs_witness_delta": max(row["max_abs_witness_delta"] for row in rows),
        "max_abs_factorized_delta": max(row["max_abs_factorized_delta"] for row in rows),
        "rows": rows,
    }
    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(output, encoding="utf-8")
    print(output, end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
