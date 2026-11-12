#!/usr/bin/env python3
"""Formal Phase P: score immutable predictions under the sealed metric contract."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from comparator_core import regenerate_mean_baselines
from metric_core import (
    MetricInputError,
    WeightVector,
    delta_r2,
    mse,
    nir,
    pearson_control_delta,
    weighted_delta_r2,
    wmse,
)


DATASETS = ("Norman", "Wessels", "Schmidt", "Replogle_exp6")
METHODS = ("WitnessCell_v14", "WitnessCell_v13", "trainMean", "baseControl")
BOOTSTRAPS = 20_000
BOOTSTRAP_SEED = 20260822


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_score(
    prediction: np.ndarray,
    truth: np.ndarray,
    training_mean: np.ndarray,
    control: np.ndarray,
    weight: WeightVector,
) -> dict[str, float]:
    functions = {
        "wmse": lambda: wmse(prediction, truth, weight),
        "weighted_r2_deltapert": lambda: weighted_delta_r2(
            prediction, truth, training_mean, weight
        ),
        "mse": lambda: mse(prediction, truth),
        "r2_deltapert": lambda: delta_r2(prediction, truth, training_mean),
        "pearson_deltactrl": lambda: pearson_control_delta(prediction, truth, control),
    }
    output = {}
    for name, function in functions.items():
        try:
            output[name] = float(function())
        except (MetricInputError, ValueError, FloatingPointError):
            output[name] = float("nan")
    return output


def comparison_point(frame: pd.DataFrame, baseline: str) -> dict:
    selected = frame[frame["method"].isin(["WitnessCell_v14", baseline])]
    wide = selected.pivot(
        index=["dataset", "seed", "condition"], columns="method", values=["wmse", "weighted_r2_deltapert"]
    ).dropna()
    split_rows = []
    for (dataset, seed), local in wide.groupby(level=[0, 1]):
        wc_wmse = float(local[("wmse", "WitnessCell_v14")].mean())
        base_wmse = float(local[("wmse", baseline)].mean())
        wc_r2 = float(local[("weighted_r2_deltapert", "WitnessCell_v14")].mean())
        base_r2 = float(local[("weighted_r2_deltapert", baseline)].mean())
        split_rows.append(
            {
                "dataset": dataset,
                "seed": int(seed),
                "baseline": baseline,
                "operations": len(local),
                "wmse_ratio": wc_wmse / base_wmse,
                "weighted_r2_difference": wc_r2 - base_r2,
            }
        )
    splits = pd.DataFrame(split_rows)
    pooled_wc = selected[selected.method == "WitnessCell_v14"].set_index(
        ["dataset", "seed", "condition"]
    )
    pooled_base = selected[selected.method == baseline].set_index(
        ["dataset", "seed", "condition"]
    )
    common = pooled_wc.index.intersection(pooled_base.index)
    return {
        "splits": splits,
        "dataset_equal_wmse_ratio": float(np.exp(np.mean(np.log(splits["wmse_ratio"])))),
        "dataset_equal_weighted_r2_difference": float(splits["weighted_r2_difference"].mean()),
        "operation_pooled_wmse_ratio": float(
            pooled_wc.loc[common, "wmse"].mean() / pooled_base.loc[common, "wmse"].mean()
        ),
        "operation_pooled_weighted_r2_difference": float(
            pooled_wc.loc[common, "weighted_r2_deltapert"].mean()
            - pooled_base.loc[common, "weighted_r2_deltapert"].mean()
        ),
        "operations": len(common),
    }


def cluster_bootstrap(frame: pd.DataFrame, baseline: str) -> tuple[np.ndarray, np.ndarray]:
    """Dataset-stratified condition-cluster bootstrap preserving seed appearances."""

    selected = frame[frame["method"].isin(["WitnessCell_v14", baseline])]
    wide = selected.pivot(
        index=["dataset", "seed", "condition"], columns="method", values=["wmse", "weighted_r2_deltapert"]
    ).dropna()
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    log_ratios = []
    r2_differences = []
    for dataset in DATASETS:
        local = wide.xs(dataset, level="dataset")
        conditions = sorted(set(local.index.get_level_values("condition")))
        n = len(conditions)
        wc_wmse = np.zeros((n, 3), dtype=np.float64)
        base_wmse = np.zeros((n, 3), dtype=np.float64)
        wc_r2 = np.zeros((n, 3), dtype=np.float64)
        base_r2 = np.zeros((n, 3), dtype=np.float64)
        counts = np.zeros((n, 3), dtype=np.float64)
        condition_index = {condition: index for index, condition in enumerate(conditions)}
        for (seed, condition), row in local.iterrows():
            i = condition_index[str(condition)]
            j = int(seed) - 1
            wc_wmse[i, j] = row[("wmse", "WitnessCell_v14")]
            base_wmse[i, j] = row[("wmse", baseline)]
            wc_r2[i, j] = row[("weighted_r2_deltapert", "WitnessCell_v14")]
            base_r2[i, j] = row[("weighted_r2_deltapert", baseline)]
            counts[i, j] = 1.0
        sampled = rng.integers(0, n, size=(BOOTSTRAPS, n))
        denominator = counts[sampled].sum(axis=1)
        if np.any(denominator <= 0):
            raise RuntimeError(f"bootstrap produced an empty split for {dataset}")
        wc_wmse_mean = wc_wmse[sampled].sum(axis=1) / denominator
        base_wmse_mean = base_wmse[sampled].sum(axis=1) / denominator
        wc_r2_mean = wc_r2[sampled].sum(axis=1) / denominator
        base_r2_mean = base_r2[sampled].sum(axis=1) / denominator
        log_ratios.append(np.log(wc_wmse_mean / base_wmse_mean))
        r2_differences.append(wc_r2_mean - base_r2_mean)
    ratio = np.exp(np.mean(np.concatenate(log_ratios, axis=1), axis=1))
    difference = np.mean(np.concatenate(r2_differences, axis=1), axis=1)
    return ratio, difference


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--contracts", type=Path, required=True)
    parser.add_argument("--contract-seal", type=Path, required=True)
    parser.add_argument("--phase-m-manifest", type=Path, required=True)
    parser.add_argument("--scope", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    args.out.mkdir(parents=True, exist_ok=False)
    seal = json.loads(args.contract_seal.read_text())
    if not seal.get("pass"):
        raise RuntimeError("Phase M contract seal did not pass")
    phase_m_verdict_path = args.phase_m_manifest.parent / "metric_validity_verdict.json"
    phase_m_verdict = json.loads(phase_m_verdict_path.read_text())
    if phase_m_verdict["status"] != "PASS_PHASE_M_EXECUTION_AND_SEAL":
        raise RuntimeError("formal Phase M did not seal successfully")
    if phase_m_verdict["phase_m_manifest_sha256"] != sha256(args.phase_m_manifest):
        raise RuntimeError("formal Phase M manifest hash drift")
    for line in args.phase_m_manifest.read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        artifact = args.phase_m_manifest.parent / name
        if sha256(artifact) != digest:
            raise RuntimeError(f"sealed Phase M output drift: {name}")
    scope = json.loads(args.scope.read_text())
    if scope["scope"] != "MANDATORY_BASELINES_ONLY" or not scope["fixed_before_prediction_scoring"]:
        raise RuntimeError("comparator scope is not frozen")

    operation_rows: list[dict] = []
    nir_rows: list[dict] = []
    for dataset in DATASETS:
        contract_path = args.contracts / f"{dataset}.phase_m_contract.npz"
        ledger_path = args.contracts / f"{dataset}.split_ledger.json"
        expected = seal["contracts"][dataset]
        if sha256(contract_path) != expected["contract_sha256"]:
            raise RuntimeError(f"contract drift: {dataset}")
        if sha256(ledger_path) != expected["split_ledger_sha256"]:
            raise RuntimeError(f"split ledger drift: {dataset}")
        ledger = json.loads(ledger_path.read_text())
        with np.load(contract_path, allow_pickle=False) as contract:
            conditions = contract["conditions"].astype(str)
            genes = contract["genes"].astype(str)
            full_counts = contract["full_counts"].astype(np.int64)
            full_means = contract["full_means"].astype(np.float64)
            first_means = contract["first_means"].astype(np.float64)
            weights = contract["weights"].astype(np.float64)
            weight_evaluable = contract["weight_evaluable"].astype(bool)
            control = contract["control"].astype(np.float64)
        condition_index = {condition: index for index, condition in enumerate(conditions)}

        for seed in (1, 2, 3):
            split = ledger["split_contract"][str(seed)]
            targets = list(split["test"])
            train_noncontrol = [value for value in split["train"] if value != "control"]
            metric_mean = np.mean([full_means[condition_index[value]] for value in train_noncontrol], axis=0)
            if hashlib.sha256(metric_mean.tobytes()).hexdigest() != split["training_condition_mean_sha256"]:
                raise RuntimeError(f"metric training mean drift: {dataset} seed {seed}")

            v14_path = (
                repo
                / "experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo"
                / dataset
                / f"seed{seed}/deploy_predictions.npz"
            )
            v13_path = (
                repo
                / "experiments/18_dual_head_evidence_gate/predictions/formal_combo"
                / dataset
                / f"seed{seed}/deploy_predictions.npz"
            )
            with np.load(v14_path, allow_pickle=False) as v14, np.load(v13_path, allow_pickle=False) as v13:
                if v14["conditions"].astype(str).tolist() != targets:
                    raise RuntimeError(f"v14 condition order drift: {dataset} seed {seed}")
                if not np.array_equal(v14["genes"].astype(str), genes):
                    raise RuntimeError(f"v14 gene order drift: {dataset} seed {seed}")
                if not np.array_equal(v13["conditions"], v14["conditions"]):
                    raise RuntimeError(f"v13/v14 condition drift: {dataset} seed {seed}")
                v14_prediction = v14["prediction"].astype(np.float64)
                v13_prediction = v13["prediction"].astype(np.float64)
                subgroups = v14["subgroups"].astype(str)

            train_mean_prediction, base_control_prediction = regenerate_mean_baselines(
                targets,
                split["train"],
                conditions,
                full_means,
                full_counts,
                control,
            )
            predictions = {
                "WitnessCell_v14": v14_prediction,
                "WitnessCell_v13": v13_prediction,
                "trainMean": train_mean_prediction,
                "baseControl": base_control_prediction,
            }
            evaluable_ids = []
            truth_matrix = []
            nir_predictions = {method: [] for method in METHODS}
            for target_position, target in enumerate(targets):
                index = condition_index[target]
                if not weight_evaluable[index]:
                    continue
                try:
                    weight = WeightVector(genes, weights[index])
                except MetricInputError:
                    continue
                truth = first_means[index]
                evaluable_ids.append(target)
                truth_matrix.append(truth)
                for method in METHODS:
                    prediction = predictions[method][target_position]
                    nir_predictions[method].append(prediction)
                    operation_rows.append(
                        {
                            "dataset": dataset,
                            "seed": seed,
                            "condition": target,
                            "subgroup": subgroups[target_position],
                            "method": method,
                            **safe_score(prediction, truth, metric_mean, control, weight),
                        }
                    )
            if len(evaluable_ids) >= 2:
                truth_stack = np.stack(truth_matrix)
                for method in METHODS:
                    scores = nir(
                        np.stack(nir_predictions[method]),
                        truth_stack,
                        evaluable_ids,
                        [dataset] * len(evaluable_ids),
                    )
                    for condition, value in scores.items():
                        nir_rows.append(
                            {
                                "dataset": dataset,
                                "seed": seed,
                                "condition": condition,
                                "method": method,
                                "nir": value,
                            }
                        )

    operations = pd.DataFrame(operation_rows)
    nir_frame = pd.DataFrame(nir_rows)
    operation_path = args.out / "per_operation_metrics.csv"
    nir_path = args.out / "per_operation_nir.csv"
    operations.to_csv(operation_path, index=False)
    nir_frame.to_csv(nir_path, index=False)
    split_summary = operations.groupby(["dataset", "seed", "method"], as_index=False).agg(
        operations=("condition", "size"),
        wmse=("wmse", "mean"),
        weighted_r2_deltapert=("weighted_r2_deltapert", "mean"),
        mse=("mse", "mean"),
        r2_deltapert=("r2_deltapert", "mean"),
        pearson_deltactrl=("pearson_deltactrl", "mean"),
    )
    split_summary = split_summary.merge(
        nir_frame.groupby(["dataset", "seed", "method"], as_index=False)["nir"].mean(),
        on=["dataset", "seed", "method"],
        how="left",
    )
    split_summary.to_csv(args.out / "split_summary.csv", index=False)

    comparison_splits = []
    comparisons = {}
    bootstrap_payload = {}
    for baseline in ("WitnessCell_v13", "trainMean", "baseControl"):
        point = comparison_point(operations, baseline)
        comparison_splits.append(point.pop("splits"))
        ratio_bootstrap, r2_bootstrap = cluster_bootstrap(operations, baseline)
        ratio_ci = np.quantile(ratio_bootstrap, [0.025, 0.975]).tolist()
        r2_ci = np.quantile(r2_bootstrap, [0.025, 0.975]).tolist()
        point.update(
            {
                "wmse_ratio_ci95": [float(value) for value in ratio_ci],
                "weighted_r2_difference_ci95": [float(value) for value in r2_ci],
                "co_primary_pass": bool(
                    point["dataset_equal_wmse_ratio"] < 1.0
                    and ratio_ci[1] < 1.0
                    and point["dataset_equal_weighted_r2_difference"] > 0.0
                    and r2_ci[0] > 0.0
                ),
            }
        )
        comparisons[baseline] = point
        bootstrap_payload[f"{baseline}__wmse_ratio"] = ratio_bootstrap
        bootstrap_payload[f"{baseline}__weighted_r2_difference"] = r2_bootstrap
    pd.concat(comparison_splits, ignore_index=True).to_csv(
        args.out / "split_comparisons.csv", index=False
    )
    np.savez_compressed(args.out / "cluster_bootstrap_arrays.npz", **bootstrap_payload)

    counts = operations.groupby("method").size().to_dict()
    expected_units = 654
    completeness = {
        method: {
            "evaluable": int(counts.get(method, 0)),
            "fraction_of_654": float(counts.get(method, 0) / expected_units),
        }
        for method in METHODS
    }
    twelve_splits_present = all(
        len(operations[(operations.dataset == dataset) & (operations.seed == seed)]) > 0
        for dataset in DATASETS
        for seed in (1, 2, 3)
    )
    completeness_pass = twelve_splits_present and all(
        value["fraction_of_654"] >= 0.90 for value in completeness.values()
    )
    pred_uninformative = bool(
        completeness_pass
        and comparisons["trainMean"]["co_primary_pass"]
        and comparisons["baseControl"]["co_primary_pass"]
    )
    v13 = comparisons["WitnessCell_v13"]
    exact_fallback_audit = json.loads(
        (repo / "experiments/22_metric_calibration_stress/audit/frozen_asset_audit.json").read_text()
    )["checks"]["inactive_archives_exact"]
    if not exact_fallback_audit or v13["dataset_equal_wmse_ratio"] >= 1 or v13["dataset_equal_weighted_r2_difference"] <= 0:
        endpoint_compatibility = "FAIL"
    elif v13["wmse_ratio_ci95"][1] < 1 and v13["weighted_r2_difference_ci95"][0] > 0:
        endpoint_compatibility = "STRONG"
    else:
        endpoint_compatibility = "DIRECTIONAL"

    verdict = {
        "status": "PASS_PHASE_P_EXECUTION",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "formal_structure": "12 dataset-by-seed splits; 20,000 dataset-stratified condition-cluster bootstraps",
        "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_replicates": BOOTSTRAPS,
        "comparator_scope": scope["scope"],
        "completeness": completeness,
        "all_12_splits_present": twelve_splits_present,
        "completeness_pass": completeness_pass,
        "comparisons": comparisons,
        "PRED_LINEAR": "NOT_ADJUDICATED",
        "PRED_UNINFORMATIVE": "PASS" if pred_uninformative else "FAIL",
        "ENDPOINT_COMPATIBILITY": endpoint_compatibility,
        "artifacts": {
            "per_operation_metrics_sha256": sha256(operation_path),
            "per_operation_nir_sha256": sha256(nir_path),
            "cluster_bootstrap_arrays_sha256": sha256(args.out / "cluster_bootstrap_arrays.npz"),
        },
    }
    (args.out / "prediction_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
