#!/usr/bin/env python3
"""Formal candidate-blind Phase M scoring from a pre-frozen metric contract."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from metric_core import MetricInputError, WeightVector, drf, nir
from phase_m_metric_validity import safe_metrics, sha256


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--split-ledger", type=Path, required=True)
    parser.add_argument("--contract-seal", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    ledger = json.loads(args.split_ledger.read_text())
    seal = json.loads(args.contract_seal.read_text())
    dataset = str(ledger["dataset"])
    actual_contract_sha = sha256(args.contract)
    expected = seal["contracts"][dataset]
    if actual_contract_sha != expected["contract_sha256"]:
        raise MetricInputError(f"contract hash drift for {dataset}")
    if sha256(args.split_ledger) != expected["split_ledger_sha256"]:
        raise MetricInputError(f"split-ledger hash drift for {dataset}")

    with np.load(args.contract, allow_pickle=False) as contract:
        conditions = contract["conditions"].astype(str).tolist()
        genes = contract["genes"].astype(str)
        full_counts = contract["full_counts"].astype(np.int64)
        full_means = contract["full_means"].astype(float)
        first_means = contract["first_means"].astype(float)
        second_means = contract["second_means"].astype(float)
        weights = contract["weights"].astype(float)
        weight_evaluable = contract["weight_evaluable"].astype(bool)
        second_pvalues = contract["second_pvalues_adjusted"].astype(float)
        control = contract["control"].astype(float)
    condition_index = {condition: index for index, condition in enumerate(conditions)}

    validity_rows = []
    nir_rows = []
    for seed_text, split in ledger["split_contract"].items():
        seed = int(seed_text)
        train_noncontrol = [value for value in split["train"] if value != "control"]
        train_mean = np.mean([full_means[condition_index[value]] for value in train_noncontrol], axis=0)
        if hashlib.sha256(train_mean.tobytes()).hexdigest() != split["training_condition_mean_sha256"]:
            raise MetricInputError(f"training mean drift for {dataset} seed {seed}")

        ids = []
        truth_rows = []
        control_rows: dict[str, list[np.ndarray]] = {
            "mean_negative": [],
            "control_negative": [],
            "technical_duplicate": [],
            "interpolated_duplicate": [],
        }
        for target in split["test"]:
            index = condition_index[target]
            if not weight_evaluable[index]:
                validity_rows.append(
                    {"dataset": dataset, "seed": seed, "condition": target, "evaluable": False}
                )
                continue
            weight = WeightVector(genes, weights[index])
            truth = first_means[index]
            duplicate = second_means[index]
            alpha = 1.0 - second_pvalues[index]
            interpolated = alpha * duplicate + (1.0 - alpha) * train_mean
            controls = {
                "mean_negative": train_mean,
                "control_negative": control,
                "technical_duplicate": duplicate,
                "interpolated_duplicate": interpolated,
            }
            scored = {
                name: safe_metrics(value, truth, train_mean, control, weight)
                for name, value in controls.items()
            }
            row = {"dataset": dataset, "seed": seed, "condition": target, "evaluable": True}
            for control_name, metrics in scored.items():
                for metric_name, value in metrics.items():
                    row[f"{control_name}__{metric_name}"] = value
            for metric_name, higher in (
                ("wmse", False),
                ("weighted_r2_deltapert", True),
                ("mse", False),
                ("r2_deltapert", True),
                ("pearson_deltactrl", True),
            ):
                negative_name = "control_negative" if metric_name == "pearson_deltactrl" else "mean_negative"
                baseline_value = scored[negative_name][metric_name]
                for positive_name in ("technical_duplicate", "interpolated_duplicate"):
                    try:
                        value = drf(
                            baseline_value,
                            scored[positive_name][metric_name],
                            higher_better=higher,
                        )
                    except MetricInputError:
                        value = float("nan")
                    row[f"drf__{negative_name}__{positive_name}__{metric_name}"] = value
            validity_rows.append(row)
            ids.append(target)
            truth_rows.append(truth)
            for name, value in controls.items():
                control_rows[name].append(value)

        if len(ids) >= 2:
            truth_matrix = np.stack(truth_rows)
            for name, values in control_rows.items():
                scores = nir(np.stack(values), truth_matrix, ids, [dataset] * len(ids))
                for target, value in scores.items():
                    nir_rows.append(
                        {"dataset": dataset, "seed": seed, "condition": target, "control": name, "nir": value}
                    )

    validity = pd.DataFrame(validity_rows)
    nir_frame = pd.DataFrame(nir_rows)
    validity_path = args.out / f"{dataset}.metric_validity.csv"
    nir_path = args.out / f"{dataset}.nir_controls.csv"
    validity.to_csv(validity_path, index=False)
    nir_frame.to_csv(nir_path, index=False)
    report = {
        "status": "PASS_FORMAL_PHASE_M_CANDIDATE_BLIND_CONTROLS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": dataset,
        "contract_sha256": actual_contract_sha,
        "rows": len(validity),
        "evaluable_rows": int(validity["evaluable"].sum()),
        "metric_validity_sha256": sha256(validity_path),
        "nir_controls_sha256": sha256(nir_path),
    }
    (args.out / f"{dataset}.phase_m_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
