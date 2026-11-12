#!/usr/bin/env python3
"""Candidate-blind Phase M: technical splits, weights and metric validity.

This entry point has no argument for, import from, or filesystem dependency on
WitnessCell/model outputs.  It creates the metric contract that Phase P may
consume only after Phase M has been sealed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from metric_core import (
    MetricInputError,
    WeightVector,
    delta_r2,
    drf,
    mse,
    nir,
    pearson_control_delta,
    source_weight_transform,
    weighted_delta_r2,
    wmse,
)


MIN_CELLS_DEG = 4
TECHNICAL_SPLIT_SEED = 42


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity_hash(values: list[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def as_dense(value) -> np.ndarray:
    if hasattr(value, "toarray"):
        value = value.toarray()
    return np.asarray(value)


def technical_split(labels: np.ndarray, cell_ids: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    """Reproduce the locked global-seed, sequential-permutation split."""

    np.random.seed(TECHNICAL_SPLIT_SEED)
    assignment = np.full(labels.size, -1, dtype=np.int8)
    ledger = []
    for condition in pd.unique(labels):
        indices = np.flatnonzero(labels == condition)
        ordered_ids = cell_ids[indices]
        if indices.size < 2:
            ledger.append(
                {
                    "condition": str(condition),
                    "cells": int(indices.size),
                    "first_count": 0,
                    "second_count": 0,
                    "first_identity_sha256": None,
                    "second_identity_sha256": None,
                }
            )
            continue
        permutation = np.random.permutation(indices.size)
        split = indices.size // 2
        first_local = permutation[:split]
        second_local = permutation[split:]
        assignment[indices[first_local]] = 0
        assignment[indices[second_local]] = 1
        first_ids = ordered_ids[first_local].astype(str).tolist()
        second_ids = ordered_ids[second_local].astype(str).tolist()
        ledger.append(
            {
                "condition": str(condition),
                "cells": int(indices.size),
                "first_count": len(first_ids),
                "second_count": len(second_ids),
                "first_identity_sha256": identity_hash(first_ids),
                "second_identity_sha256": identity_hash(second_ids),
            }
        )
    return assignment, ledger


def group_centroids(
    matrix: np.ndarray,
    labels: np.ndarray,
    conditions: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    means = np.zeros((len(conditions), matrix.shape[1]), dtype=np.float64)
    variances = np.zeros((len(conditions), matrix.shape[1]), dtype=np.float64)
    counts = np.zeros(len(conditions), dtype=np.int64)
    for index, condition in enumerate(conditions):
        selected = matrix[labels == condition]
        counts[index] = selected.shape[0]
        if selected.shape[0]:
            means[index] = np.mean(selected, axis=0, dtype=np.float64)
            variances[index] = np.var(selected, axis=0, dtype=np.float64)
    return means, variances, counts


def rank_tables(adata_half: ad.AnnData, condition_key: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = adata_half.obs[condition_key].astype(str).value_counts()
    valid = counts[counts >= MIN_CELLS_DEG].index.tolist()
    valid = [value for value in valid if "control" not in value.lower() and "ctrl" not in value.lower()]
    subset = adata_half[adata_half.obs[condition_key].astype(str).isin(valid)].copy()
    subset.obs[condition_key] = subset.obs[condition_key].astype("category")
    with warnings.catch_warnings():
        # Scanpy also computes unused log-fold-change display columns. Formal
        # weights use only the t scores/p-values; suppressing display warnings
        # does not alter those arrays.
        warnings.simplefilter("ignore", RuntimeWarning)
        warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
        sc.tl.rank_genes_groups(
            subset,
            condition_key,
            method="t-test_overestim_var",
            reference="rest",
        )
    result = subset.uns["rank_genes_groups"]
    return pd.DataFrame(result["names"]), pd.DataFrame(result["scores"]), pd.DataFrame(result["pvals_adj"])


def aligned_weights(
    names: pd.DataFrame, scores: pd.DataFrame, genes: np.ndarray, conditions: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.zeros((len(conditions), genes.size), dtype=np.float64)
    evaluable = np.zeros(len(conditions), dtype=bool)
    for index, condition in enumerate(conditions):
        if condition not in scores.columns:
            continue
        try:
            weight = source_weight_transform(scores[condition].to_numpy(), names[condition].astype(str), genes)
        except MetricInputError:
            continue
        matrix[index] = weight.values
        evaluable[index] = True
    return matrix, evaluable


def aligned_adjusted_pvalues(
    names: pd.DataFrame, pvalues: pd.DataFrame, genes: np.ndarray, conditions: list[str]
) -> np.ndarray:
    matrix = np.ones((len(conditions), genes.size), dtype=np.float64)
    for index, condition in enumerate(conditions):
        if condition not in pvalues.columns:
            continue
        names_array = names[condition].astype(str).to_numpy()
        if len(set(names_array.tolist())) != names_array.size:
            raise MetricInputError(f"duplicate second-half ranked gene for {condition}")
        lookup = pd.Series(pvalues[condition].to_numpy(dtype=float), index=names_array)
        aligned = lookup.reindex(genes, fill_value=1.0).to_numpy(dtype=float)
        matrix[index] = np.nan_to_num(aligned, nan=1.0)
    return matrix


def load_official_split_module(module_dir: Path):
    sys.path.insert(0, str(module_dir))
    from official_split import (  # type: ignore
        filter_gears_supported_conditions,
        load_gears_supported_genes,
        make_official_split,
    )

    return filter_gears_supported_conditions, load_gears_supported_genes, make_official_split


def safe_metrics(pred: np.ndarray, truth: np.ndarray, train_mean: np.ndarray, control: np.ndarray, weight: WeightVector) -> dict:
    output = {}
    functions = {
        "wmse": lambda: wmse(pred, truth, weight),
        "weighted_r2_deltapert": lambda: weighted_delta_r2(pred, truth, train_mean, weight),
        "mse": lambda: mse(pred, truth),
        "r2_deltapert": lambda: delta_r2(pred, truth, train_mean),
        "pearson_deltactrl": lambda: pearson_control_delta(pred, truth, control),
    }
    for name, function in functions.items():
        try:
            output[name] = function()
        except (MetricInputError, ValueError, FloatingPointError):
            output[name] = float("nan")
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--official-module", type=Path, required=True)
    parser.add_argument("--gears-assets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    source_sha = sha256(args.data)
    adata = ad.read_h5ad(args.data)
    condition_key = "perturbation"
    if condition_key not in adata.obs:
        raise MetricInputError(f"{condition_key} missing from {args.data}")
    labels = adata.obs[condition_key].astype(str).to_numpy()
    cell_ids = adata.obs_names.astype(str).to_numpy()
    # AnnData/Pandas can expose string indices as an object-dtype array.  A
    # formal contract must remain loadable with allow_pickle=False, so force a
    # fixed-width Unicode dtype at the serialization boundary.
    genes = np.asarray(adata.var_names.astype(str).tolist(), dtype=str)
    if len(set(genes.tolist())) != genes.size:
        raise MetricInputError("formal h5ad gene identifiers are not unique")
    matrix = as_dense(adata.X)
    if not np.all(np.isfinite(matrix)):
        raise MetricInputError("formal cell matrix contains non-finite values")

    conditions = pd.unique(labels).astype(str).tolist()
    assignment, ledger = technical_split(labels, cell_ids)
    if np.any(assignment < 0):
        raise MetricInputError("one or more cells could not enter a technical duplicate split")

    full_means, full_variances, full_counts = group_centroids(matrix, labels, conditions)
    first_labels = labels[assignment == 0]
    second_labels = labels[assignment == 1]
    first_means, _, first_counts = group_centroids(matrix[assignment == 0], first_labels, conditions)
    second_means, _, second_counts = group_centroids(matrix[assignment == 1], second_labels, conditions)

    first = ad.AnnData(X=matrix[assignment == 0], obs=adata.obs.iloc[np.flatnonzero(assignment == 0)].copy(), var=adata.var.copy())
    second = ad.AnnData(X=matrix[assignment == 1], obs=adata.obs.iloc[np.flatnonzero(assignment == 1)].copy(), var=adata.var.copy())
    first_names, first_scores, _ = rank_tables(first, condition_key)
    second_names, _, second_pvalues_adjusted = rank_tables(second, condition_key)
    weights, weight_evaluable = aligned_weights(first_names, first_scores, genes, conditions)
    second_pvalues = aligned_adjusted_pvalues(second_names, second_pvalues_adjusted, genes, conditions)

    filter_supported, load_supported, make_split = load_official_split_module(args.official_module)
    supported_genes = load_supported(args.gears_assets)
    supported_conditions = filter_supported(conditions, supported_genes)
    condition_index = {condition: index for index, condition in enumerate(conditions)}
    if "control" not in condition_index:
        raise MetricInputError("formal dataset has no condition named control")
    control = full_means[condition_index["control"]]

    split_contract = {}
    validity_rows = []
    nir_rows = []
    for seed in (1, 2, 3):
        split = make_split(supported_conditions, seed)
        train_noncontrol = [value for value in split.train if value != "control"]
        train_mean = np.mean([full_means[condition_index[value]] for value in train_noncontrol], axis=0)
        split_contract[str(seed)] = {
            "train": list(split.train),
            "validation": list(split.validation),
            "test": list(split.test),
            "test_subgroup": {key: list(value) for key, value in split.test_subgroup.items()},
            "training_condition_mean_sha256": hashlib.sha256(train_mean.tobytes()).hexdigest(),
        }

        if args.prepare_only:
            continue

        eligible_for_nir = []
        nir_truth = []
        nir_mean = []
        nir_duplicate = []
        nir_interpolated = []
        nir_control = []
        for target in split.test:
            index = condition_index[target]
            if not weight_evaluable[index]:
                validity_rows.append(
                    {"dataset": args.dataset, "seed": seed, "condition": target, "evaluable": False}
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
            scored = {name: safe_metrics(value, truth, train_mean, control, weight) for name, value in controls.items()}
            row = {"dataset": args.dataset, "seed": seed, "condition": target, "evaluable": True}
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
                    positive_value = scored[positive_name][metric_name]
                    try:
                        row[f"drf__{negative_name}__{positive_name}__{metric_name}"] = drf(
                            baseline_value, positive_value, higher_better=higher
                        )
                    except MetricInputError:
                        row[f"drf__{negative_name}__{positive_name}__{metric_name}"] = float("nan")
            validity_rows.append(row)

            eligible_for_nir.append(target)
            nir_truth.append(truth)
            nir_mean.append(train_mean)
            nir_duplicate.append(duplicate)
            nir_interpolated.append(interpolated)
            nir_control.append(control)

        if len(eligible_for_nir) >= 2:
            truth_matrix = np.stack(nir_truth)
            for control_name, values in (
                ("mean_negative", nir_mean),
                ("control_negative", nir_control),
                ("technical_duplicate", nir_duplicate),
                ("interpolated_duplicate", nir_interpolated),
            ):
                scores = nir(np.stack(values), truth_matrix, eligible_for_nir, [args.dataset] * len(eligible_for_nir))
                for target, value in scores.items():
                    nir_rows.append(
                        {"dataset": args.dataset, "seed": seed, "condition": target, "control": control_name, "nir": value}
                    )

    np.savez_compressed(
        args.out / f"{args.dataset}.phase_m_contract.npz",
        dataset=np.asarray(args.dataset),
        source_sha256=np.asarray(source_sha),
        conditions=np.asarray(conditions, dtype=str),
        genes=genes,
        full_counts=full_counts,
        full_variances=full_variances,
        first_counts=first_counts,
        second_counts=second_counts,
        full_means=full_means,
        first_means=first_means,
        second_means=second_means,
        weights=weights,
        weight_evaluable=weight_evaluable,
        second_pvalues_adjusted=second_pvalues,
        control=control,
    )
    if not args.prepare_only:
        pd.DataFrame(validity_rows).to_csv(args.out / f"{args.dataset}.metric_validity.csv", index=False)
        pd.DataFrame(nir_rows).to_csv(args.out / f"{args.dataset}.nir_controls.csv", index=False)
    ledger_payload = {
        "dataset": args.dataset,
        "source_h5ad": str(args.data.resolve()),
        "source_sha256": source_sha,
        "technical_split_seed": TECHNICAL_SPLIT_SEED,
        "algorithm": "np.random.seed(42) then sequential np.random.permutation per pd.unique condition",
        "condition_order": conditions,
        "cell_split_ledger": ledger,
        "split_contract": split_contract,
    }
    (args.out / f"{args.dataset}.split_ledger.json").write_text(json.dumps(ledger_payload, indent=2) + "\n")
    report = {
        "status": (
            "PASS_PREFREEZE_CANDIDATE_BLIND_CONTRACT"
            if args.prepare_only
            else "PASS_PHASE_M_CANDIDATE_BLIND_METRIC_CONTRACT"
        ),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": args.dataset,
        "source_sha256": source_sha,
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "conditions": len(conditions),
        "supported_conditions": len(supported_conditions),
        "evaluable_weight_conditions": int(weight_evaluable.sum()),
        "prepare_only": bool(args.prepare_only),
        "formal_rows": None if args.prepare_only else len(validity_rows),
        "formal_evaluable_rows": None
        if args.prepare_only
        else int(sum(bool(row["evaluable"]) for row in validity_rows)),
        "contract_sha256": sha256(args.out / f"{args.dataset}.phase_m_contract.npz"),
        "metric_validity_sha256": None
        if args.prepare_only
        else sha256(args.out / f"{args.dataset}.metric_validity.csv"),
        "split_ledger_sha256": sha256(args.out / f"{args.dataset}.split_ledger.json"),
    }
    (args.out / f"{args.dataset}.phase_m_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
