#!/usr/bin/env python3
"""Run frozen WitnessCell v1 and incremental-gated WitnessCell v4."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pickle
from pathlib import Path

import numpy as np

from official_split import (
    filter_gears_supported_conditions,
    load_gears_supported_genes,
    make_official_split,
)
import witnesscell_v1_core
import witnesscell_v4_core


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_records(path: Path, records: tuple[dict, ...]) -> None:
    if not records:
        path.write_text("condition\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--dataset")
    parser.add_argument("--seed", type=int, choices=(1, 2, 3), required=True)
    parser.add_argument("--gears-assets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--anchor-mode", choices=("add", "replace"), default="add")
    parser.add_argument("--self-summary", choices=("mean", "median"), default="mean")
    parser.add_argument("--weight-mode", choices=("fixed", "closed_form"), default="closed_form")
    parser.add_argument("--gate-metric", choices=("top100_mse", "all_gene_mse"), default="all_gene_mse")
    parser.add_argument("--identity-mode", choices=("multihead_v14",), default="multihead_v14")
    parser.add_argument("--go-top-k", type=int, default=10)
    parser.add_argument(
        "--dual-dense-mode",
        choices=("auto", "mean_only", "mean_go", "mean_self", "mean_go_self"),
        default="mean_go",
    )
    parser.add_argument(
        "--dual-sparse-mode",
        choices=("self_only", "joint_self", "go_only", "go_self", "scale_dense"),
        default="joint_self",
    )
    parser.add_argument(
        "--dual-support-mode",
        choices=(
            "go_residual", "go_program", "global_activity", "hybrid",
            "go_de_vote", "global_de_frequency", "de_hybrid",
        ),
        default="go_residual",
    )
    parser.add_argument("--dual-support-k", type=int, default=1)
    parser.add_argument(
        "--dual-direction-gate-metric",
        choices=("top100_pcc", "top100_mse"),
        default="top100_pcc",
    )
    nesting = parser.add_mutually_exclusive_group()
    nesting.add_argument("--nested-gate", dest="nested_gate", action="store_true")
    nesting.add_argument("--no-nested-gate", dest="nested_gate", action="store_false")
    parser.set_defaults(nested_gate=False)
    parser.add_argument("--require-incremental-direction", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    data = np.load(args.cache, allow_pickle=False)
    dataset = args.dataset or str(data["dataset"])
    conditions = data["conditions"].astype(str).tolist()
    genes = data["genes"].astype(str).tolist()
    matrix = data["means"].astype(float)
    variance_matrix = data["variances"].astype(float)
    count_array = data["counts"].astype(np.int64)
    index = {condition: position for position, condition in enumerate(conditions)}
    means = {condition: matrix[position] for position, condition in enumerate(conditions)}
    variances = {
        condition: variance_matrix[position]
        for position, condition in enumerate(conditions)
    }
    counts = {
        condition: int(count_array[position])
        for position, condition in enumerate(conditions)
    }

    supported_genes = load_gears_supported_genes(args.gears_assets)
    with (args.gears_assets / "gene2go_all.pkl").open("rb") as stream:
        gene2go = pickle.load(stream)
    official_conditions = filter_gears_supported_conditions(conditions, supported_genes)
    split = make_official_split(official_conditions, args.seed)
    fit_names = {"control", *split.train, *split.validation}
    fit_means = {condition: means[condition] for condition in fit_names}
    fit_variances = {condition: variances[condition] for condition in fit_names}
    fit_counts = {condition: counts[condition] for condition in fit_names}

    baseline = witnesscell_v1_core.fit_predict(
        fit_means,
        list(split.train),
        list(split.validation),
        list(split.test),
    )
    upgraded = witnesscell_v4_core.fit_predict(
        fit_means,
        fit_variances,
        fit_counts,
        genes,
        list(split.train),
        list(split.validation),
        list(split.test),
        anchor_mode=args.anchor_mode,
        self_summary=args.self_summary,
        weight_mode=args.weight_mode,
        gate_metric=args.gate_metric,
        identity_mode=args.identity_mode,
        gene2go=gene2go,
        go_top_k=args.go_top_k,
        dual_dense_mode=args.dual_dense_mode,
        dual_sparse_mode=args.dual_sparse_mode,
        dual_support_mode=args.dual_support_mode,
        dual_support_k=args.dual_support_k,
        dual_direction_gate_metric=args.dual_direction_gate_metric,
        dual_nested_gate=args.nested_gate,
        dual_require_incremental_direction=args.require_incremental_direction,
    )

    prediction = np.stack([upgraded.predictions[c] for c in split.test]).astype(np.float32)
    factorized = np.stack([
        upgraded.factorized_predictions[c] for c in split.test
    ]).astype(np.float32)
    baseline_prediction = np.stack([
        baseline.predictions[c] for c in split.test
    ]).astype(np.float32)
    baseline_factorized = np.stack([
        baseline.factorized_predictions[c] for c in split.test
    ]).astype(np.float32)
    truth = np.stack([means[c] for c in split.test]).astype(np.float32)
    truth_variance = np.stack([variances[c] for c in split.test]).astype(np.float32)
    test_counts = np.asarray([counts[c] for c in split.test], dtype=np.int64)
    subgroup = {
        condition: name
        for name, values in split.test_subgroup.items()
        for condition in values
    }
    subgroup_array = np.asarray([subgroup[c] for c in split.test], dtype=str)
    control = means["control"].astype(np.float32)
    common = {
        "dataset": np.asarray(dataset),
        "seed": np.asarray(args.seed, dtype=np.int64),
        "conditions": np.asarray(split.test, dtype=str),
        "subgroups": subgroup_array,
        "genes": np.asarray(genes, dtype=str),
        "prediction": prediction,
        "factorized_prediction": factorized,
        "baseline_prediction": baseline_prediction,
        "baseline_factorized_prediction": baseline_factorized,
        "control": control,
        "control_variance": variances["control"].astype(np.float32),
        "control_count": np.asarray(counts["control"], dtype=np.int64),
        "test_cell_counts": test_counts,
    }
    np.savez_compressed(
        args.out / "predictions.npz",
        **common,
        truth=truth,
        truth_variance=truth_variance,
    )
    # Target-free deployment archive: no held-out means or variances.
    np.savez_compressed(args.out / "deploy_predictions.npz", **common)
    write_records(args.out / "identity_loo_records.csv", upgraded.identity_head.records)

    head = upgraded.identity_head
    max_baseline_factorized_delta = float(np.max(np.abs(
        baseline_factorized - factorized
    ))) if not head.active else None
    manifest = {
        "status": "PASS_WITNESSCELL_V4_INCREMENTAL_AMPLITUDE_EVIDENCE_PREDICTION",
        "method": "WitnessCell_endpoint_v4_incremental_amplitude",
        "dataset": dataset,
        "seed": args.seed,
        "cache": str(args.cache.resolve()),
        "cache_sha256": sha256(args.cache),
        "official_gears_filter": {
            "raw_conditions": len(conditions),
            "supported_conditions": len(official_conditions),
            "removed_conditions": sorted(set(conditions) - set(official_conditions)),
        },
        "split_counts": {
            "train": len(split.train),
            "validation": len(split.validation),
            "test": len(split.test),
        },
        "test_subgroup_counts": {
            key: len(value) for key, value in split.test_subgroup.items()
        },
        "identity_head": {
            "active": head.active,
            "dense_active": getattr(head, "dense_active", head.active),
            "sparse_active": getattr(head, "sparse_active", False),
            "identity_mode": args.identity_mode,
            "anchor_mode": getattr(head, "anchor_mode", None),
            "self_summary": getattr(head, "self_summary", "mean"),
            "go_top_k": getattr(head, "go_top_k", getattr(head, "top_k", None)),
            "go_residual_weight": getattr(head, "go_residual_weight", None),
            "dense_mode": getattr(head, "dense_mode", None),
            "selected_recipe": getattr(head, "selected_recipe", "fixed"),
            "candidate_training_scores": getattr(head, "candidate_training_scores", ()),
            "sparse_mode": getattr(head, "sparse_mode", None),
            "support_mode": getattr(head, "support_mode", None),
            "support_k": getattr(head, "support_k", None),
            "dense_weights": getattr(head, "dense_weights", None),
            "sparse_weights": getattr(head, "sparse_weights", None),
            "correction_active": getattr(head, "correction_active", False),
            "correction_weights": getattr(head, "correction_weights", None),
            "all_gene_upgrade_mean": getattr(head, "all_gene_upgrade_mean", None),
            "all_gene_upgrade_lower": getattr(head, "all_gene_upgrade_lower", None),
            "top100_mse_upgrade_mean": getattr(head, "top100_mse_upgrade_mean", None),
            "top100_mse_upgrade_lower": getattr(head, "top100_mse_upgrade_lower", None),
            "top100_pcc_upgrade_mean": getattr(head, "top100_pcc_upgrade_mean", None),
            "top100_pcc_upgrade_lower": getattr(head, "top100_pcc_upgrade_lower", None),
            "combined_dense_weights": getattr(head, "combined_dense_weights", None),
            "direction_gate_metric": getattr(head, "direction_gate_metric", None),
            "incremental_direction_required": getattr(head, "incremental_direction_required", False),
            "dense_all_gene_gain_lower": getattr(head, "dense_all_gene_gain_lower", None),
            "final_all_gene_gain_lower": getattr(head, "final_all_gene_gain_lower", None),
            "final_top100_pcc_delta_lower": getattr(head, "final_top100_pcc_delta_lower", None),
            "sparse_incremental_top100_pcc_delta_lower": getattr(head, "sparse_incremental_top100_pcc_delta_lower", None),
            "weight_mode": args.weight_mode,
            "gate_metric": getattr(head, "gate_metric", "dual"),
            "gate_lower_bound": getattr(head, "gate_lower_bound", None),
            "background_weight": getattr(head, "background_weight", None),
            "self_weight": getattr(head, "self_weight", None),
            "loo_all_gene_gain_mean": getattr(head, "loo_all_gene_gain_mean", getattr(head, "final_all_gene_gain_mean", None)),
            "loo_all_gene_gain_lower": getattr(head, "loo_all_gene_gain_lower", getattr(head, "final_all_gene_gain_lower", None)),
            "loo_top100_gain_mean": getattr(head, "loo_top100_gain_mean", getattr(head, "final_top100_gain_mean", None)),
            "loo_top100_gain_lower": getattr(head, "loo_top100_gain_lower", getattr(head, "final_top100_gain_lower", None)),
            "loo_pcc_delta_mean": getattr(head, "loo_pcc_delta_mean", getattr(head, "final_top100_pcc_delta_mean", None)),
            "known_single_count": head.known_single_count,
        },
        "selected_v1": {
            "alpha": baseline.alpha,
            "noise_ratio": baseline.noise_ratio,
            "gamma": baseline.gamma,
        },
        "selected_v2": {
            "alpha": upgraded.alpha,
            "noise_ratio": upgraded.noise_ratio,
            "gamma": upgraded.gamma,
        },
        "inactive_gate_max_factorized_delta_vs_v1": max_baseline_factorized_delta,
        "leakage_contract": "fit receives only control, official train and official validation; endpoint weights and activation are fit only from training singles; deploy archive contains no test truth",
    }
    (args.out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
