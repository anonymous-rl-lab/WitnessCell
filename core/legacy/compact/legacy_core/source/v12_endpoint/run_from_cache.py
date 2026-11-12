#!/usr/bin/env python3
"""Run frozen WitnessCell v1 and endpoint-gated v2 from one moment cache."""
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
import witnesscell_v2_core


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
    parser.add_argument("--identity-mode", choices=("go_shrinkage", "self_only"), default="go_shrinkage")
    parser.add_argument("--go-top-k", type=int, default=40)
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
    upgraded = witnesscell_v2_core.fit_predict(
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
        "status": "PASS_WITNESSCELL_V2_ENDPOINT_EVIDENCE_PREDICTION",
        "method": "WitnessCell_endpoint_v2",
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
            "identity_mode": args.identity_mode,
            "anchor_mode": getattr(head, "anchor_mode", None),
            "self_summary": getattr(head, "self_summary", "mean"),
            "go_top_k": getattr(head, "top_k", None),
            "go_residual_weight": getattr(head, "go_residual_weight", None),
            "weight_mode": args.weight_mode,
            "gate_metric": head.gate_metric,
            "gate_lower_bound": head.gate_lower_bound,
            "background_weight": head.background_weight,
            "self_weight": head.self_weight,
            "loo_all_gene_gain_mean": head.loo_all_gene_gain_mean,
            "loo_all_gene_gain_lower": head.loo_all_gene_gain_lower,
            "loo_top100_gain_mean": head.loo_top100_gain_mean,
            "loo_top100_gain_lower": head.loo_top100_gain_lower,
            "loo_pcc_delta_mean": head.loo_pcc_delta_mean,
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
