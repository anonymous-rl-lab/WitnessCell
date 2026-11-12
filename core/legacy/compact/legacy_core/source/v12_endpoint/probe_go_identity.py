#!/usr/bin/env python3
"""Probe the GO-shrunk background plus self-response endpoint head.

Every coefficient and activation decision is fitted from known-single LOO
rows. Test unseen singles are scored only after the fitted recipe is frozen.
"""
from __future__ import annotations

import argparse
import csv
import json
import pickle
from pathlib import Path

import numpy as np

from identity_head import EPS, one_sided_lower, pearson, self_statistic, welch_scores
from official_split import (
    filter_gears_supported_conditions,
    load_gears_supported_genes,
    make_official_split,
)


TOP_KS = (3, 5, 10, 20, 40, -1)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def jaccard(left: set, right: set) -> float:
    return len(left & right) / max(len(left | right), 1)


def features(
    node: str,
    effects: dict[str, np.ndarray],
    gene_index: dict[str, int],
    gene2go: dict,
    top_k: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    background = np.mean(np.stack(list(effects.values())), axis=0)
    target_go = set(gene2go.get(node, ()))
    neighbors = []
    for known, effect in effects.items():
        similarity = jaccard(target_go, set(gene2go.get(known, ())))
        neighbors.append((similarity, known, effect))
    neighbors.sort(key=lambda row: (row[0], row[1]), reverse=True)
    selected = neighbors if top_k < 0 else neighbors[:top_k]
    weights = np.asarray([row[0] for row in selected], dtype=float)
    if weights.sum() <= EPS:
        go_program = background
    else:
        go_program = np.average(
            np.stack([row[2] for row in selected]), axis=0, weights=weights
        )
    go_residual = go_program - background
    anchor = np.zeros_like(background)
    if node in gene_index:
        anchor[gene_index[node]] = self_statistic(effects, gene_index, "mean")
    return background, go_residual, anchor


def fit_weights(
    feature_rows: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    targets: list[np.ndarray],
    mode: str,
) -> np.ndarray:
    if mode == "full_ols":
        gram = np.zeros((3, 3), dtype=np.float64)
        cross = np.zeros(3, dtype=np.float64)
        for row, target in zip(feature_rows, targets, strict=True):
            design = np.stack(row, axis=1).astype(np.float64)
            gram += design.T @ design
            cross += design.T @ np.asarray(target, dtype=np.float64)
        weight = np.linalg.solve(gram + 1e-8 * np.eye(3), cross)
        return np.clip(weight, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0])
    if mode == "residual_ols":
        gram = np.zeros((2, 2), dtype=np.float64)
        cross = np.zeros(2, dtype=np.float64)
        for (background, go_residual, anchor), target in zip(
            feature_rows, targets, strict=True
        ):
            design = np.stack((go_residual, anchor), axis=1).astype(np.float64)
            residual = np.asarray(target - background, dtype=np.float64)
            gram += design.T @ design
            cross += design.T @ residual
        weight = np.linalg.solve(gram + 1e-8 * np.eye(2), cross)
        weight = np.clip(weight, [0.0, 0.0], [2.0, 2.0])
        return np.asarray([1.0, weight[0], weight[1]])
    raise ValueError(mode)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, action="append", required=True)
    parser.add_argument("--gears-assets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=(1, 2, 3))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.gears_assets / "gene2go_all.pkl").open("rb") as stream:
        gene2go = pickle.load(stream)
    supported = load_gears_supported_genes(args.gears_assets)
    fit_rows: list[dict] = []
    condition_rows: list[dict] = []

    for cache_path in args.cache:
        data = np.load(cache_path, allow_pickle=False)
        dataset = str(data["dataset"])
        conditions = data["conditions"].astype(str).tolist()
        genes = data["genes"].astype(str).tolist()
        matrix = data["means"].astype(float)
        variance_matrix = data["variances"].astype(float)
        count_array = data["counts"].astype(np.int64)
        means = dict(zip(conditions, matrix))
        variances = dict(zip(conditions, variance_matrix))
        counts = dict(zip(conditions, map(int, count_array)))
        official = filter_gears_supported_conditions(conditions, supported)
        gene_index = {gene: index for index, gene in enumerate(genes)}
        control = means["control"]

        for seed in args.seeds:
            split = make_official_split(official, seed)
            effects = {
                condition: means[condition] - control
                for condition in split.train
                if condition != "control" and "+" not in condition
            }
            background = np.mean(np.stack(list(effects.values())), axis=0)
            for top_k in TOP_KS:
                loo_features = []
                loo_targets = []
                held_nodes = sorted(effects)
                for held in held_nodes:
                    remaining = {node: effect for node, effect in effects.items() if node != held}
                    loo_features.append(features(held, remaining, gene_index, gene2go, top_k))
                    loo_targets.append(effects[held])
                for weight_mode in ("full_ols", "residual_ols"):
                    weight = fit_weights(loo_features, loo_targets, weight_mode)
                    loo_all_gain = []
                    loo_top_gain = []
                    loo_pcc_delta = []
                    for held, row, truth in zip(
                        held_nodes, loo_features, loo_targets, strict=True
                    ):
                        baseline = row[0]
                        prediction = sum(value * feature for value, feature in zip(weight, row, strict=True))
                        score = welch_scores(
                            means[held], variances[held], counts[held],
                            control, variances["control"], counts["control"],
                        )
                        order = np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]
                        base_all = float(np.mean(np.square(baseline - truth)))
                        pred_all = float(np.mean(np.square(prediction - truth)))
                        base_top = float(np.mean(np.square(baseline[order] - truth[order])))
                        pred_top = float(np.mean(np.square(prediction[order] - truth[order])))
                        loo_all_gain.append((base_all - pred_all) / max(base_all, EPS))
                        loo_top_gain.append((base_top - pred_top) / max(base_top, EPS))
                        loo_pcc_delta.append(
                            pearson(prediction[order], truth[order])
                            - pearson(baseline[order], truth[order])
                        )
                    all_lower = one_sided_lower(np.asarray(loo_all_gain))
                    top_lower = one_sided_lower(np.asarray(loo_top_gain))
                    # The recorded safety gate used the one-sided lower bound
                    # of the all-gene LOO gain. The top-100 panel remains an
                    # independent directional audit rather than a selector.
                    active = all_lower > 0.0
                    fit_rows.append({
                        "dataset": dataset,
                        "seed": seed,
                        "top_k": top_k,
                        "weight_mode": weight_mode,
                        "active": int(active),
                        "background_weight": weight[0],
                        "go_residual_weight": weight[1],
                        "self_weight": weight[2],
                        "loo_all_gene_gain_mean": np.mean(loo_all_gain),
                        "loo_all_gene_gain_lower": all_lower,
                        "loo_top100_gain_mean": np.mean(loo_top_gain),
                        "loo_top100_gain_lower": top_lower,
                        "loo_pcc_delta_mean": np.mean(loo_pcc_delta),
                    })
                    for condition in split.test_subgroup["unseen_single"]:
                        row = features(condition, effects, gene_index, gene2go, top_k)
                        raw_prediction = sum(value * feature for value, feature in zip(weight, row, strict=True))
                        prediction = raw_prediction if active else background
                        truth = means[condition] - control
                        score = welch_scores(
                            means[condition], variances[condition], counts[condition],
                            control, variances["control"], counts["control"],
                        )
                        order = np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]
                        base_all = float(np.mean(np.square(background - truth)))
                        pred_all = float(np.mean(np.square(prediction - truth)))
                        base_top = float(np.mean(np.square(background[order] - truth[order])))
                        pred_top = float(np.mean(np.square(prediction[order] - truth[order])))
                        condition_rows.append({
                            "dataset": dataset,
                            "seed": seed,
                            "condition": condition,
                            "top_k": top_k,
                            "weight_mode": weight_mode,
                            "active": int(active),
                            "baseline_all_gene_mse": base_all,
                            "candidate_all_gene_mse": pred_all,
                            "relative_all_gene_mse_gain": (base_all - pred_all) / max(base_all, EPS),
                            "baseline_top100_mse": base_top,
                            "candidate_top100_mse": pred_top,
                            "relative_top100_mse_gain": (base_top - pred_top) / max(base_top, EPS),
                            "top100_pcc_delta": pearson(prediction[order], truth[order]) - pearson(background[order], truth[order]),
                        })

    write_csv(args.out / "go_identity_fit_by_split.csv", fit_rows)
    write_csv(args.out / "go_unseen_single_by_condition.csv", condition_rows)
    summary_rows = []
    for top_k in TOP_KS:
        for weight_mode in ("full_ols", "residual_ols"):
            block = [
                row for row in condition_rows
                if row["top_k"] == top_k and row["weight_mode"] == weight_mode
            ]
            summary_rows.append({
                "top_k": top_k,
                "weight_mode": weight_mode,
                "condition_seed_units": len(block),
                "active_split_count": len({(row["dataset"], row["seed"]) for row in block if row["active"]}),
                "pooled_all_gene_mse_gain": 1.0 - sum(row["candidate_all_gene_mse"] for row in block) / max(sum(row["baseline_all_gene_mse"] for row in block), EPS),
                "pooled_top100_mse_gain": 1.0 - sum(row["candidate_top100_mse"] for row in block) / max(sum(row["baseline_top100_mse"] for row in block), EPS),
                "mean_top100_pcc_delta": np.mean([row["top100_pcc_delta"] for row in block]),
            })
    write_csv(args.out / "go_unseen_single_summary.csv", summary_rows)
    verdict = {
        "status": "PASS_GO_SHRINKAGE_ENDPOINT_PROBE",
        "selection_evidence": "go_identity_fit_by_split.csv only",
        "candidate_count": len(summary_rows),
        "test_rows_are_blind_audit_not_selection": True,
    }
    (args.out / "verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))
    for row in summary_rows:
        print(row)


if __name__ == "__main__":
    main()
