#!/usr/bin/env python3
"""Training-only multi-head screen for WitnessCell v14.

This probe never reads held-out combination means for fitting, routing, or
gating.  It uses official training singles and pseudo-holds out one known
single at a time.  The broad screen reuses one cross-fitted feature bank per
held endpoint; ``--nested`` performs an honest outer LOO in which the outer
endpoint is removed before every inner coefficient fit.

The dense background is always mean + GO residual and is fit on all genes.
Sparse heads are fit only on each training single's Welch top-100 residual.
"""
from __future__ import annotations

import argparse
import csv
import json
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from dual_head import go_program
from identity_head import EPS, one_sided_lower, pearson, self_statistic, welch_scores
from official_split import (
    filter_gears_supported_conditions,
    load_gears_supported_genes,
    make_official_split,
)


@dataclass(frozen=True)
class Bank:
    background: np.ndarray
    go_residual: np.ndarray
    self_anchor: np.ndarray
    go_sparse: Mapping[int, np.ndarray]
    response_sparse: Mapping[tuple[int, int], np.ndarray]
    coreg_sparse: Mapping[int, np.ndarray]
    coreg_go_sparse: Mapping[int, np.ndarray]
    response_coreg_sparse: Mapping[tuple[int, int], np.ndarray]
    go_amplitude_anchor: np.ndarray
    response_amplitude_anchor: Mapping[int, np.ndarray]
    fingerprint_tail_anchor: np.ndarray
    fingerprint_rms_anchor: np.ndarray


@dataclass(frozen=True)
class Recipe:
    name: str
    family: str
    support_k: int = 0
    neighbor_k: int = 0
    ridge: float = 0.0

    @property
    def incremental(self) -> bool:
        return self.family.startswith("inc_")


def top_mask(score: np.ndarray, size: int, exclude: int | None = None) -> np.ndarray:
    value = np.asarray(score, dtype=np.float64).copy()
    if exclude is not None:
        value[exclude] = 0.0
    order = np.argsort(-np.abs(value), kind="stable")[: min(size, len(value))]
    mask = np.zeros(len(value), dtype=np.float64)
    mask[order] = 1.0
    if exclude is not None:
        mask[exclude] = 0.0
    return mask


def safe_correlations(matrix: np.ndarray, target_index: int) -> np.ndarray:
    """Column correlations across observed perturbation rows."""
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    target = centered[:, target_index]
    denominator = np.sqrt(
        np.sum(centered * centered, axis=0) * np.sum(target * target)
    )
    return np.divide(
        centered.T @ target,
        denominator,
        out=np.zeros(matrix.shape[1], dtype=np.float64),
        where=denominator > EPS,
    )


def fingerprint_program(
    target: str,
    effects: Mapping[str, np.ndarray],
    gene_index: Mapping[str, int],
    neighbor_k: int,
) -> np.ndarray:
    """Stable kNN replacement for the rejected high-dimensional ridge.

    Endpoint similarity is computed from how the endpoint genes respond across
    the remaining known perturbations.  The privileged on-target coordinate
    in each known endpoint fingerprint is replaced by its off-diagonal mean.
    Only positive cosine neighbors contribute.
    """
    nodes = sorted(effects)
    matrix = np.stack([effects[node] for node in nodes]).astype(np.float64)
    background = matrix.mean(axis=0)
    if target not in gene_index:
        return background
    target_vector = matrix[:, gene_index[target]].copy()
    target_vector -= target_vector.mean()
    target_norm = float(np.linalg.norm(target_vector))
    similarities: list[tuple[float, str, np.ndarray]] = []
    for row_index, node in enumerate(nodes):
        if node not in gene_index:
            continue
        vector = matrix[:, gene_index[node]].copy()
        if len(vector) > 1:
            vector[row_index] = float(
                (vector.sum() - vector[row_index]) / (len(vector) - 1)
            )
        vector -= vector.mean()
        denominator = target_norm * float(np.linalg.norm(vector))
        similarity = float(target_vector @ vector / denominator) if denominator > EPS else 0.0
        similarities.append((max(similarity, 0.0), node, effects[node]))
    similarities.sort(key=lambda row: (row[0], row[1]), reverse=True)
    selected = similarities[: min(neighbor_k, len(similarities))]
    weights = np.asarray([row[0] for row in selected], dtype=np.float64)
    if weights.sum() <= EPS:
        return background
    return np.average(
        np.stack([row[2] for row in selected]), axis=0, weights=weights
    )


def build_bank(
    target: str,
    effects: Mapping[str, np.ndarray],
    genes: Sequence[str],
    gene2go: Mapping,
    go_top_k: int,
    support_sizes: Sequence[int],
    neighbor_sizes: Sequence[int],
) -> Bank:
    gene_index = {str(gene): index for index, gene in enumerate(genes)}
    matrix = np.stack(list(effects.values())).astype(np.float64)
    background = matrix.mean(axis=0)
    go_residual = go_program(target, effects, gene2go, go_top_k) - background
    target_index = gene_index.get(target)
    self_anchor = np.zeros(len(genes), dtype=np.float64)
    if target_index is not None:
        self_anchor[target_index] = self_statistic(effects, gene_index, "mean")

    go_amplitude_anchor = np.zeros(len(genes), dtype=np.float64)
    fingerprint_tail_anchor = np.zeros(len(genes), dtype=np.float64)
    fingerprint_rms_anchor = np.zeros(len(genes), dtype=np.float64)
    response_amplitude_anchor: dict[int, np.ndarray] = {}
    if target_index is not None:
        fingerprint = matrix[:, target_index]
        go_amplitude_anchor[target_index] = background[target_index] + go_residual[target_index]
        tail_quantile = 0.10 if self_anchor[target_index] < 0.0 else 0.90
        fingerprint_tail_anchor[target_index] = float(
            np.quantile(fingerprint, tail_quantile)
        )
        fingerprint_rms_anchor[target_index] = float(
            np.sign(self_anchor[target_index]) * np.sqrt(np.mean(np.square(fingerprint)))
        )

    correlations = (
        safe_correlations(matrix, target_index)
        if target_index is not None
        else np.zeros(len(genes), dtype=np.float64)
    )
    self_scale = self_statistic(effects, gene_index, "mean")
    coreg_direction = self_scale * correlations
    if target_index is not None:
        coreg_direction[target_index] = 0.0

    go_sparse: dict[int, np.ndarray] = {}
    coreg_sparse: dict[int, np.ndarray] = {}
    coreg_go_sparse: dict[int, np.ndarray] = {}
    response_sparse: dict[tuple[int, int], np.ndarray] = {}
    response_coreg_sparse: dict[tuple[int, int], np.ndarray] = {}
    for support_k in support_sizes:
        go_sparse[support_k] = go_residual * top_mask(
            go_residual, support_k, target_index
        )
        coreg_sparse[support_k] = coreg_direction * top_mask(
            correlations, support_k, target_index
        )
        consensus_score = np.abs(correlations) * np.sqrt(
            np.abs(go_residual) / (np.median(np.abs(go_residual)) + EPS) + EPS
        )
        coreg_go_sparse[support_k] = coreg_direction * top_mask(
            consensus_score, support_k, target_index
        )
        for neighbor_k in neighbor_sizes:
            response_program = fingerprint_program(
                target, effects, gene_index, neighbor_k
            )
            response_residual = response_program - background
            anchor = np.zeros(len(genes), dtype=np.float64)
            if target_index is not None:
                anchor[target_index] = response_program[target_index]
            response_amplitude_anchor[neighbor_k] = anchor
            response_sparse[(neighbor_k, support_k)] = response_residual * top_mask(
                response_residual, support_k, target_index
            )
            response_coreg_sparse[(neighbor_k, support_k)] = response_residual * top_mask(
                np.abs(response_residual) * (np.abs(correlations) + 0.05),
                support_k,
                target_index,
            )
    return Bank(
        background=background,
        go_residual=go_residual,
        self_anchor=self_anchor,
        go_sparse=go_sparse,
        response_sparse=response_sparse,
        coreg_sparse=coreg_sparse,
        coreg_go_sparse=coreg_go_sparse,
        response_coreg_sparse=response_coreg_sparse,
        go_amplitude_anchor=go_amplitude_anchor,
        response_amplitude_anchor=response_amplitude_anchor,
        fingerprint_tail_anchor=fingerprint_tail_anchor,
        fingerprint_rms_anchor=fingerprint_rms_anchor,
    )


def sparse_features(bank: Bank, recipe: Recipe) -> tuple[np.ndarray, ...]:
    if recipe.family == "self":
        return (bank.self_anchor,)
    if recipe.family == "go":
        return (bank.self_anchor, bank.go_sparse[recipe.support_k])
    if recipe.family == "response":
        return (
            bank.self_anchor,
            bank.response_sparse[(recipe.neighbor_k, recipe.support_k)],
        )
    if recipe.family == "coreg":
        return (bank.self_anchor, bank.coreg_sparse[recipe.support_k])
    if recipe.family == "coreg_go":
        return (bank.self_anchor, bank.coreg_go_sparse[recipe.support_k])
    if recipe.family == "response_coreg":
        return (
            bank.self_anchor,
            bank.response_coreg_sparse[(recipe.neighbor_k, recipe.support_k)],
        )
    if recipe.family == "amp_go":
        return (bank.self_anchor, bank.go_amplitude_anchor)
    if recipe.family == "amp_response":
        return (
            bank.self_anchor,
            bank.response_amplitude_anchor[recipe.neighbor_k],
        )
    if recipe.family == "amp_fingerprint":
        return (
            bank.self_anchor,
            bank.fingerprint_tail_anchor,
            bank.fingerprint_rms_anchor,
        )
    if recipe.family == "amp_joint":
        return (
            bank.self_anchor,
            bank.go_amplitude_anchor,
            bank.response_amplitude_anchor[recipe.neighbor_k],
            bank.fingerprint_tail_anchor,
            bank.fingerprint_rms_anchor,
        )
    if recipe.family == "amp_fp_go":
        return (
            bank.self_anchor,
            bank.fingerprint_tail_anchor,
            bank.fingerprint_rms_anchor,
            bank.go_sparse[recipe.support_k],
        )
    if recipe.family == "inc_amp_fp":
        return (
            bank.self_anchor,
            bank.fingerprint_tail_anchor,
            bank.fingerprint_rms_anchor,
        )
    if recipe.family == "inc_go":
        return (bank.go_sparse[recipe.support_k],)
    if recipe.family == "inc_amp_go":
        return (
            bank.self_anchor,
            bank.fingerprint_tail_anchor,
            bank.fingerprint_rms_anchor,
            bank.go_sparse[recipe.support_k],
        )
    if recipe.family == "inc_response":
        return (bank.response_sparse[(recipe.neighbor_k, recipe.support_k)],)
    raise ValueError(recipe.family)


def solve(
    feature_rows: Sequence[tuple[np.ndarray, ...]],
    targets: Sequence[np.ndarray],
    indices: Sequence[np.ndarray] | None,
    lower: float,
    upper: float,
    ridge: float = 1e-6,
    relative_ridge: float = 0.0,
) -> np.ndarray:
    width = len(feature_rows[0])
    gram = np.zeros((width, width), dtype=np.float64)
    cross = np.zeros(width, dtype=np.float64)
    for row_index, (features, target) in enumerate(zip(feature_rows, targets, strict=True)):
        design = np.stack(features, axis=1).astype(np.float64)
        response = np.asarray(target, dtype=np.float64)
        if indices is not None:
            keep = indices[row_index]
            design = design[keep]
            response = response[keep]
        gram += design.T @ design
        cross += design.T @ response
    penalty = ridge * np.eye(width)
    if relative_ridge > 0.0:
        penalty += relative_ridge * np.diag(np.maximum(np.diag(gram), EPS))
    weight = np.linalg.solve(gram + penalty, cross)
    return np.clip(weight, lower, upper)


def dense_fit(banks: Sequence[Bank], targets: Sequence[np.ndarray]) -> np.ndarray:
    return solve(
        [(bank.background, bank.go_residual) for bank in banks],
        targets,
        indices=None,
        lower=0.0,
        upper=2.0,
    )


def dense_predict(bank: Bank, weight: np.ndarray) -> np.ndarray:
    return weight[0] * bank.background + weight[1] * bank.go_residual


def v13_fit(banks: Sequence[Bank], targets: Sequence[np.ndarray]) -> np.ndarray:
    """Frozen v13 joint mean+GO+self all-gene fit under the same LOO view."""
    return solve(
        [(bank.background, bank.go_residual, bank.self_anchor) for bank in banks],
        targets,
        indices=None,
        lower=0.0,
        upper=2.0,
    )


def v13_predict(bank: Bank, weight: np.ndarray) -> np.ndarray:
    return (
        weight[0] * bank.background
        + weight[1] * bank.go_residual
        + weight[2] * bank.self_anchor
    )


def sparse_fit(
    banks: Sequence[Bank],
    targets: Sequence[np.ndarray],
    top_indices: Sequence[np.ndarray],
    dense_weight: np.ndarray,
    recipe: Recipe,
    frozen_v13_weight: np.ndarray,
) -> np.ndarray:
    if recipe.incremental:
        residuals = [
            target - v13_predict(bank, frozen_v13_weight)
            for bank, target in zip(banks, targets, strict=True)
        ]
        lower, upper = -1.5, 1.5
    else:
        residuals = [
            target - dense_predict(bank, dense_weight)
            for bank, target in zip(banks, targets, strict=True)
        ]
        lower, upper = 0.0, 2.5
    return solve(
        [sparse_features(bank, recipe) for bank in banks],
        residuals,
        indices=top_indices,
        lower=lower,
        upper=upper,
        ridge=1e-4,
        relative_ridge=recipe.ridge,
    )


def recipes(support_sizes: Sequence[int], neighbor_sizes: Sequence[int]) -> list[Recipe]:
    result = [
        Recipe("self", "self"),
        Recipe("amp_go", "amp_go"),
        Recipe("amp_fingerprint", "amp_fingerprint"),
    ]
    for neighbor_k in neighbor_sizes:
        result.extend([
            Recipe(f"amp_response_n{neighbor_k}", "amp_response", neighbor_k=neighbor_k),
            Recipe(f"amp_joint_n{neighbor_k}", "amp_joint", neighbor_k=neighbor_k),
        ])
    for ridge_name, ridge in (("r001", 0.01), ("r01", 0.1), ("r1", 1.0), ("r10", 10.0)):
        result.append(Recipe(f"amp_fp_{ridge_name}", "amp_fingerprint", ridge=ridge))
    for ridge_name, ridge in (
        ("r003", 0.03), ("r005", 0.05), ("r01", 0.1),
        ("r02", 0.2), ("r03", 0.3), ("r05", 0.5),
        ("r1", 1.0), ("r3", 3.0), ("r10", 10.0), ("r30", 30.0),
    ):
        result.append(Recipe(f"inc_amp_fp_{ridge_name}", "inc_amp_fp", ridge=ridge))
    for support_k in support_sizes:
        result.extend([
            Recipe(f"go_s{support_k}", "go", support_k),
            Recipe(f"coreg_s{support_k}", "coreg", support_k),
            Recipe(f"coreg_go_s{support_k}", "coreg_go", support_k),
        ])
        for ridge_name, ridge in (("r001", 0.01), ("r01", 0.1), ("r1", 1.0)):
            result.append(Recipe(
                f"amp_fp_go_s{support_k}_{ridge_name}",
                "amp_fp_go",
                support_k=support_k,
                ridge=ridge,
            ))
        for ridge_name, ridge in (("r01", 0.1), ("r1", 1.0), ("r3", 3.0), ("r10", 10.0), ("r30", 30.0)):
            result.extend([
                Recipe(
                    f"inc_go_s{support_k}_{ridge_name}", "inc_go",
                    support_k=support_k, ridge=ridge,
                ),
                Recipe(
                    f"inc_amp_go_s{support_k}_{ridge_name}", "inc_amp_go",
                    support_k=support_k, ridge=ridge,
                ),
            ])
        for neighbor_k in neighbor_sizes:
            for ridge_name, ridge in (("r1", 1.0), ("r3", 3.0), ("r10", 10.0)):
                result.append(Recipe(
                    f"inc_response_n{neighbor_k}_s{support_k}_{ridge_name}",
                    "inc_response", support_k=support_k,
                    neighbor_k=neighbor_k, ridge=ridge,
                ))
        for neighbor_k in neighbor_sizes:
            result.extend([
                Recipe(
                    f"response_n{neighbor_k}_s{support_k}",
                    "response",
                    support_k,
                    neighbor_k,
                ),
                Recipe(
                    f"response_coreg_n{neighbor_k}_s{support_k}",
                    "response_coreg",
                    support_k,
                    neighbor_k,
                ),
            ])
    return result


def top_indices_for(
    node: str,
    means: Mapping[str, np.ndarray],
    variances: Mapping[str, np.ndarray],
    counts: Mapping[str, int],
) -> np.ndarray:
    score = welch_scores(
        means[node], variances[node], counts[node],
        means["control"], variances["control"], counts["control"],
    )
    return np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]


def evaluate_one(
    truth: np.ndarray,
    top: np.ndarray,
    bank: Bank,
    dense_weight: np.ndarray,
    sparse_weight: np.ndarray,
    recipe: Recipe,
    v13_weight: np.ndarray,
) -> dict:
    baseline = bank.background
    dense = dense_predict(bank, dense_weight)
    sparse = sum(
        weight * feature
        for weight, feature in zip(sparse_weight, sparse_features(bank, recipe), strict=True)
    )
    v13 = v13_predict(bank, v13_weight)
    final = v13 + sparse if recipe.incremental else dense + sparse
    base_all = float(np.mean(np.square(baseline - truth)))
    dense_all = float(np.mean(np.square(dense - truth)))
    final_all = float(np.mean(np.square(final - truth)))
    v13_all = float(np.mean(np.square(v13 - truth)))
    base_top = float(np.mean(np.square(baseline[top] - truth[top])))
    dense_top = float(np.mean(np.square(dense[top] - truth[top])))
    final_top = float(np.mean(np.square(final[top] - truth[top])))
    v13_top = float(np.mean(np.square(v13[top] - truth[top])))
    v13_pcc = pearson(v13[top], truth[top])
    final_pcc = pearson(final[top], truth[top])
    return {
        "baseline_all_gene_mse": base_all,
        "dense_all_gene_mse": dense_all,
        "final_all_gene_mse": final_all,
        "v13_all_gene_mse": v13_all,
        "final_all_gene_gain": (base_all - final_all) / max(base_all, EPS),
        "upgrade_all_gene_gain": (v13_all - final_all) / max(v13_all, EPS),
        "baseline_top100_mse": base_top,
        "dense_top100_mse": dense_top,
        "final_top100_mse": final_top,
        "v13_top100_mse": v13_top,
        "final_top100_gain": (base_top - final_top) / max(base_top, EPS),
        "upgrade_top100_gain": (v13_top - final_top) / max(v13_top, EPS),
        "incremental_top100_gain": (dense_top - final_top) / max(dense_top, EPS),
        "baseline_top100_pcc": pearson(baseline[top], truth[top]),
        "dense_top100_pcc": pearson(dense[top], truth[top]),
        "v13_top100_pcc": v13_pcc,
        "final_top100_pcc": final_pcc,
        "final_top100_pcc_delta": final_pcc - pearson(baseline[top], truth[top]),
        "incremental_top100_pcc_delta": final_pcc - pearson(dense[top], truth[top]),
        "upgrade_top100_pcc_delta": final_pcc - v13_pcc,
    }


def screen_split(payload: dict, recipe_list: Sequence[Recipe], nested: bool) -> list[dict]:
    nodes = payload["nodes"]
    effects = payload["effects"]
    genes = payload["genes"]
    gene2go = payload["gene2go"]
    support_sizes = payload["support_sizes"]
    neighbor_sizes = payload["neighbor_sizes"]
    rows: list[dict] = []

    if not nested:
        banks = []
        truths = []
        tops = []
        for held in nodes:
            remaining = {node: effects[node] for node in nodes if node != held}
            banks.append(build_bank(
                held, remaining, genes, gene2go, 10, support_sizes, neighbor_sizes
            ))
            truths.append(effects[held])
            tops.append(payload["tops"][held])
        dense_weight = dense_fit(banks, truths)
        frozen_v13_weight = v13_fit(banks, truths)
        for recipe in recipe_list:
            sparse_weight = sparse_fit(
                banks, truths, tops, dense_weight, recipe, frozen_v13_weight
            )
            for held, bank, truth, top in zip(nodes, banks, truths, tops, strict=True):
                rows.append({
                    "condition": held,
                    "recipe": recipe.name,
                    "dense_weights": ";".join(map(str, dense_weight)),
                    "sparse_weights": ";".join(map(str, sparse_weight)),
                    **evaluate_one(
                        truth, top, bank, dense_weight, sparse_weight, recipe,
                        frozen_v13_weight,
                    ),
                })
        return rows

    # Honest outer LOO: the outer target is unavailable even inside the inner
    # coefficient-fitting feature rows.
    for outer in nodes:
        remaining_nodes = [node for node in nodes if node != outer]
        inner_banks = []
        inner_truths = []
        inner_tops = []
        for inner in remaining_nodes:
            inner_effects = {
                node: effects[node]
                for node in remaining_nodes
                if node != inner
            }
            inner_banks.append(build_bank(
                inner, inner_effects, genes, gene2go, 10,
                support_sizes, neighbor_sizes,
            ))
            inner_truths.append(effects[inner])
            inner_tops.append(payload["tops"][inner])
        outer_effects = {node: effects[node] for node in remaining_nodes}
        outer_bank = build_bank(
            outer, outer_effects, genes, gene2go, 10,
            support_sizes, neighbor_sizes,
        )
        dense_weight = dense_fit(inner_banks, inner_truths)
        frozen_v13_weight = v13_fit(inner_banks, inner_truths)
        for recipe in recipe_list:
            sparse_weight = sparse_fit(
                inner_banks, inner_truths, inner_tops, dense_weight, recipe,
                frozen_v13_weight,
            )
            rows.append({
                "condition": outer,
                "recipe": recipe.name,
                "dense_weights": ";".join(map(str, dense_weight)),
                "sparse_weights": ";".join(map(str, sparse_weight)),
                **evaluate_one(
                    effects[outer], payload["tops"][outer], outer_bank,
                    dense_weight, sparse_weight, recipe, frozen_v13_weight,
                ),
            })
    return rows


def summarize(rows: Sequence[dict]) -> list[dict]:
    groups: dict[tuple[str, str, int], list[dict]] = {}
    for row in rows:
        groups.setdefault((row["dataset"], row["recipe"], row["seed"]), []).append(row)
    result = []
    for (dataset, recipe, seed), block in sorted(groups.items()):
        def values(key: str) -> np.ndarray:
            return np.asarray([row[key] for row in block], dtype=np.float64)
        base_all = values("baseline_all_gene_mse").sum()
        final_all = values("final_all_gene_mse").sum()
        v13_all = values("v13_all_gene_mse").sum()
        base_top = values("baseline_top100_mse").sum()
        final_top = values("final_top100_mse").sum()
        v13_top = values("v13_top100_mse").sum()
        result.append({
            "dataset": dataset,
            "seed": seed,
            "recipe": recipe,
            "condition_count": len(block),
            "pooled_all_gene_gain": 1.0 - final_all / max(base_all, EPS),
            "all_gene_gain_lower": one_sided_lower(values("final_all_gene_gain"), 0.95),
            "pooled_upgrade_all_gene_gain": 1.0 - final_all / max(v13_all, EPS),
            "upgrade_all_gene_gain_lower": one_sided_lower(values("upgrade_all_gene_gain"), 0.95),
            "pooled_top100_gain": 1.0 - final_top / max(base_top, EPS),
            "top100_gain_lower": one_sided_lower(values("final_top100_gain"), 0.95),
            "pooled_upgrade_top100_gain": 1.0 - final_top / max(v13_top, EPS),
            "upgrade_top100_gain_lower": one_sided_lower(values("upgrade_top100_gain"), 0.95),
            "mean_top100_pcc_delta": float(values("final_top100_pcc_delta").mean()),
            "top100_pcc_delta_lower": one_sided_lower(values("final_top100_pcc_delta"), 0.95),
            "mean_upgrade_top100_pcc_delta": float(values("upgrade_top100_pcc_delta").mean()),
            "upgrade_top100_pcc_delta_lower": one_sided_lower(values("upgrade_top100_pcc_delta"), 0.95),
            "mean_incremental_top100_gain": float(values("incremental_top100_gain").mean()),
            "incremental_top100_gain_lower": one_sided_lower(values("incremental_top100_gain"), 0.95),
            "mean_incremental_top100_pcc_delta": float(values("incremental_top100_pcc_delta").mean()),
            "incremental_top100_pcc_delta_lower": one_sided_lower(values("incremental_top100_pcc_delta"), 0.95),
        })
    return result


def aggregate(summary: Sequence[dict]) -> list[dict]:
    recipes_found = sorted({row["recipe"] for row in summary})
    result = []
    for recipe in recipes_found:
        block = [row for row in summary if row["recipe"] == recipe]
        result.append({
            "recipe": recipe,
            "split_count": len(block),
            "mean_all_gene_gain": float(np.mean([row["pooled_all_gene_gain"] for row in block])),
            "min_all_gene_lcb": float(min(row["all_gene_gain_lower"] for row in block)),
            "mean_upgrade_all_gene_gain": float(np.mean([row["pooled_upgrade_all_gene_gain"] for row in block])),
            "min_upgrade_all_gene_lcb": float(min(row["upgrade_all_gene_gain_lower"] for row in block)),
            "mean_top100_gain": float(np.mean([row["pooled_top100_gain"] for row in block])),
            "min_top100_lcb": float(min(row["top100_gain_lower"] for row in block)),
            "mean_upgrade_top100_gain": float(np.mean([row["pooled_upgrade_top100_gain"] for row in block])),
            "min_upgrade_top100_gain_lcb": float(min(row["upgrade_top100_gain_lower"] for row in block)),
            "mean_top100_pcc_delta": float(np.mean([row["mean_top100_pcc_delta"] for row in block])),
            "min_top100_pcc_lcb": float(min(row["top100_pcc_delta_lower"] for row in block)),
            "mean_upgrade_top100_pcc_delta": float(np.mean([row["mean_upgrade_top100_pcc_delta"] for row in block])),
            "min_upgrade_top100_pcc_lcb": float(min(row["upgrade_top100_pcc_delta_lower"] for row in block)),
            "mean_incremental_top100_gain": float(np.mean([row["mean_incremental_top100_gain"] for row in block])),
            "min_incremental_top100_gain_lcb": float(min(row["incremental_top100_gain_lower"] for row in block)),
            "mean_incremental_top100_pcc_delta": float(np.mean([row["mean_incremental_top100_pcc_delta"] for row in block])),
            "min_incremental_top100_pcc_lcb": float(min(row["incremental_top100_pcc_delta_lower"] for row in block)),
        })
    return sorted(
        result,
        key=lambda row: (
            -row["mean_top100_gain"],
            -row["mean_top100_pcc_delta"],
            -row["mean_all_gene_gain"],
            row["recipe"],
        ),
    )


def write_csv(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, action="append", required=True)
    parser.add_argument("--gears-assets", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=(1,))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--nested", action="store_true")
    parser.add_argument("--recipes", nargs="*", default=())
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    support_sizes = (5, 10, 20, 40)
    neighbor_sizes = (3, 5, 10)
    all_recipes = recipes(support_sizes, neighbor_sizes)
    if args.recipes:
        requested = set(args.recipes)
        all_recipes = [recipe for recipe in all_recipes if recipe.name in requested]
        missing = requested - {recipe.name for recipe in all_recipes}
        if missing:
            raise ValueError(f"unknown recipes: {sorted(missing)}")

    with (args.gears_assets / "gene2go_all.pkl").open("rb") as stream:
        gene2go = pickle.load(stream)
    supported = load_gears_supported_genes(args.gears_assets)
    detail_rows: list[dict] = []
    for cache_path in args.cache:
        data = np.load(cache_path, allow_pickle=False)
        dataset = str(data["dataset"])
        conditions = data["conditions"].astype(str).tolist()
        genes = data["genes"].astype(str).tolist()
        means = dict(zip(conditions, data["means"].astype(np.float64)))
        variances = dict(zip(conditions, data["variances"].astype(np.float64)))
        counts = dict(zip(conditions, map(int, data["counts"].astype(np.int64))))
        official = filter_gears_supported_conditions(conditions, supported)
        for seed in args.seeds:
            split = make_official_split(official, seed)
            effects = {
                condition: means[condition] - means["control"]
                for condition in split.train
                if condition != "control" and "+" not in condition
            }
            nodes = sorted(node for node in effects if node in set(genes))
            payload = {
                "dataset": dataset,
                "seed": seed,
                "nodes": nodes,
                "effects": effects,
                "means": means,
                "variances": variances,
                "counts": counts,
                "genes": genes,
                "gene2go": gene2go,
                "support_sizes": support_sizes,
                "neighbor_sizes": neighbor_sizes,
                "tops": {
                    node: top_indices_for(node, means, variances, counts)
                    for node in nodes
                },
            }
            block = screen_split(payload, all_recipes, args.nested)
            for row in block:
                detail_rows.append({
                    "dataset": dataset,
                    "seed": seed,
                    "nested": int(args.nested),
                    **row,
                })

    summary_rows = summarize(detail_rows)
    aggregate_rows = aggregate(summary_rows)
    write_csv(args.out / "by_condition.csv", detail_rows)
    write_csv(args.out / "by_split.csv", summary_rows)
    write_csv(args.out / "aggregate.csv", aggregate_rows)
    verdict = {
        "status": "PASS_TRAINING_ONLY_MULTHEAD_PROBE",
        "nested_outer_loo": args.nested,
        "test_means_used_for_fit_or_selection": False,
        "recipe_count": len(all_recipes),
        "split_count": len({(row["dataset"], row["seed"]) for row in detail_rows}),
        "top_candidates": aggregate_rows[:10],
    }
    (args.out / "verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
