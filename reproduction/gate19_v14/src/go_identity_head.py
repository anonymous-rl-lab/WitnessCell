#!/usr/bin/env python3
"""GO-shrunk endpoint program with a training-LOO self-response anchor."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from identity_head import EPS, one_sided_lower, pearson, self_statistic, welch_scores


def jaccard(left: set, right: set) -> float:
    return len(left & right) / max(len(left | right), 1)


def features_for_node(
    node: str,
    effects: Mapping[str, np.ndarray],
    gene_index: Mapping[str, int],
    gene2go: Mapping,
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
    self_anchor = np.zeros_like(background)
    if node in gene_index:
        self_anchor[gene_index[node]] = self_statistic(
            effects, gene_index, "mean"
        )
    return background, go_residual, self_anchor


def closed_form_weights(
    rows: Sequence[tuple[np.ndarray, np.ndarray, np.ndarray]],
    targets: Sequence[np.ndarray],
) -> tuple[float, float, float]:
    gram = np.zeros((3, 3), dtype=np.float64)
    cross = np.zeros(3, dtype=np.float64)
    for row, target in zip(rows, targets, strict=True):
        design = np.stack(row, axis=1).astype(np.float64)
        gram += design.T @ design
        cross += design.T @ np.asarray(target, dtype=np.float64)
    weight = np.linalg.solve(gram + 1e-8 * np.eye(3), cross)
    weight = np.clip(weight, [0.0, 0.0, 0.0], [2.0, 2.0, 2.0])
    return float(weight[0]), float(weight[1]), float(weight[2])


@dataclass(frozen=True)
class GOIdentityHeadFit:
    active: bool
    background_weight: float
    go_residual_weight: float
    self_weight: float
    top_k: int
    gate_metric: str
    gate_lower_bound: float
    loo_all_gene_gain_mean: float
    loo_all_gene_gain_lower: float
    loo_top100_gain_mean: float
    loo_top100_gain_lower: float
    loo_pcc_delta_mean: float
    known_single_count: int
    records: tuple[dict, ...]

    def predict(
        self,
        node: str,
        effects: Mapping[str, np.ndarray],
        gene_index: Mapping[str, int],
        gene2go: Mapping,
    ) -> np.ndarray:
        background, go_residual, self_anchor = features_for_node(
            node, effects, gene_index, gene2go, self.top_k
        )
        if not self.active:
            return background
        return (
            self.background_weight * background
            + self.go_residual_weight * go_residual
            + self.self_weight * self_anchor
        )


def fit_go_identity_head(
    single_effects: Mapping[str, np.ndarray],
    single_means: Mapping[str, np.ndarray],
    single_variances: Mapping[str, np.ndarray],
    single_counts: Mapping[str, int],
    control_mean: np.ndarray,
    control_variance: np.ndarray,
    control_count: int,
    genes: Sequence[str],
    gene2go: Mapping,
    top_k: int = 40,
    gate_confidence: float = 0.95,
) -> GOIdentityHeadFit:
    gene_index = {str(gene): index for index, gene in enumerate(genes)}
    nodes = sorted(node for node in single_effects if node in gene_index)
    if len(nodes) < 3:
        return GOIdentityHeadFit(
            active=False,
            background_weight=1.0,
            go_residual_weight=0.0,
            self_weight=0.0,
            top_k=top_k,
            gate_metric="all_gene_mse",
            gate_lower_bound=float("-inf"),
            loo_all_gene_gain_mean=0.0,
            loo_all_gene_gain_lower=float("-inf"),
            loo_top100_gain_mean=0.0,
            loo_top100_gain_lower=float("-inf"),
            loo_pcc_delta_mean=0.0,
            known_single_count=len(nodes),
            records=(),
        )

    rows = []
    targets = []
    for held in nodes:
        remaining = {
            node: single_effects[node] for node in nodes if node != held
        }
        rows.append(features_for_node(held, remaining, gene_index, gene2go, top_k))
        targets.append(single_effects[held])
    background_weight, go_weight, self_weight = closed_form_weights(rows, targets)
    records = []
    for held, row, truth in zip(nodes, rows, targets, strict=True):
        background = row[0]
        prediction = (
            background_weight * row[0]
            + go_weight * row[1]
            + self_weight * row[2]
        )
        score = welch_scores(
            single_means[held], single_variances[held], single_counts[held],
            control_mean, control_variance, control_count,
        )
        order = np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]
        baseline_all = float(np.mean(np.square(background - truth)))
        enhanced_all = float(np.mean(np.square(prediction - truth)))
        baseline_top = float(np.mean(np.square(background[order] - truth[order])))
        enhanced_top = float(np.mean(np.square(prediction[order] - truth[order])))
        baseline_pcc = pearson(background[order], truth[order])
        enhanced_pcc = pearson(prediction[order], truth[order])
        records.append({
            "condition": held,
            "baseline_all_gene_mse": baseline_all,
            "enhanced_all_gene_mse": enhanced_all,
            "relative_all_gene_mse_gain": (baseline_all - enhanced_all) / max(baseline_all, EPS),
            "baseline_top100_mse": baseline_top,
            "enhanced_top100_mse": enhanced_top,
            "relative_top100_mse_gain": (baseline_top - enhanced_top) / max(baseline_top, EPS),
            "baseline_top100_pcc": baseline_pcc,
            "enhanced_top100_pcc": enhanced_pcc,
            "top100_pcc_delta": enhanced_pcc - baseline_pcc,
        })
    all_gain = np.asarray([row["relative_all_gene_mse_gain"] for row in records])
    top_gain = np.asarray([row["relative_top100_mse_gain"] for row in records])
    pcc_delta = np.asarray([row["top100_pcc_delta"] for row in records])
    all_lower = one_sided_lower(all_gain, gate_confidence)
    top_lower = one_sided_lower(top_gain, gate_confidence)
    return GOIdentityHeadFit(
        active=bool(all_lower > 0.0),
        background_weight=background_weight,
        go_residual_weight=go_weight,
        self_weight=self_weight,
        top_k=top_k,
        gate_metric="all_gene_mse",
        gate_lower_bound=all_lower,
        loo_all_gene_gain_mean=float(all_gain.mean()),
        loo_all_gene_gain_lower=all_lower,
        loo_top100_gain_mean=float(top_gain.mean()),
        loo_top100_gain_lower=top_lower,
        loo_pcc_delta_mean=float(pcc_delta.mean()),
        known_single_count=len(nodes),
        records=tuple(records),
    )
