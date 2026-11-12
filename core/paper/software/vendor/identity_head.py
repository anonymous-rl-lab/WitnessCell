#!/usr/bin/env python3
"""Training-single-only endpoint identity head for WitnessCell.

For an unseen perturbation endpoint, the v1 fallback is the mean known-single
effect. The new head decomposes that fallback into a background program and a
sparse on-target self-response anchor. All weights and the activation gate are
estimated from leave-one-known-single-out predictions; target outcomes are
never used.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from scipy import stats


EPS = 1e-12


def welch_scores(
    x_mean: np.ndarray,
    x_var: np.ndarray,
    x_n: int,
    y_mean: np.ndarray,
    y_var: np.ndarray,
    y_n: int,
) -> np.ndarray:
    x_unbiased = x_var * x_n / max(x_n - 1, 1)
    y_unbiased = y_var * y_n / max(y_n - 1, 1)
    denominator = np.sqrt(
        x_unbiased / max(x_n, 1) + y_unbiased / max(y_n, 1)
    )
    return (x_mean - y_mean) / np.maximum(denominator, EPS)


def pearson(left: np.ndarray, right: np.ndarray) -> float:
    left = left - left.mean()
    right = right - right.mean()
    return float(np.sum(left * right) / max(np.linalg.norm(left) * np.linalg.norm(right), EPS))


def one_sided_lower(values: np.ndarray, confidence: float = 0.95) -> float:
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return float("-inf")
    standard_error = float(values.std(ddof=1) / np.sqrt(len(values)))
    critical = float(stats.t.ppf(confidence, len(values) - 1))
    return float(values.mean() - critical * standard_error)


def self_statistic(
    effects: Mapping[str, np.ndarray],
    gene_index: Mapping[str, int],
    statistic: str,
) -> float:
    values = [
        float(effect[gene_index[node]])
        for node, effect in effects.items()
        if node in gene_index
    ]
    if not values:
        return 0.0
    if statistic == "mean":
        return float(np.mean(values))
    if statistic == "median":
        return float(np.median(values))
    raise ValueError(f"unknown self statistic: {statistic}")


def features_for_node(
    node: str,
    effects: Mapping[str, np.ndarray],
    gene_index: Mapping[str, int],
    anchor_mode: str,
    self_summary: str,
) -> tuple[np.ndarray, np.ndarray]:
    background = np.mean(np.stack(list(effects.values())), axis=0)
    anchor = np.zeros_like(background)
    if node in gene_index:
        index = gene_index[node]
        self_value = self_statistic(effects, gene_index, self_summary)
        if anchor_mode == "add":
            anchor[index] = self_value
        elif anchor_mode == "replace":
            anchor[index] = self_value - background[index]
        else:
            raise ValueError(f"unknown anchor mode: {anchor_mode}")
    return background, anchor


def closed_form_weights(
    feature_rows: Sequence[tuple[np.ndarray, np.ndarray]],
    targets: Sequence[np.ndarray],
    ridge: float = 1e-8,
) -> tuple[float, float]:
    gram = np.zeros((2, 2), dtype=np.float64)
    cross = np.zeros(2, dtype=np.float64)
    for (background, anchor), target in zip(feature_rows, targets, strict=True):
        x1 = np.asarray(background, dtype=np.float64)
        x2 = np.asarray(anchor, dtype=np.float64)
        y = np.asarray(target, dtype=np.float64)
        gram[0, 0] += np.dot(x1, x1)
        gram[0, 1] += np.dot(x1, x2)
        gram[1, 0] += np.dot(x2, x1)
        gram[1, 1] += np.dot(x2, x2)
        cross[0] += np.dot(x1, y)
        cross[1] += np.dot(x2, y)
    gram += ridge * np.eye(2)
    weights = np.linalg.solve(gram, cross)
    return float(np.clip(weights[0], 0.0, 2.0)), float(np.clip(weights[1], 0.0, 2.0))


@dataclass(frozen=True)
class IdentityHeadFit:
    active: bool
    background_weight: float
    self_weight: float
    anchor_mode: str
    self_summary: str
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
    ) -> np.ndarray:
        background, anchor = features_for_node(
            node, effects, gene_index, self.anchor_mode, self.self_summary
        )
        if not self.active:
            return background
        return self.background_weight * background + self.self_weight * anchor


def fit_identity_head(
    single_effects: Mapping[str, np.ndarray],
    single_means: Mapping[str, np.ndarray],
    single_variances: Mapping[str, np.ndarray],
    single_counts: Mapping[str, int],
    control_mean: np.ndarray,
    control_variance: np.ndarray,
    control_count: int,
    genes: Sequence[str],
    anchor_mode: str = "add",
    self_summary: str = "mean",
    weight_mode: str = "closed_form",
    gate_metric: str = "top100_mse",
    gate_confidence: float = 0.95,
) -> IdentityHeadFit:
    gene_index = {str(gene): index for index, gene in enumerate(genes)}
    nodes = sorted(node for node in single_effects if node in gene_index)
    if len(nodes) < 3:
        return IdentityHeadFit(
            active=False,
            background_weight=1.0,
            self_weight=0.0,
            anchor_mode=anchor_mode,
            self_summary=self_summary,
            gate_metric=gate_metric,
            gate_lower_bound=float("-inf"),
            loo_all_gene_gain_mean=0.0,
            loo_all_gene_gain_lower=float("-inf"),
            loo_top100_gain_mean=0.0,
            loo_top100_gain_lower=float("-inf"),
            loo_pcc_delta_mean=0.0,
            known_single_count=len(nodes),
            records=(),
        )

    feature_rows = []
    targets = []
    held_rows = []
    for held in nodes:
        remaining = {node: single_effects[node] for node in nodes if node != held}
        features = features_for_node(
            held, remaining, gene_index, anchor_mode, self_summary
        )
        feature_rows.append(features)
        targets.append(single_effects[held])
        held_rows.append((held, remaining, features))

    if weight_mode == "fixed":
        background_weight, self_weight = 1.0, 1.0
    elif weight_mode == "closed_form":
        background_weight, self_weight = closed_form_weights(feature_rows, targets)
    else:
        raise ValueError(f"unknown weight mode: {weight_mode}")

    records = []
    for held, _, (background, anchor) in held_rows:
        prediction = background_weight * background + self_weight * anchor
        truth = single_effects[held]
        baseline_mse = float(np.mean(np.square(background - truth)))
        enhanced_mse = float(np.mean(np.square(prediction - truth)))
        score = welch_scores(
            single_means[held],
            single_variances[held],
            single_counts[held],
            control_mean,
            control_variance,
            control_count,
        )
        order = np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]
        baseline_top = float(np.mean(np.square(background[order] - truth[order])))
        enhanced_top = float(np.mean(np.square(prediction[order] - truth[order])))
        baseline_pcc = pearson(background[order], truth[order])
        enhanced_pcc = pearson(prediction[order], truth[order])
        records.append({
            "condition": held,
            "baseline_all_gene_mse": baseline_mse,
            "enhanced_all_gene_mse": enhanced_mse,
            "relative_all_gene_mse_gain": (baseline_mse - enhanced_mse) / max(baseline_mse, EPS),
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
    lower = top_lower if gate_metric == "top100_mse" else all_lower
    return IdentityHeadFit(
        active=bool(lower > 0.0),
        background_weight=background_weight,
        self_weight=self_weight,
        anchor_mode=anchor_mode,
        self_summary=self_summary,
        gate_metric=gate_metric,
        gate_lower_bound=lower,
        loo_all_gene_gain_mean=float(all_gain.mean()),
        loo_all_gene_gain_lower=all_lower,
        loo_top100_gain_mean=float(top_gain.mean()),
        loo_top100_gain_lower=top_lower,
        loo_pcc_delta_mean=float(pcc_delta.mean()),
        known_single_count=len(nodes),
        records=tuple(records),
    )
