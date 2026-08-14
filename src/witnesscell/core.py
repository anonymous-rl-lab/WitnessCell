"""Frozen WitnessCell v14 numerical core.

This module is a productized, dependency-light transcription of the Gate 19
reference implementation.  The default path intentionally preserves the
published training-only evidence gates, factorized backbone, endpoint-
incidence interaction witness, and exact fallback behavior.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import stats

EPS = 1e-12
DEFAULT_ALPHA_GRID = (0.0, 0.02, 0.05, 0.1, 0.2, 0.5)
DEFAULT_NOISE_GRID = (0.01, 0.03, 0.1, 0.3, 1.0)


def parse_condition(condition: str, separator: str = "+") -> tuple[str, ...]:
    """Return one or two perturbation endpoints from a condition label."""
    if not isinstance(condition, str) or not condition.strip():
        raise ValueError("condition labels must be non-empty strings")
    parts = tuple(part.strip() for part in condition.split(separator))
    if any(not part for part in parts):
        raise ValueError(f"invalid condition label: {condition!r}")
    if len(parts) > 2:
        raise ValueError(
            "WitnessCell 0.1 supports singles and two-endpoint combinations only"
        )
    if len(set(parts)) != len(parts):
        raise ValueError("a combination must contain two distinct endpoints")
    return parts


def saturate(effect: np.ndarray, strength: float) -> np.ndarray:
    """Apply the frozen element-wise saturating factorized map."""
    value = np.asarray(effect, dtype=float)
    return value / (1.0 + float(strength) * np.abs(value))


def closed_form_gamma(prediction: np.ndarray, truth: np.ndarray) -> float:
    """Return the validation-calibrated residual amplitude in ``[0, 1]``."""
    prediction = np.asarray(prediction, dtype=float)
    truth = np.asarray(truth, dtype=float)
    denominator = float(np.sum(prediction * prediction))
    if denominator <= EPS:
        return 0.0
    return float(np.clip(np.sum(prediction * truth) / denominator, 0.0, 1.0))


def incidence(conditions: Sequence[str], nodes: Sequence[str]) -> np.ndarray:
    """Build the unsigned endpoint-incidence design used by Gate 19."""
    node_id = {node: index for index, node in enumerate(nodes)}
    result = np.zeros((len(conditions), len(nodes)), dtype=float)
    for row, condition in enumerate(conditions):
        for node in parse_condition(condition):
            result[row, node_id[node]] = 1.0
    return result


def kernel_predict(
    train_design: np.ndarray,
    target_design: np.ndarray,
    response: np.ndarray,
    noise_ratio: float,
) -> np.ndarray:
    """Transfer residuals through the frozen endpoint-incidence kernel."""
    train_kernel = train_design @ train_design.T / 2.0
    cross_kernel = train_design @ target_design.T / 2.0
    covariance = train_kernel + float(noise_ratio) * np.eye(len(train_kernel))
    return np.linalg.solve(covariance, cross_kernel).T @ response


def welch_scores(
    x_mean: np.ndarray,
    x_var: np.ndarray,
    x_n: int,
    y_mean: np.ndarray,
    y_var: np.ndarray,
    y_n: int,
) -> np.ndarray:
    """Compute the frozen absolute-Welch ranking score."""
    x_unbiased = np.asarray(x_var, float) * x_n / max(x_n - 1, 1)
    y_unbiased = np.asarray(y_var, float) * y_n / max(y_n - 1, 1)
    denominator = np.sqrt(
        x_unbiased / max(x_n, 1) + y_unbiased / max(y_n, 1)
    )
    return (np.asarray(x_mean, float) - np.asarray(y_mean, float)) / np.maximum(
        denominator, EPS
    )


def pearson(left: np.ndarray, right: np.ndarray) -> float:
    """Stable Pearson correlation used by endpoint evidence gates."""
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    left = left - left.mean()
    right = right - right.mean()
    return float(
        np.sum(left * right)
        / max(np.linalg.norm(left) * np.linalg.norm(right), EPS)
    )


def one_sided_lower(values: np.ndarray, confidence: float = 0.95) -> float:
    """One-sided Student-t lower confidence bound from the frozen contract."""
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return float("-inf")
    standard_error = float(values.std(ddof=1) / np.sqrt(len(values)))
    critical = float(stats.t.ppf(confidence, len(values) - 1))
    return float(values.mean() - critical * standard_error)


def _self_statistic(
    effects: Mapping[str, np.ndarray], gene_index: Mapping[str, int]
) -> float:
    values = [
        float(effect[gene_index[node]])
        for node, effect in effects.items()
        if node in gene_index
    ]
    return float(np.mean(values)) if values else 0.0


def _jaccard(left: set[Any], right: set[Any]) -> float:
    return len(left & right) / max(len(left | right), 1)


def _go_program(
    node: str,
    effects: Mapping[str, np.ndarray],
    gene2go: Mapping[str, Sequence[Any]],
    top_k: int,
) -> np.ndarray:
    background = np.mean(np.stack(list(effects.values())), axis=0)
    target_go = set(gene2go.get(node, ()))
    neighbors: list[tuple[float, str, np.ndarray]] = []
    for known, effect in effects.items():
        similarity = _jaccard(target_go, set(gene2go.get(known, ())))
        neighbors.append((similarity, known, effect))
    neighbors.sort(key=lambda row: (row[0], row[1]), reverse=True)
    selected = neighbors[:top_k]
    weights = np.asarray([row[0] for row in selected], dtype=np.float64)
    if weights.sum() <= EPS:
        return background
    return np.average(
        np.stack([row[2] for row in selected]), axis=0, weights=weights
    )


@dataclass(frozen=True)
class _EndpointFeatures:
    background: np.ndarray
    go_residual: np.ndarray
    self_anchor: np.ndarray


def _endpoint_features(
    node: str,
    effects: Mapping[str, np.ndarray],
    gene_index: Mapping[str, int],
    gene2go: Mapping[str, Sequence[Any]],
    go_top_k: int,
) -> _EndpointFeatures:
    background = np.mean(np.stack(list(effects.values())), axis=0)
    go_residual = _go_program(node, effects, gene2go, go_top_k) - background
    self_anchor = np.zeros_like(background)
    if node in gene_index:
        self_anchor[gene_index[node]] = _self_statistic(effects, gene_index)
    return _EndpointFeatures(background, go_residual, self_anchor)


def _solve_nonnegative_ols(
    feature_rows: Sequence[Sequence[np.ndarray]],
    targets: Sequence[np.ndarray],
    upper: float = 2.0,
) -> np.ndarray:
    width = len(feature_rows[0])
    gram = np.zeros((width, width), dtype=np.float64)
    cross = np.zeros(width, dtype=np.float64)
    for row, target in zip(feature_rows, targets, strict=True):
        design = np.stack(row, axis=1).astype(np.float64)
        response = np.asarray(target, dtype=np.float64)
        gram += design.T @ design
        cross += design.T @ response
    weight = np.linalg.solve(gram + 1e-8 * np.eye(width), cross)
    return np.clip(weight, 0.0, upper)


@dataclass(frozen=True)
class DenseSparseEndpointHead:
    """Frozen v13 dense-background and sparse-self endpoint head."""

    dense_active: bool
    sparse_active: bool
    dense_weights: tuple[float, float]
    joint_weights: tuple[float, float, float]
    dense_all_gene_gain_mean: float
    dense_all_gene_gain_lower: float
    final_all_gene_gain_mean: float
    final_all_gene_gain_lower: float
    final_top100_pcc_delta_mean: float
    final_top100_pcc_delta_lower: float
    known_single_count: int
    records: tuple[dict[str, float | str], ...]
    go_top_k: int = 10

    @property
    def active(self) -> bool:
        return self.dense_active

    def predict(
        self,
        node: str,
        effects: Mapping[str, np.ndarray],
        gene_index: Mapping[str, int],
        gene2go: Mapping[str, Sequence[Any]],
    ) -> np.ndarray:
        row = _endpoint_features(
            node, effects, gene_index, gene2go, self.go_top_k
        )
        if not self.dense_active:
            return row.background
        if self.sparse_active:
            weight = np.asarray(self.joint_weights)
            return (
                weight[0] * row.background
                + weight[1] * row.go_residual
                + weight[2] * row.self_anchor
            )
        weight = np.asarray(self.dense_weights)
        return weight[0] * row.background + weight[1] * row.go_residual


def fit_dense_sparse_head(
    single_effects: Mapping[str, np.ndarray],
    single_means: Mapping[str, np.ndarray],
    single_variances: Mapping[str, np.ndarray],
    single_counts: Mapping[str, int],
    control_mean: np.ndarray,
    control_variance: np.ndarray,
    control_count: int,
    genes: Sequence[str],
    gene2go: Mapping[str, Sequence[Any]],
    gate_confidence: float = 0.95,
) -> DenseSparseEndpointHead:
    """Fit the exact Gate 18/v13 endpoint program used by WitnessCell v14."""
    gene_index = {str(gene): index for index, gene in enumerate(genes)}
    nodes = sorted(node for node in single_effects if node in gene_index)
    if len(nodes) < 4:
        return DenseSparseEndpointHead(
            False,
            False,
            (1.0, 0.0),
            (1.0, 0.0, 0.0),
            0.0,
            float("-inf"),
            0.0,
            float("-inf"),
            0.0,
            float("-inf"),
            len(nodes),
            (),
        )

    rows: list[_EndpointFeatures] = []
    targets: list[np.ndarray] = []
    top_indices: list[np.ndarray] = []
    for held in nodes:
        remaining = {node: single_effects[node] for node in nodes if node != held}
        rows.append(_endpoint_features(held, remaining, gene_index, gene2go, 10))
        targets.append(np.asarray(single_effects[held], dtype=float))
        score = welch_scores(
            single_means[held],
            single_variances[held],
            single_counts[held],
            control_mean,
            control_variance,
            control_count,
        )
        top_indices.append(
            np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]
        )

    dense_weight = _solve_nonnegative_ols(
        [(row.background, row.go_residual) for row in rows], targets
    )
    joint_weight = _solve_nonnegative_ols(
        [(row.background, row.go_residual, row.self_anchor) for row in rows],
        targets,
    )

    records: list[dict[str, float | str]] = []
    for held, row, truth, top in zip(
        nodes, rows, targets, top_indices, strict=True
    ):
        baseline = row.background
        dense = dense_weight[0] * row.background + dense_weight[1] * row.go_residual
        final = (
            joint_weight[0] * row.background
            + joint_weight[1] * row.go_residual
            + joint_weight[2] * row.self_anchor
        )
        base_all = float(np.mean(np.square(baseline - truth)))
        dense_all = float(np.mean(np.square(dense - truth)))
        final_all = float(np.mean(np.square(final - truth)))
        base_top = float(np.mean(np.square(baseline[top] - truth[top])))
        dense_top = float(np.mean(np.square(dense[top] - truth[top])))
        final_top = float(np.mean(np.square(final[top] - truth[top])))
        records.append(
            {
                "condition": held,
                "baseline_all_gene_mse": base_all,
                "dense_all_gene_mse": dense_all,
                "final_all_gene_mse": final_all,
                "dense_all_gene_gain": (base_all - dense_all) / max(base_all, EPS),
                "final_all_gene_gain": (base_all - final_all) / max(base_all, EPS),
                "baseline_top100_mse": base_top,
                "dense_top100_mse": dense_top,
                "final_top100_mse": final_top,
                "final_top100_gain": (base_top - final_top) / max(base_top, EPS),
                "baseline_top100_pcc": pearson(baseline[top], truth[top]),
                "dense_top100_pcc": pearson(dense[top], truth[top]),
                "final_top100_pcc": pearson(final[top], truth[top]),
            }
        )

    def values(key: str) -> np.ndarray:
        return np.asarray([float(record[key]) for record in records])

    dense_all_gain = (
        values("baseline_all_gene_mse") - values("dense_all_gene_mse")
    ) / np.maximum(values("baseline_all_gene_mse"), EPS)
    final_all_gain = (
        values("baseline_all_gene_mse") - values("final_all_gene_mse")
    ) / np.maximum(values("baseline_all_gene_mse"), EPS)
    final_pcc_delta = values("final_top100_pcc") - values(
        "baseline_top100_pcc"
    )
    dense_lower = one_sided_lower(dense_all_gain, gate_confidence)
    final_all_lower = one_sided_lower(final_all_gain, gate_confidence)
    final_pcc_lower = one_sided_lower(final_pcc_delta, gate_confidence)
    sparse_active = bool(final_all_lower > 0.0 and final_pcc_lower > 0.0)
    dense_active = bool(dense_lower > 0.0 or sparse_active)
    return DenseSparseEndpointHead(
        dense_active=dense_active,
        sparse_active=sparse_active,
        dense_weights=(float(dense_weight[0]), float(dense_weight[1])),
        joint_weights=(
            float(joint_weight[0]),
            float(joint_weight[1]),
            float(joint_weight[2]),
        ),
        dense_all_gene_gain_mean=float(dense_all_gain.mean()),
        dense_all_gene_gain_lower=float(dense_lower),
        final_all_gene_gain_mean=float(final_all_gain.mean()),
        final_all_gene_gain_lower=float(final_all_lower),
        final_top100_pcc_delta_mean=float(final_pcc_delta.mean()),
        final_top100_pcc_delta_lower=float(final_pcc_lower),
        known_single_count=len(nodes),
        records=tuple(records),
    )


def _amplitude_features(
    node: str,
    effects: Mapping[str, np.ndarray],
    gene_index: Mapping[str, int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    width = len(next(iter(effects.values())))
    self_anchor = np.zeros(width, dtype=np.float64)
    tail_anchor = np.zeros(width, dtype=np.float64)
    rms_anchor = np.zeros(width, dtype=np.float64)
    target_index = gene_index.get(node)
    if target_index is None:
        return self_anchor, tail_anchor, rms_anchor
    matrix = np.stack(list(effects.values())).astype(np.float64)
    self_anchor[target_index] = _self_statistic(effects, gene_index)
    fingerprint = matrix[:, target_index]
    tail_quantile = 0.10 if self_anchor[target_index] < 0.0 else 0.90
    tail_anchor[target_index] = float(np.quantile(fingerprint, tail_quantile))
    rms_anchor[target_index] = float(
        np.sign(self_anchor[target_index]) * np.sqrt(np.mean(np.square(fingerprint)))
    )
    return self_anchor, tail_anchor, rms_anchor


def _solve_amplitude(
    feature_rows: Sequence[tuple[np.ndarray, np.ndarray, np.ndarray]],
    targets: Sequence[np.ndarray],
    indices: Sequence[np.ndarray],
) -> np.ndarray:
    width = 3
    gram = np.zeros((width, width), dtype=np.float64)
    cross = np.zeros(width, dtype=np.float64)
    for features, target, keep in zip(
        feature_rows, targets, indices, strict=True
    ):
        design = np.stack(features, axis=1).astype(np.float64)[keep]
        response = np.asarray(target, dtype=np.float64)[keep]
        gram += design.T @ design
        cross += design.T @ response
    penalty = 1e-4 * np.eye(width)
    penalty += 0.2 * np.diag(np.maximum(np.diag(gram), EPS))
    weight = np.linalg.solve(gram + penalty, cross)
    return np.clip(weight, -1.5, 1.5)


@dataclass(frozen=True)
class IncrementalAmplitudeHead:
    """Frozen Gate 19 one-coordinate response-fingerprint correction."""

    baseline_head: DenseSparseEndpointHead
    correction_active: bool
    correction_weights: tuple[float, float, float]
    all_gene_upgrade_mean: float
    all_gene_upgrade_lower: float
    top100_mse_upgrade_mean: float
    top100_mse_upgrade_lower: float
    top100_pcc_upgrade_mean: float
    top100_pcc_upgrade_lower: float
    known_single_count: int
    records: tuple[dict[str, float | str], ...]
    ridge: float = 0.2

    def predict(
        self,
        node: str,
        effects: Mapping[str, np.ndarray],
        gene_index: Mapping[str, int],
        gene2go: Mapping[str, Sequence[Any]],
    ) -> np.ndarray:
        baseline = self.baseline_head.predict(node, effects, gene_index, gene2go)
        if not self.correction_active:
            return baseline
        features = _amplitude_features(node, effects, gene_index)
        correction = sum(
            weight * feature
            for weight, feature in zip(
                self.correction_weights, features, strict=True
            )
        )
        return baseline + correction


def _fit_amplitude_weights(
    nodes: Sequence[str],
    effects: Mapping[str, np.ndarray],
    baseline: DenseSparseEndpointHead,
    gene_index: Mapping[str, int],
    gene2go: Mapping[str, Sequence[Any]],
    top_indices: Mapping[str, np.ndarray],
) -> np.ndarray:
    feature_rows: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    residuals: list[np.ndarray] = []
    tops: list[np.ndarray] = []
    for held in nodes:
        remaining = {node: effects[node] for node in nodes if node != held}
        feature_rows.append(_amplitude_features(held, remaining, gene_index))
        residuals.append(
            np.asarray(effects[held], float)
            - baseline.predict(held, remaining, gene_index, gene2go)
        )
        tops.append(top_indices[held])
    return _solve_amplitude(feature_rows, residuals, tops)


def fit_incremental_amplitude_head(
    single_effects: Mapping[str, np.ndarray],
    single_means: Mapping[str, np.ndarray],
    single_variances: Mapping[str, np.ndarray],
    single_counts: Mapping[str, int],
    control_mean: np.ndarray,
    control_variance: np.ndarray,
    control_count: int,
    genes: Sequence[str],
    gene2go: Mapping[str, Sequence[Any]],
    gate_confidence: float = 0.95,
) -> IncrementalAmplitudeHead:
    """Fit the complete frozen v14 endpoint head without target outcomes."""
    gene_index = {str(gene): index for index, gene in enumerate(genes)}
    nodes = sorted(node for node in single_effects if node in gene_index)
    baseline = fit_dense_sparse_head(
        single_effects,
        single_means,
        single_variances,
        single_counts,
        control_mean,
        control_variance,
        control_count,
        genes,
        gene2go,
        gate_confidence,
    )
    if len(nodes) < 6:
        return IncrementalAmplitudeHead(
            baseline,
            False,
            (0.0, 0.0, 0.0),
            0.0,
            float("-inf"),
            0.0,
            float("-inf"),
            0.0,
            float("-inf"),
            len(nodes),
            (),
        )

    top_indices = {
        node: np.argsort(
            -np.abs(
                welch_scores(
                    single_means[node],
                    single_variances[node],
                    single_counts[node],
                    control_mean,
                    control_variance,
                    control_count,
                )
            ),
            kind="stable",
        )[: min(100, len(genes))]
        for node in nodes
    }

    records: list[dict[str, float | str]] = []
    for outer in nodes:
        remaining_nodes = [node for node in nodes if node != outer]
        remaining_effects = {
            node: single_effects[node] for node in remaining_nodes
        }
        outer_baseline = fit_dense_sparse_head(
            remaining_effects,
            single_means,
            single_variances,
            single_counts,
            control_mean,
            control_variance,
            control_count,
            genes,
            gene2go,
            gate_confidence,
        )
        correction_weight = _fit_amplitude_weights(
            remaining_nodes,
            remaining_effects,
            outer_baseline,
            gene_index,
            gene2go,
            top_indices,
        )
        base = outer_baseline.predict(
            outer, remaining_effects, gene_index, gene2go
        )
        correction = sum(
            weight * feature
            for weight, feature in zip(
                correction_weight,
                _amplitude_features(outer, remaining_effects, gene_index),
                strict=True,
            )
        )
        final = base + correction
        truth = np.asarray(single_effects[outer], float)
        top = top_indices[outer]
        base_all = float(np.mean(np.square(base - truth)))
        final_all = float(np.mean(np.square(final - truth)))
        base_top = float(np.mean(np.square(base[top] - truth[top])))
        final_top = float(np.mean(np.square(final[top] - truth[top])))
        base_pcc = pearson(base[top], truth[top])
        final_pcc = pearson(final[top], truth[top])
        records.append(
            {
                "condition": outer,
                "v13_all_gene_mse": base_all,
                "v14_all_gene_mse": final_all,
                "upgrade_all_gene_gain": (base_all - final_all)
                / max(base_all, EPS),
                "v13_top100_mse": base_top,
                "v14_top100_mse": final_top,
                "upgrade_top100_mse_gain": (base_top - final_top)
                / max(base_top, EPS),
                "v13_top100_pcc": base_pcc,
                "v14_top100_pcc": final_pcc,
                "upgrade_top100_pcc_delta": final_pcc - base_pcc,
            }
        )

    all_gain = np.asarray([float(row["upgrade_all_gene_gain"]) for row in records])
    top_gain = np.asarray(
        [float(row["upgrade_top100_mse_gain"]) for row in records]
    )
    pcc_delta = np.asarray(
        [float(row["upgrade_top100_pcc_delta"]) for row in records]
    )
    all_lower = one_sided_lower(all_gain, gate_confidence)
    top_lower = one_sided_lower(top_gain, gate_confidence)
    pcc_lower = one_sided_lower(pcc_delta, gate_confidence)
    active = bool(all_lower > 0.0 and pcc_lower > 0.0)
    final_weight = _fit_amplitude_weights(
        nodes,
        single_effects,
        baseline,
        gene_index,
        gene2go,
        top_indices,
    )
    return IncrementalAmplitudeHead(
        baseline_head=baseline,
        correction_active=active,
        correction_weights=(
            float(final_weight[0]),
            float(final_weight[1]),
            float(final_weight[2]),
        ),
        all_gene_upgrade_mean=float(all_gain.mean()),
        all_gene_upgrade_lower=float(all_lower),
        top100_mse_upgrade_mean=float(top_gain.mean()),
        top100_mse_upgrade_lower=float(top_lower),
        top100_pcc_upgrade_mean=float(pcc_delta.mean()),
        top100_pcc_upgrade_lower=float(pcc_lower),
        known_single_count=len(nodes),
        records=tuple(records),
    )


@dataclass(frozen=True)
class FittedWitnessCore:
    """Immutable fitted state for the frozen condition-mean predictor."""

    genes: tuple[str, ...]
    control_label: str
    control_mean: np.ndarray
    known_single_effects: Mapping[str, np.ndarray]
    gene2go: Mapping[str, tuple[str, ...]]
    endpoint_head: IncrementalAmplitudeHead
    alpha: float
    noise_ratio: float
    gamma: float
    train_doubles: tuple[str, ...]
    train_residual: np.ndarray
    validation_doubles: tuple[str, ...]
    validation_mse_factorized: float
    validation_mse_witness: float

    @property
    def gene_index(self) -> dict[str, int]:
        return {gene: index for index, gene in enumerate(self.genes)}

    def endpoint_effect(self, node: str, corrected: bool = True) -> np.ndarray:
        if node in self.known_single_effects:
            return np.asarray(self.known_single_effects[node], dtype=float)
        if not self.known_single_effects:
            return np.zeros_like(self.control_mean)
        if corrected:
            return self.endpoint_head.predict(
                node,
                self.known_single_effects,
                self.gene_index,
                self.gene2go,
            )
        return self.endpoint_head.baseline_head.predict(
            node,
            self.known_single_effects,
            self.gene_index,
            self.gene2go,
        )

    def factorized_effect(self, condition: str, corrected: bool = True) -> np.ndarray:
        if condition == self.control_label:
            return np.zeros_like(self.control_mean)
        parts = parse_condition(condition)
        if len(parts) == 1:
            return self.endpoint_effect(parts[0], corrected=corrected)
        additive = np.sum(
            [self.endpoint_effect(node, corrected=corrected) for node in parts],
            axis=0,
        )
        return saturate(additive, self.alpha)

    def predict_effects(
        self, conditions: Sequence[str]
    ) -> tuple[np.ndarray, np.ndarray]:
        conditions = tuple(conditions)
        factorized = np.stack(
            [self.factorized_effect(condition, corrected=True) for condition in conditions]
        )
        corrected = factorized.copy()
        target_rows = [
            index
            for index, condition in enumerate(conditions)
            if condition != self.control_label and len(parse_condition(condition)) == 2
        ]
        if self.train_doubles and target_rows and self.gamma != 0.0:
            target_doubles = [conditions[index] for index in target_rows]
            nodes = sorted(
                {
                    node
                    for condition in self.train_doubles + tuple(target_doubles)
                    for node in parse_condition(condition)
                }
            )
            raw = kernel_predict(
                incidence(self.train_doubles, nodes),
                incidence(target_doubles, nodes),
                self.train_residual,
                self.noise_ratio,
            )
            for local, output_index in enumerate(target_rows):
                corrected[output_index] += self.gamma * raw[local]
        return corrected, factorized


def fit_witness_core(
    means: Mapping[str, np.ndarray],
    variances: Mapping[str, np.ndarray],
    counts: Mapping[str, int],
    genes: Sequence[str],
    train_conditions: Sequence[str],
    validation_conditions: Sequence[str],
    gene2go: Mapping[str, Sequence[Any]],
    control_label: str = "control",
    alpha_grid: Sequence[float] = DEFAULT_ALPHA_GRID,
    noise_grid: Sequence[float] = DEFAULT_NOISE_GRID,
    gate_confidence: float = 0.95,
) -> FittedWitnessCore:
    """Fit WitnessCell v14 using training singles and validation doubles only."""
    control = np.asarray(means[control_label], dtype=float)
    train_conditions = tuple(train_conditions)
    validation_conditions = tuple(validation_conditions)
    known_single = {
        condition: np.asarray(means[condition], dtype=float) - control
        for condition in train_conditions
        if condition != control_label and len(parse_condition(condition)) == 1
    }
    normalized_go = {
        str(gene): tuple(sorted(map(str, terms))) for gene, terms in gene2go.items()
    }
    endpoint_head = fit_incremental_amplitude_head(
        single_effects=known_single,
        single_means={condition: np.asarray(means[condition], float) for condition in known_single},
        single_variances={condition: np.asarray(variances[condition], float) for condition in known_single},
        single_counts={condition: int(counts[condition]) for condition in known_single},
        control_mean=control,
        control_variance=np.asarray(variances[control_label], float),
        control_count=int(counts[control_label]),
        genes=tuple(map(str, genes)),
        gene2go=normalized_go,
        gate_confidence=gate_confidence,
    )

    gene_index = {str(gene): index for index, gene in enumerate(genes)}

    def endpoint_effect(node: str, corrected: bool) -> np.ndarray:
        if node in known_single:
            return known_single[node]
        if not known_single:
            return np.zeros_like(control)
        if corrected:
            return endpoint_head.predict(node, known_single, gene_index, normalized_go)
        return endpoint_head.baseline_head.predict(
            node, known_single, gene_index, normalized_go
        )

    def additive(condition: str, corrected: bool) -> np.ndarray:
        return np.sum(
            [endpoint_effect(node, corrected) for node in parse_condition(condition)],
            axis=0,
        )

    train_doubles = tuple(
        condition
        for condition in train_conditions
        if condition != control_label and len(parse_condition(condition)) == 2
    )
    validation_doubles = tuple(
        condition
        for condition in validation_conditions
        if condition != control_label and len(parse_condition(condition)) == 2
    )

    if not train_doubles:
        if validation_doubles:
            factorized_candidates: list[tuple[float, float]] = []
            truth = np.stack([means[condition] - control for condition in validation_doubles])
            for alpha_value in alpha_grid:
                base = np.stack(
                    [saturate(additive(condition, False), alpha_value) for condition in validation_doubles]
                )
                factorized_candidates.append(
                    (float(np.mean((base - truth) ** 2)), float(alpha_value))
                )
            factorized_mse, alpha = min(factorized_candidates)
        else:
            alpha = 0.0
            factorized_mse = float("nan")
        return FittedWitnessCore(
            genes=tuple(map(str, genes)),
            control_label=control_label,
            control_mean=control,
            known_single_effects=known_single,
            gene2go=normalized_go,
            endpoint_head=endpoint_head,
            alpha=float(alpha),
            noise_ratio=0.0,
            gamma=0.0,
            train_doubles=(),
            train_residual=np.zeros((0, len(control))),
            validation_doubles=validation_doubles,
            validation_mse_factorized=float(factorized_mse),
            validation_mse_witness=float(factorized_mse),
        )

    nodes = sorted(
        {
            node
            for condition in train_doubles + validation_doubles
            for node in parse_condition(condition)
        }
    )
    train_design = incidence(train_doubles, nodes)
    validation_design = (
        incidence(validation_doubles, nodes)
        if validation_doubles
        else np.zeros((0, len(nodes)))
    )
    candidates: list[tuple[float, float, float, float, float]] = []
    for alpha_value in alpha_grid:
        train_base = np.stack(
            [saturate(additive(condition, False), alpha_value) for condition in train_doubles]
        )
        train_truth = np.stack([means[condition] - control for condition in train_doubles])
        residual = train_truth - train_base
        if validation_doubles:
            val_base = np.stack(
                [saturate(additive(condition, False), alpha_value) for condition in validation_doubles]
            )
            val_truth = np.stack([means[condition] - control for condition in validation_doubles])
            for noise in noise_grid:
                raw = kernel_predict(train_design, validation_design, residual, noise)
                gamma = closed_form_gamma(raw, val_truth - val_base)
                witness_mse = float(np.mean((val_base + gamma * raw - val_truth) ** 2))
                factorized_mse = float(np.mean((val_base - val_truth) ** 2))
                candidates.append(
                    (witness_mse, float(alpha_value), float(noise), gamma, factorized_mse)
                )
        else:
            candidates.append((float("nan"), float(alpha_value), 0.1, 1.0, float("nan")))

    if validation_doubles:
        candidates.sort(key=lambda row: (row[0], row[1], row[2]))
    witness_mse, alpha, noise, gamma, factorized_mse = candidates[0]
    train_base = np.stack(
        [saturate(additive(condition, False), alpha) for condition in train_doubles]
    )
    train_truth = np.stack([means[condition] - control for condition in train_doubles])
    return FittedWitnessCore(
        genes=tuple(map(str, genes)),
        control_label=control_label,
        control_mean=control,
        known_single_effects=known_single,
        gene2go=normalized_go,
        endpoint_head=endpoint_head,
        alpha=float(alpha),
        noise_ratio=float(noise),
        gamma=float(gamma),
        train_doubles=train_doubles,
        train_residual=train_truth - train_base,
        validation_doubles=validation_doubles,
        validation_mse_factorized=float(factorized_mse),
        validation_mse_witness=float(witness_mse),
    )
