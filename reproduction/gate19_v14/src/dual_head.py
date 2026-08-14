#!/usr/bin/env python3
"""Training-only dense background plus sparse directional endpoint head.

The head is deliberately split into two sequential decisions:

1. the final endpoint program may replace the frozen v1 mean fallback only
   when cross-fitted all-gene LOO evidence has a positive one-sided 95% LCB;
2. a sparse direction correction may be added only when it has positive
   cross-fitted top-100 LOO evidence both versus the v1 fallback and versus
   the dense head.  If that directional gate rejects, the dense background
   must pass the all-gene gate on its own or the model falls back to v1.

All masks, coefficients and gates use known training singles only.  Held-out
test means and validation combinations never enter this module.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping, Sequence

import numpy as np

from identity_head import EPS, one_sided_lower, pearson, self_statistic, welch_scores


def jaccard(left: set, right: set) -> float:
    return len(left & right) / max(len(left | right), 1)


def go_program(
    node: str,
    effects: Mapping[str, np.ndarray],
    gene2go: Mapping,
    top_k: int,
) -> np.ndarray:
    background = np.mean(np.stack(list(effects.values())), axis=0)
    target_go = set(gene2go.get(node, ()))
    neighbors = []
    for known, effect in effects.items():
        similarity = jaccard(target_go, set(gene2go.get(known, ())))
        neighbors.append((similarity, known, effect))
    neighbors.sort(key=lambda row: (row[0], row[1]), reverse=True)
    selected = neighbors if top_k < 0 else neighbors[:top_k]
    weights = np.asarray([row[0] for row in selected], dtype=np.float64)
    if weights.sum() <= EPS:
        return background
    return np.average(
        np.stack([row[2] for row in selected]), axis=0, weights=weights
    )


def top_support(score: np.ndarray, size: int) -> np.ndarray:
    size = min(max(int(size), 1), len(score))
    order = np.argsort(-np.abs(score), kind="stable")[:size]
    mask = np.zeros(len(score), dtype=np.float64)
    mask[order] = 1.0
    return mask


def support_score(
    mode: str,
    background: np.ndarray,
    neighbor_program: np.ndarray,
    effects: Mapping[str, np.ndarray],
) -> np.ndarray:
    go_residual = neighbor_program - background
    if mode == "go_residual":
        return np.abs(go_residual)
    if mode == "go_program":
        return np.abs(neighbor_program)
    activity = np.sqrt(np.mean(np.square(np.stack(list(effects.values()))), axis=0))
    if mode == "global_activity":
        return activity
    if mode == "hybrid":
        # Endpoint-specific GO disagreement, softly weighted by whether a gene
        # is responsive at all in the known training-single environment.
        return np.abs(go_residual) * np.sqrt(activity + EPS)
    raise ValueError(f"unknown support mode: {mode}")


def direction_vote(
    node: str,
    effects: Mapping[str, np.ndarray],
    direction_masks: Mapping[str, np.ndarray],
    gene2go: Mapping,
    top_k: int,
    weighted: bool,
) -> np.ndarray:
    target_go = set(gene2go.get(node, ()))
    neighbors = []
    for known in effects:
        similarity = jaccard(target_go, set(gene2go.get(known, ())))
        neighbors.append((similarity, known))
    neighbors.sort(key=lambda row: (row[0], row[1]), reverse=True)
    selected = neighbors if top_k < 0 else neighbors[:top_k]
    if not selected:
        return np.zeros_like(next(iter(effects.values())))
    if weighted:
        weights = np.asarray([row[0] for row in selected], dtype=np.float64)
        if weights.sum() <= EPS:
            weights = np.ones(len(selected), dtype=np.float64)
    else:
        weights = np.ones(len(selected), dtype=np.float64)
    return np.average(
        np.stack([direction_masks[row[1]] for row in selected]),
        axis=0,
        weights=weights,
    )


@dataclass(frozen=True)
class FeatureRow:
    background: np.ndarray
    go_residual: np.ndarray
    sparse_go: np.ndarray
    self_anchor: np.ndarray
    support_mask: np.ndarray


def features_for_node(
    node: str,
    effects: Mapping[str, np.ndarray],
    gene_index: Mapping[str, int],
    gene2go: Mapping,
    go_top_k: int,
    support_mode: str,
    support_k: int,
    direction_masks: Mapping[str, np.ndarray] | None = None,
) -> FeatureRow:
    background = np.mean(np.stack(list(effects.values())), axis=0)
    neighbor = go_program(node, effects, gene2go, go_top_k)
    residual = neighbor - background
    if support_mode in ("go_de_vote", "global_de_frequency", "de_hybrid"):
        if direction_masks is None:
            raise ValueError(f"{support_mode} requires training direction masks")
        if support_mode == "global_de_frequency":
            score = np.mean(np.stack(list(direction_masks.values())), axis=0)
        else:
            vote = direction_vote(
                node,
                effects,
                direction_masks,
                gene2go,
                go_top_k,
                weighted=True,
            )
            score = vote if support_mode == "go_de_vote" else vote * (np.abs(residual) + EPS)
    else:
        score = support_score(support_mode, background, neighbor, effects)
    mask = top_support(score, support_k)
    anchor = np.zeros_like(background)
    if node in gene_index:
        anchor[gene_index[node]] = self_statistic(effects, gene_index, "mean")
        mask[gene_index[node]] = 1.0
    return FeatureRow(
        background=background,
        go_residual=residual,
        sparse_go=residual * mask,
        self_anchor=anchor,
        support_mask=mask,
    )


def solve_ols(
    features: Sequence[Sequence[np.ndarray]],
    targets: Sequence[np.ndarray],
    indices: Sequence[np.ndarray] | None = None,
    upper: float = 2.0,
) -> np.ndarray:
    width = len(features[0])
    gram = np.zeros((width, width), dtype=np.float64)
    cross = np.zeros(width, dtype=np.float64)
    for row_index, (row, target) in enumerate(zip(features, targets, strict=True)):
        design = np.stack(row, axis=1).astype(np.float64)
        y = np.asarray(target, dtype=np.float64)
        if indices is not None:
            keep = indices[row_index]
            design = design[keep]
            y = y[keep]
        gram += design.T @ design
        cross += design.T @ y
    weight = np.linalg.solve(gram + 1e-8 * np.eye(width), cross)
    return np.clip(weight, 0.0, upper)


def solve_signed_ols(
    features: Sequence[Sequence[np.ndarray]],
    targets: Sequence[np.ndarray],
    indices: Sequence[np.ndarray] | None,
    lower: float,
    upper: float,
) -> np.ndarray:
    width = len(features[0])
    gram = np.zeros((width, width), dtype=np.float64)
    cross = np.zeros(width, dtype=np.float64)
    for row_index, (row, target) in enumerate(zip(features, targets, strict=True)):
        design = np.stack(row, axis=1).astype(np.float64)
        y = np.asarray(target, dtype=np.float64)
        if indices is not None:
            keep = indices[row_index]
            design = design[keep]
            y = y[keep]
        gram += design.T @ design
        cross += design.T @ y
    weight = np.linalg.solve(gram + 1e-8 * np.eye(width), cross)
    return np.clip(weight, lower, upper)


def fit_dense_weights(
    rows: Sequence[FeatureRow],
    targets: Sequence[np.ndarray],
    dense_mode: str,
) -> np.ndarray:
    if dense_mode == "mean_only":
        return solve_ols(
            [(row.background,) for row in rows], targets
        )
    if dense_mode == "mean_go":
        return solve_ols(
            [(row.background, row.go_residual) for row in rows], targets
        )
    if dense_mode == "mean_self":
        return solve_ols(
            [(row.background, row.self_anchor) for row in rows], targets
        )
    if dense_mode == "mean_go_self":
        return solve_ols(
            [(row.background, row.go_residual, row.self_anchor) for row in rows],
            targets,
        )
    raise ValueError(f"unknown dense mode: {dense_mode}")


def dense_prediction(row: FeatureRow, weight: np.ndarray, dense_mode: str) -> np.ndarray:
    if dense_mode == "mean_only":
        return weight[0] * row.background
    if dense_mode == "mean_go":
        return weight[0] * row.background + weight[1] * row.go_residual
    if dense_mode == "mean_self":
        return weight[0] * row.background + weight[1] * row.self_anchor
    if dense_mode == "mean_go_self":
        return (
            weight[0] * row.background
            + weight[1] * row.go_residual
            + weight[2] * row.self_anchor
        )
    raise ValueError(f"unknown dense mode: {dense_mode}")


def fit_sparse_weights(
    rows: Sequence[FeatureRow],
    targets: Sequence[np.ndarray],
    dense_weight: np.ndarray,
    dense_mode: str,
    top_indices: Sequence[np.ndarray],
    sparse_mode: str,
) -> np.ndarray:
    dense = [dense_prediction(row, dense_weight, dense_mode) for row in rows]
    residual_targets = [target - pred for target, pred in zip(targets, dense, strict=True)]
    if sparse_mode == "self_only":
        sparse_features = [(row.self_anchor,) for row in rows]
    elif sparse_mode == "go_only":
        sparse_features = [(row.sparse_go,) for row in rows]
    elif sparse_mode == "go_self":
        sparse_features = [(row.sparse_go, row.self_anchor) for row in rows]
    elif sparse_mode == "scale_dense":
        sparse_features = [
            (row.support_mask * prediction,)
            for row, prediction in zip(rows, dense, strict=True)
        ]
        return solve_signed_ols(
            sparse_features,
            residual_targets,
            indices=top_indices,
            lower=-1.0,
            upper=1.0,
        )
    else:
        raise ValueError(f"unknown sparse mode: {sparse_mode}")
    return solve_ols(
        sparse_features,
        residual_targets,
        indices=top_indices,
        upper=2.5,
    )


def combined_prediction(
    row: FeatureRow,
    dense_weight: np.ndarray,
    sparse_weight: np.ndarray,
    dense_mode: str,
    sparse_mode: str,
    combined_dense_weight: np.ndarray | None = None,
) -> np.ndarray:
    value = dense_prediction(
        row,
        dense_weight if combined_dense_weight is None else combined_dense_weight,
        dense_mode,
    )
    if sparse_mode == "joint_self":
        return value + sparse_weight[0] * row.self_anchor
    if sparse_mode == "self_only":
        return value + sparse_weight[0] * row.self_anchor
    if sparse_mode == "go_only":
        return value + sparse_weight[0] * row.sparse_go
    if sparse_mode == "scale_dense":
        return value + sparse_weight[0] * row.support_mask * value
    return value + sparse_weight[0] * row.sparse_go + sparse_weight[1] * row.self_anchor


@dataclass(frozen=True)
class DualHeadFit:
    dense_active: bool
    sparse_active: bool
    dense_mode: str
    sparse_mode: str
    support_mode: str
    support_k: int
    go_top_k: int
    dense_weights: tuple[float, ...]
    sparse_weights: tuple[float, ...]
    dense_all_gene_gain_mean: float
    dense_all_gene_gain_lower: float
    final_all_gene_gain_mean: float
    final_all_gene_gain_lower: float
    final_top100_gain_mean: float
    final_top100_gain_lower: float
    sparse_incremental_top100_gain_mean: float
    sparse_incremental_top100_gain_lower: float
    final_top100_pcc_delta_mean: float
    final_top100_pcc_delta_lower: float
    sparse_incremental_top100_pcc_delta_mean: float
    sparse_incremental_top100_pcc_delta_lower: float
    direction_gate_metric: str
    known_single_count: int
    records: tuple[dict, ...]
    training_direction_masks: Mapping[str, np.ndarray]
    selected_recipe: str = "fixed"
    candidate_training_scores: tuple[dict, ...] = ()
    combined_dense_weights: tuple[float, ...] = ()
    incremental_direction_required: bool = False

    @property
    def active(self) -> bool:
        return self.dense_active

    def predict(
        self,
        node: str,
        effects: Mapping[str, np.ndarray],
        gene_index: Mapping[str, int],
        gene2go: Mapping,
    ) -> np.ndarray:
        row = features_for_node(
            node,
            effects,
            gene_index,
            gene2go,
            self.go_top_k,
            self.support_mode,
            self.support_k,
            self.training_direction_masks,
        )
        if not self.dense_active:
            return row.background
        dense_weight = np.asarray(self.dense_weights)
        dense = dense_prediction(row, dense_weight, self.dense_mode)
        if not self.sparse_active:
            return dense
        combined_dense_weight = (
            np.asarray(self.combined_dense_weights)
            if self.combined_dense_weights
            else None
        )
        return combined_prediction(
            row,
            dense_weight,
            np.asarray(self.sparse_weights),
            self.dense_mode,
            self.sparse_mode,
            combined_dense_weight,
        )


def fit_dual_head(
    single_effects: Mapping[str, np.ndarray],
    single_means: Mapping[str, np.ndarray],
    single_variances: Mapping[str, np.ndarray],
    single_counts: Mapping[str, int],
    control_mean: np.ndarray,
    control_variance: np.ndarray,
    control_count: int,
    genes: Sequence[str],
    gene2go: Mapping,
    dense_mode: str = "mean_go",
    sparse_mode: str = "go_self",
    support_mode: str = "hybrid",
    support_k: int = 100,
    go_top_k: int = 20,
    gate_confidence: float = 0.95,
    nested_gate: bool = True,
    direction_gate_metric: str = "top100_pcc",
    require_incremental_direction: bool = False,
) -> DualHeadFit:
    if dense_mode == "auto":
        # A deliberately tiny training-only router.  The two recipes are
        # frozen from the background ablation: the robust global mean and the
        # GO-neighbour dense background.  Each candidate runs its own nested
        # LOO gates; routing uses only those pseudo-held-out records.
        candidate_specs = (
            ("mean_only_k40", "mean_only", 40),
            ("mean_go_k20", "mean_go", 20),
        )
        candidates = []
        scores = []
        for recipe_name, recipe_dense_mode, recipe_go_top_k in candidate_specs:
            candidate = fit_dual_head(
                single_effects=single_effects,
                single_means=single_means,
                single_variances=single_variances,
                single_counts=single_counts,
                control_mean=control_mean,
                control_variance=control_variance,
                control_count=control_count,
                genes=genes,
                gene2go=gene2go,
                dense_mode=recipe_dense_mode,
                sparse_mode=sparse_mode,
                support_mode=support_mode,
                support_k=support_k,
                go_top_k=recipe_go_top_k,
                gate_confidence=gate_confidence,
                nested_gate=nested_gate,
                direction_gate_metric=direction_gate_metric,
                require_incremental_direction=require_incremental_direction,
            )
            baseline_all = sum(
                row["baseline_all_gene_mse"] for row in candidate.records
            )
            baseline_top = sum(
                row["baseline_top100_mse"] for row in candidate.records
            )
            if not candidate.dense_active:
                deployed_all = baseline_all
                deployed_top = baseline_top
                deployed_pcc_delta = 0.0
            elif candidate.sparse_active:
                deployed_all = sum(
                    row["final_all_gene_mse"] for row in candidate.records
                )
                deployed_top = sum(
                    row["final_top100_mse"] for row in candidate.records
                )
                deployed_pcc_delta = float(np.mean([
                    row["final_top100_pcc"] - row["baseline_top100_pcc"]
                    for row in candidate.records
                ]))
            else:
                deployed_all = sum(
                    row["dense_all_gene_mse"] for row in candidate.records
                )
                deployed_top = sum(
                    row["dense_top100_mse"] for row in candidate.records
                )
                deployed_pcc_delta = float(np.mean([
                    row["dense_top100_pcc"] - row["baseline_top100_pcc"]
                    for row in candidate.records
                ]))
            score = {
                "recipe": recipe_name,
                "dense_active": candidate.dense_active,
                "sparse_active": candidate.sparse_active,
                "deployed_all_gene_gain": 1.0 - deployed_all / max(baseline_all, EPS),
                "deployed_top100_gain": 1.0 - deployed_top / max(baseline_top, EPS),
                "deployed_top100_pcc_delta": deployed_pcc_delta,
            }
            candidates.append(candidate)
            scores.append(score)
        selected_index = max(
            range(len(scores)),
            key=lambda index: (
                scores[index]["deployed_top100_gain"],
                scores[index]["deployed_all_gene_gain"],
                scores[index]["deployed_top100_pcc_delta"],
                -index,
            ),
        )
        return replace(
            candidates[selected_index],
            selected_recipe=scores[selected_index]["recipe"],
            candidate_training_scores=tuple(scores),
        )

    gene_index = {str(gene): index for index, gene in enumerate(genes)}
    nodes = sorted(node for node in single_effects if node in gene_index)
    if len(nodes) < 4:
        background_width = {
            "mean_only": 1,
            "mean_go": 2,
            "mean_self": 2,
            "mean_go_self": 3,
        }[dense_mode]
        sparse_width = 1 if sparse_mode in (
            "self_only", "joint_self", "go_only", "scale_dense"
        ) else 2
        return DualHeadFit(
            False, False, dense_mode, sparse_mode, support_mode, support_k, go_top_k,
            tuple([1.0] + [0.0] * (background_width - 1)),
            tuple([0.0] * sparse_width),
            0.0, float("-inf"), 0.0, float("-inf"), 0.0, float("-inf"),
            0.0, float("-inf"), 0.0, float("-inf"), 0.0, float("-inf"),
            direction_gate_metric, len(nodes), (), {},
        )

    direction_masks: dict[str, np.ndarray] = {}
    for node in nodes:
        score = welch_scores(
            single_means[node], single_variances[node], single_counts[node],
            control_mean, control_variance, control_count,
        )
        mask = np.zeros(len(score), dtype=np.float64)
        mask[np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]] = 1.0
        direction_masks[node] = mask

    rows: list[FeatureRow] = []
    targets: list[np.ndarray] = []
    top_indices: list[np.ndarray] = []
    for held in nodes:
        remaining = {node: single_effects[node] for node in nodes if node != held}
        remaining_masks = {
            node: direction_masks[node] for node in nodes if node != held
        }
        rows.append(features_for_node(
            held, remaining, gene_index, gene2go, go_top_k, support_mode, support_k,
            remaining_masks,
        ))
        targets.append(single_effects[held])
        score = welch_scores(
            single_means[held], single_variances[held], single_counts[held],
            control_mean, control_variance, control_count,
        )
        top_indices.append(
            np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]
        )

    final_dense_weight = fit_dense_weights(rows, targets, dense_mode)
    if sparse_mode == "joint_self":
        if dense_mode != "mean_go":
            raise ValueError("joint_self requires dense_mode=mean_go")
        joint_weight = fit_dense_weights(rows, targets, "mean_go_self")
        final_combined_dense_weight = joint_weight[:2]
        final_sparse_weight = joint_weight[2:]
    else:
        final_combined_dense_weight = final_dense_weight
        final_sparse_weight = fit_sparse_weights(
            rows, targets, final_dense_weight, dense_mode, top_indices, sparse_mode
        )

    records = []
    for held_index, (held, row, truth, top) in enumerate(
        zip(nodes, rows, targets, top_indices, strict=True)
    ):
        if nested_gate:
            keep = [index for index in range(len(rows)) if index != held_index]
            train_rows = [rows[index] for index in keep]
            train_targets = [targets[index] for index in keep]
            train_top = [top_indices[index] for index in keep]
            dense_weight = fit_dense_weights(train_rows, train_targets, dense_mode)
            if sparse_mode == "joint_self":
                joint_weight = fit_dense_weights(
                    train_rows, train_targets, "mean_go_self"
                )
                combined_dense_weight = joint_weight[:2]
                sparse_weight = joint_weight[2:]
            else:
                combined_dense_weight = dense_weight
                sparse_weight = fit_sparse_weights(
                    train_rows, train_targets, dense_weight, dense_mode,
                    train_top, sparse_mode
                )
        else:
            dense_weight = final_dense_weight
            combined_dense_weight = final_combined_dense_weight
            sparse_weight = final_sparse_weight
        baseline = row.background
        dense = dense_prediction(row, dense_weight, dense_mode)
        combined = combined_prediction(
            row, dense_weight, sparse_weight, dense_mode, sparse_mode,
            combined_dense_weight,
        )
        base_all = float(np.mean(np.square(baseline - truth)))
        dense_all = float(np.mean(np.square(dense - truth)))
        final_all = float(np.mean(np.square(combined - truth)))
        base_top = float(np.mean(np.square(baseline[top] - truth[top])))
        dense_top = float(np.mean(np.square(dense[top] - truth[top])))
        final_top = float(np.mean(np.square(combined[top] - truth[top])))
        records.append({
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
            "sparse_incremental_top100_gain": (dense_top - final_top) / max(dense_top, EPS),
            "baseline_top100_pcc": pearson(baseline[top], truth[top]),
            "dense_top100_pcc": pearson(dense[top], truth[top]),
            "final_top100_pcc": pearson(combined[top], truth[top]),
            "final_top100_pcc_delta": pearson(combined[top], truth[top]) - pearson(baseline[top], truth[top]),
        })

    def values(key: str) -> np.ndarray:
        return np.asarray([record[key] for record in records], dtype=np.float64)

    dense_all = values("dense_all_gene_gain")
    final_all = values("final_all_gene_gain")
    final_top = values("final_top100_gain")
    incremental_top = values("sparse_incremental_top100_gain")
    final_pcc = values("final_top100_pcc") - values("baseline_top100_pcc")
    incremental_pcc = values("final_top100_pcc") - values("dense_top100_pcc")
    dense_all_lower = one_sided_lower(dense_all, gate_confidence)
    final_all_lower = one_sided_lower(final_all, gate_confidence)
    final_top_lower = one_sided_lower(final_top, gate_confidence)
    incremental_top_lower = one_sided_lower(incremental_top, gate_confidence)
    final_pcc_lower = one_sided_lower(final_pcc, gate_confidence)
    incremental_pcc_lower = one_sided_lower(incremental_pcc, gate_confidence)
    dense_standalone_active = bool(dense_all_lower > 0.0)
    if direction_gate_metric == "top100_pcc":
        sparse_active = bool(
            final_all_lower > 0.0
            and final_pcc_lower > 0.0
            and (
                not require_incremental_direction
                or incremental_pcc_lower > 0.0
            )
        )
    elif direction_gate_metric == "top100_mse":
        sparse_active = bool(
            final_all_lower > 0.0
            and final_top_lower > 0.0
            and (
                not require_incremental_direction
                or incremental_top_lower > 0.0
            )
        )
    else:
        raise ValueError(f"unknown direction gate metric: {direction_gate_metric}")
    # Joint activation is intentional: the all-gene gate protects the final
    # expression program, while the top-100 gate owns the sparse correction.
    # When the sparse correction is rejected, only a stand-alone-safe dense
    # head is deployable.
    dense_active = bool(dense_standalone_active or sparse_active)
    return DualHeadFit(
        dense_active=dense_active,
        sparse_active=sparse_active,
        dense_mode=dense_mode,
        sparse_mode=sparse_mode,
        support_mode=support_mode,
        support_k=support_k,
        go_top_k=go_top_k,
        dense_weights=tuple(map(float, final_dense_weight)),
        sparse_weights=tuple(map(float, final_sparse_weight)),
        dense_all_gene_gain_mean=float(dense_all.mean()),
        dense_all_gene_gain_lower=float(dense_all_lower),
        final_all_gene_gain_mean=float(final_all.mean()),
        final_all_gene_gain_lower=float(final_all_lower),
        final_top100_gain_mean=float(final_top.mean()),
        final_top100_gain_lower=float(final_top_lower),
        sparse_incremental_top100_gain_mean=float(incremental_top.mean()),
        sparse_incremental_top100_gain_lower=float(incremental_top_lower),
        final_top100_pcc_delta_mean=float(final_pcc.mean()),
        final_top100_pcc_delta_lower=float(final_pcc_lower),
        sparse_incremental_top100_pcc_delta_mean=float(incremental_pcc.mean()),
        sparse_incremental_top100_pcc_delta_lower=float(incremental_pcc_lower),
        direction_gate_metric=direction_gate_metric,
        known_single_count=len(nodes),
        records=tuple(records),
        training_direction_masks=direction_masks,
        combined_dense_weights=tuple(map(float, final_combined_dense_weight)),
        incremental_direction_required=require_incremental_direction,
    )
