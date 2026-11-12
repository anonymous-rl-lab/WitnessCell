#!/usr/bin/env python3
"""WitnessCell v2: frozen v1 interaction kernel plus endpoint identity head."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np

from identity_head import IdentityHeadFit, fit_identity_head
from go_identity_head import GOIdentityHeadFit, fit_go_identity_head


EPS = 1e-12


def saturate(effect: np.ndarray, strength: float) -> np.ndarray:
    return effect / (1.0 + strength * np.abs(effect))


def endpoints(condition: str) -> tuple[str, ...]:
    return tuple(condition.split("+")) if "+" in condition else (condition,)


def closed_form_gamma(prediction: np.ndarray, truth: np.ndarray) -> float:
    denominator = float(np.sum(prediction * prediction))
    if denominator <= EPS:
        return 0.0
    return float(np.clip(np.sum(prediction * truth) / denominator, 0.0, 1.0))


def incidence(conditions: list[str], nodes: list[str]) -> np.ndarray:
    node_id = {node: index for index, node in enumerate(nodes)}
    result = np.zeros((len(conditions), len(nodes)), dtype=float)
    for row, condition in enumerate(conditions):
        for node in endpoints(condition):
            result[row, node_id[node]] = 1.0
    return result


def kernel_predict(
    train_design: np.ndarray,
    target_design: np.ndarray,
    response: np.ndarray,
    noise_ratio: float,
) -> np.ndarray:
    train_kernel = train_design @ train_design.T / 2.0
    cross_kernel = train_design @ target_design.T / 2.0
    covariance = train_kernel + noise_ratio * np.eye(len(train_kernel))
    return np.linalg.solve(covariance, cross_kernel).T @ response


@dataclass(frozen=True)
class WitnessFitV2:
    predictions: dict[str, np.ndarray]
    factorized_predictions: dict[str, np.ndarray]
    alpha: float
    noise_ratio: float
    gamma: float
    known_single_genes: tuple[str, ...]
    training_doubles: tuple[str, ...]
    validation_doubles: tuple[str, ...]
    validation_mse_factorized: float
    validation_mse_witness: float
    identity_head: IdentityHeadFit | GOIdentityHeadFit


def fit_predict(
    means: Mapping[str, np.ndarray],
    variances: Mapping[str, np.ndarray],
    counts: Mapping[str, int],
    genes: Sequence[str],
    train_conditions: list[str],
    validation_conditions: list[str],
    target_conditions: list[str],
    alpha_grid: tuple[float, ...] = (0.0, 0.02, 0.05, 0.1, 0.2, 0.5),
    noise_grid: tuple[float, ...] = (0.01, 0.03, 0.1, 0.3, 1.0),
    anchor_mode: str = "add",
    self_summary: str = "mean",
    weight_mode: str = "closed_form",
    gate_metric: str = "all_gene_mse",
    identity_mode: str = "go_shrinkage",
    gene2go: Mapping[str, Any] | None = None,
    go_top_k: int = 40,
) -> WitnessFitV2:
    control = np.asarray(means["control"], dtype=float)
    known_single = {
        condition: np.asarray(means[condition], dtype=float) - control
        for condition in train_conditions
        if condition != "control" and "+" not in condition
    }
    if known_single:
        default_single = np.mean(np.stack(list(known_single.values())), axis=0)
    else:
        default_single = np.zeros_like(control)
    identity_arguments = dict(
        single_effects=known_single,
        single_means={condition: np.asarray(means[condition], dtype=float) for condition in known_single},
        single_variances={condition: np.asarray(variances[condition], dtype=float) for condition in known_single},
        single_counts={condition: int(counts[condition]) for condition in known_single},
        control_mean=control,
        control_variance=np.asarray(variances["control"], dtype=float),
        control_count=int(counts["control"]),
        genes=genes,
    )
    if identity_mode == "go_shrinkage":
        if gene2go is None:
            raise ValueError("gene2go is required for go_shrinkage identity mode")
        identity = fit_go_identity_head(
            **identity_arguments,
            gene2go=gene2go,
            top_k=go_top_k,
        )
    elif identity_mode == "self_only":
        identity = fit_identity_head(
            **identity_arguments,
            anchor_mode=anchor_mode,
            self_summary=self_summary,
            weight_mode=weight_mode,
            gate_metric=gate_metric,
        )
    else:
        raise ValueError(f"unknown identity mode: {identity_mode}")
    gene_index = {str(gene): index for index, gene in enumerate(genes)}

    all_conditions = list(dict.fromkeys(
        train_conditions + validation_conditions + target_conditions
    ))
    nodes = sorted({
        node
        for condition in all_conditions
        if condition != "control"
        for node in endpoints(condition)
    })
    train_doubles = [condition for condition in train_conditions if "+" in condition]
    validation_doubles = [condition for condition in validation_conditions if "+" in condition]

    def endpoint_effect(node: str) -> np.ndarray:
        if node in known_single:
            return known_single[node]
        if not known_single:
            return default_single
        if identity_mode == "go_shrinkage":
            return identity.predict(node, known_single, gene_index, gene2go)
        return identity.predict(node, known_single, gene_index)

    def additive(condition: str) -> np.ndarray:
        return np.sum([endpoint_effect(node) for node in endpoints(condition)], axis=0)

    if not train_doubles:
        if validation_doubles:
            alpha_scores = []
            for candidate_alpha in alpha_grid:
                val_base = np.stack([
                    saturate(additive(c), candidate_alpha)
                    for c in validation_doubles
                ])
                val_truth = np.stack([means[c] - control for c in validation_doubles])
                alpha_scores.append((
                    float(np.mean((val_base - val_truth) ** 2)),
                    candidate_alpha,
                ))
            factorized_val_mse, alpha = min(alpha_scores)
        else:
            alpha = 0.0
            factorized_val_mse = float("nan")
        factorized = {}
        for condition in target_conditions:
            if condition == "control":
                base_effect = np.zeros_like(control)
            elif "+" in condition:
                base_effect = saturate(additive(condition), alpha)
            else:
                base_effect = endpoint_effect(condition)
            factorized[condition] = control + base_effect
        return WitnessFitV2(
            predictions=dict(factorized),
            factorized_predictions=factorized,
            alpha=float(alpha),
            noise_ratio=0.0,
            gamma=0.0,
            known_single_genes=tuple(sorted(known_single)),
            training_doubles=(),
            validation_doubles=tuple(validation_doubles),
            validation_mse_factorized=float(factorized_val_mse),
            validation_mse_witness=float(factorized_val_mse),
            identity_head=identity,
        )

    train_design = incidence(train_doubles, nodes)
    validation_design = (
        incidence(validation_doubles, nodes)
        if validation_doubles
        else np.zeros((0, len(nodes)))
    )
    candidates = []
    for alpha in alpha_grid:
        train_base = np.stack([saturate(additive(c), alpha) for c in train_doubles])
        train_truth = np.stack([means[c] - control for c in train_doubles])
        residual = train_truth - train_base
        if validation_doubles:
            val_base = np.stack([saturate(additive(c), alpha) for c in validation_doubles])
            val_truth = np.stack([means[c] - control for c in validation_doubles])
            for noise in noise_grid:
                raw = kernel_predict(train_design, validation_design, residual, noise)
                gamma = closed_form_gamma(raw, val_truth - val_base)
                witness_mse = float(np.mean((val_base + gamma * raw - val_truth) ** 2))
                factorized_mse = float(np.mean((val_base - val_truth) ** 2))
                candidates.append((witness_mse, alpha, noise, gamma, factorized_mse))
        else:
            candidates.append((float("nan"), alpha, 0.1, 1.0, float("nan")))

    if validation_doubles:
        candidates.sort(key=lambda row: (row[0], row[1], row[2]))
    val_mse, alpha, noise, gamma, factorized_val_mse = candidates[0]
    train_base = np.stack([saturate(additive(c), alpha) for c in train_doubles])
    train_truth = np.stack([means[c] - control for c in train_doubles])
    train_residual = train_truth - train_base
    target_doubles = [condition for condition in target_conditions if "+" in condition]
    target_design = (
        incidence(target_doubles, nodes)
        if target_doubles
        else np.zeros((0, len(nodes)))
    )
    target_raw = (
        kernel_predict(train_design, target_design, train_residual, noise)
        if target_doubles
        else np.zeros((0, len(control)))
    )
    correction = {
        condition: target_raw[index]
        for index, condition in enumerate(target_doubles)
    }
    witness = {}
    factorized = {}
    for condition in target_conditions:
        if condition == "control":
            base_effect = np.zeros_like(control)
        elif "+" in condition:
            base_effect = saturate(additive(condition), alpha)
        else:
            base_effect = endpoint_effect(condition)
        factorized[condition] = control + base_effect
        witness[condition] = control + base_effect + gamma * correction.get(
            condition, np.zeros_like(control)
        )
    return WitnessFitV2(
        predictions=witness,
        factorized_predictions=factorized,
        alpha=float(alpha),
        noise_ratio=float(noise),
        gamma=float(gamma),
        known_single_genes=tuple(sorted(known_single)),
        training_doubles=tuple(train_doubles),
        validation_doubles=tuple(validation_doubles),
        validation_mse_factorized=float(factorized_val_mse),
        validation_mse_witness=float(val_mse),
        identity_head=identity,
    )
