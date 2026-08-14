#!/usr/bin/env python3
"""Core of the scPerturBench WitnessCell adapter.

The model uses training single perturbations as a factorized backbone and a
regularized incidence-kernel predictor over training-double residuals.  Only
official validation outcomes select saturation, kernel noise, and the clipped
Witness amplitude.  Test outcomes are never passed to this module.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


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
    # Divide by two so double-double self-similarity is one.  For a target
    # sharing one endpoint with a training pair, the kernel entry is 1/2.
    train_kernel = train_design @ train_design.T / 2.0
    cross_kernel = train_design @ target_design.T / 2.0
    covariance = train_kernel + noise_ratio * np.eye(len(train_kernel))
    return np.linalg.solve(covariance, cross_kernel).T @ response


@dataclass(frozen=True)
class WitnessFit:
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


def fit_predict(
    means: dict[str, np.ndarray],
    train_conditions: list[str],
    validation_conditions: list[str],
    target_conditions: list[str],
    alpha_grid: tuple[float, ...] = (0.0, 0.02, 0.05, 0.1, 0.2, 0.5),
    noise_grid: tuple[float, ...] = (0.01, 0.03, 0.1, 0.3, 1.0),
) -> WitnessFit:
    control = np.asarray(means["control"], dtype=float)
    train_set = set(train_conditions)
    known_single = {
        condition: np.asarray(means[condition], dtype=float) - control
        for condition in train_conditions
        if condition != "control" and "+" not in condition
    }
    if known_single:
        default_single = np.mean(np.stack(list(known_single.values())), axis=0)
    else:
        default_single = np.zeros_like(control)

    all_conditions = list(dict.fromkeys(
        train_conditions + validation_conditions + target_conditions
    ))
    nodes = sorted({node for condition in all_conditions if condition != "control" for node in endpoints(condition)})
    train_doubles = [condition for condition in train_conditions if "+" in condition]
    validation_doubles = [condition for condition in validation_conditions if "+" in condition]

    def endpoint_effect(node: str) -> np.ndarray:
        return known_single.get(node, default_single)

    def additive(condition: str) -> np.ndarray:
        return np.sum([endpoint_effect(node) for node in endpoints(condition)], axis=0)

    # A valid official split can contain no training double (Replogle_exp6,
    # seed 1).  With no intervention witness, the identifiable member of this
    # model family is its factorized backbone: gamma must be zero.  Alpha may
    # still be selected from validation doubles if present; otherwise the
    # predeclared neutral default is the additive alpha=0 model.
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

        factorized: dict[str, np.ndarray] = {}
        for condition in target_conditions:
            if condition == "control":
                base_effect = np.zeros_like(control)
            elif "+" in condition:
                base_effect = saturate(additive(condition), alpha)
            else:
                base_effect = endpoint_effect(condition)
            factorized[condition] = control + base_effect
        return WitnessFit(
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
        )

    train_design = incidence(train_doubles, nodes)
    validation_design = incidence(validation_doubles, nodes) if validation_doubles else np.zeros((0, len(nodes)))

    candidates: list[tuple[float, float, float, float, float]] = []
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
            # No validation double is an allowed but explicitly degraded mode.
            candidates.append((float("nan"), alpha, 0.1, 1.0, float("nan")))

    if validation_doubles:
        candidates.sort(key=lambda row: (row[0], row[1], row[2]))
        val_mse, alpha, noise, gamma, factorized_val_mse = candidates[0]
    else:
        val_mse, alpha, noise, gamma, factorized_val_mse = candidates[0]

    train_base = np.stack([saturate(additive(c), alpha) for c in train_doubles])
    train_truth = np.stack([means[c] - control for c in train_doubles])
    train_residual = train_truth - train_base
    target_doubles = [condition for condition in target_conditions if "+" in condition]
    target_design = incidence(target_doubles, nodes) if target_doubles else np.zeros((0, len(nodes)))
    target_raw = (
        kernel_predict(train_design, target_design, train_residual, noise)
        if target_doubles
        else np.zeros((0, len(control)))
    )
    correction = {condition: target_raw[index] for index, condition in enumerate(target_doubles)}

    witness: dict[str, np.ndarray] = {}
    factorized: dict[str, np.ndarray] = {}
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

    return WitnessFit(
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
    )
