#!/usr/bin/env python3
"""Shared utilities for estimated Witness Risk experiments.

The deployable model is a two-kernel Gaussian linear predictor.  The node
incidence kernel represents the factorized mechanism, while an inductive RBF
kernel over target-safe pair descriptors represents non-factorized
discrepancy.  Every hyperparameter is selected using only outer-training
double perturbations.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial.distance import cdist


EPS = 1e-12


def incidence(edges: np.ndarray, n_nodes: int) -> np.ndarray:
    matrix = np.zeros((len(edges), n_nodes), dtype=float)
    matrix[np.arange(len(edges)), edges[:, 0]] = 1.0
    matrix[np.arange(len(edges)), edges[:, 1]] = 1.0
    return matrix


def safe_edge_split(
    indices: np.ndarray,
    edges: np.ndarray,
    n_nodes: int,
    fraction: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Hold out edges while leaving every held-out endpoint in the fit graph."""
    local_edges = edges[indices]
    degree = np.bincount(local_edges.ravel(), minlength=n_nodes)
    held_local: list[int] = []
    requested = max(1, round(fraction * len(indices)))
    for local in rng.permutation(len(indices)):
        i, j = local_edges[local]
        if degree[i] > 1 and degree[j] > 1:
            held_local.append(int(local))
            degree[i] -= 1
            degree[j] -= 1
        if len(held_local) >= requested:
            break
    held_local_array = np.asarray(sorted(held_local), dtype=int)
    fit_local = np.setdiff1d(np.arange(len(indices)), held_local_array)
    return indices[fit_local], indices[held_local_array]


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    clean = np.maximum(np.asarray(weights, dtype=float), 0.0)
    if clean.sum() <= EPS:
        clean = np.ones_like(clean)
    return float(np.sum(clean * values) / np.sum(clean))


def row_cosine(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    numerator = np.sum(prediction * truth, axis=1)
    denominator = np.linalg.norm(prediction, axis=1) * np.linalg.norm(truth, axis=1)
    return numerator / (denominator + EPS)


@dataclass(frozen=True)
class KernelMatrices:
    geometry_train: np.ndarray
    geometry_cross: np.ndarray
    discrepancy_train: np.ndarray
    discrepancy_cross: np.ndarray


def standardized_squared_distances(
    train_features: np.ndarray,
    target_features: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    mean = train_features.mean(axis=0)
    scale = train_features.std(axis=0)
    scale[scale <= 1e-10] = 1.0
    train = (train_features - mean) / scale
    target = (target_features - mean) / scale
    within = cdist(train, train, metric="sqeuclidean")
    upper = within[np.triu_indices(len(train), 1)]
    reference = float(np.median(upper[upper > 0])) if np.any(upper > 0) else 1.0
    cross = cdist(train, target, metric="sqeuclidean")
    return within, cross, max(reference, EPS)


def make_kernel_matrices(
    train: np.ndarray,
    target: np.ndarray,
    design: np.ndarray,
    pair_features: np.ndarray,
    length_factor: float,
) -> KernelMatrices:
    within, cross, reference = standardized_squared_distances(
        pair_features[train], pair_features[target]
    )
    discrepancy_train = np.exp(-within / (length_factor * reference))
    discrepancy_cross = np.exp(-cross / (length_factor * reference))
    # Every unsigned pair-incidence row has squared norm two.  Division by two
    # gives unit diagonal and makes rho a genuine signal-mixture fraction.
    geometry_train = design[train] @ design[train].T / 2.0
    geometry_cross = design[train] @ design[target].T / 2.0
    return KernelMatrices(
        geometry_train=geometry_train,
        geometry_cross=geometry_cross,
        discrepancy_train=discrepancy_train,
        discrepancy_cross=discrepancy_cross,
    )


@dataclass(frozen=True)
class FitResult:
    prediction: np.ndarray
    risk: np.ndarray
    weights: np.ndarray
    scale: float
    discrepancy_variance: float
    noise_variance: float


def fit_predict(
    train_response: np.ndarray,
    kernels: KernelMatrices,
    rho: float,
    noise_ratio: float,
) -> FitResult:
    """Posterior mean and exact squared-loss predictive risk.

    With unit-diagonal component kernels, the scale MLE is

      tr(Y^T C0^{-1}Y) / (m p).

    The returned discrepancy quantities correspond to
    K_hat = scale*rho*K_phi, k_t = scale*rho*k_phi,t, and
    k_tt = scale*rho.  The total risk also contains the factorized mechanism
    kernel and target observation noise.
    """
    signal_train = (
        (1.0 - rho) * kernels.geometry_train
        + rho * kernels.discrepancy_train
    )
    signal_cross = (
        (1.0 - rho) * kernels.geometry_cross
        + rho * kernels.discrepancy_cross
    )
    covariance = signal_train + noise_ratio * np.eye(len(signal_train))
    solved_cross = np.linalg.solve(covariance, signal_cross)
    solved_response = np.linalg.solve(covariance, train_response)
    scale = float(np.sum(train_response * solved_response) / train_response.size)
    scale = max(scale, EPS)
    unit_risk = 1.0 + noise_ratio - np.sum(signal_cross * solved_cross, axis=0)
    unit_risk = np.maximum(unit_risk, EPS)
    return FitResult(
        prediction=solved_cross.T @ train_response,
        risk=scale * unit_risk,
        weights=solved_cross,
        scale=scale,
        discrepancy_variance=scale * rho,
        noise_variance=scale * noise_ratio,
    )


def weighted_mse(
    prediction: np.ndarray,
    truth: np.ndarray,
    pair_weights: np.ndarray,
) -> float:
    per_pair = np.mean((prediction - truth) ** 2, axis=1)
    return weighted_mean(per_pair, pair_weights)


def weighted_cosine(
    prediction: np.ndarray,
    truth: np.ndarray,
    pair_weights: np.ndarray,
) -> float:
    return weighted_mean(row_cosine(prediction, truth), pair_weights)


def oracle_conditional(
    full_response: np.ndarray,
    train: np.ndarray,
    target: np.ndarray,
) -> FitResult:
    """Leakage-allowed empirical oracle used only as a diagnostic ceiling.

    C = YY^T/p makes the quadratic risk exactly equal to the mean squared
    error across output genes.  A tiny eigenvalue floor is added only for
    numerical stability.
    """
    second_moment = full_response @ full_response.T / full_response.shape[1]
    train_covariance = second_moment[np.ix_(train, train)]
    floor = max(float(np.mean(np.diag(train_covariance))) * 1e-10, 1e-12)
    train_covariance = train_covariance + floor * np.eye(len(train))
    cross = second_moment[np.ix_(train, target)]
    target_variance = np.diag(second_moment[np.ix_(target, target)]) + floor
    solved_cross = np.linalg.solve(train_covariance, cross)
    risk = target_variance - np.sum(cross * solved_cross, axis=0)
    return FitResult(
        prediction=solved_cross.T @ full_response[train],
        risk=np.maximum(risk, EPS),
        weights=solved_cross,
        scale=1.0,
        discrepancy_variance=float("nan"),
        noise_variance=floor,
    )

