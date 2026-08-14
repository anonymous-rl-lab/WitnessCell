#!/usr/bin/env python3
"""Monte Carlo verification of the exact geometry–adequacy BLUP risk."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def incidence(edges: np.ndarray, n_nodes: int) -> np.ndarray:
    x = np.zeros((len(edges), n_nodes))
    x[np.arange(len(edges)), edges[:, 0]] = 1.0
    x[np.arange(len(edges)), edges[:, 1]] = 1.0
    return x


def exact_blup(X: np.ndarray, target: np.ndarray, sigma: np.ndarray,
               cross: np.ndarray, target_variance: float):
    inverse = np.linalg.inv(sigma)
    information = X.T @ inverse @ X
    information_inverse = np.linalg.inv(information)
    effective_target = target - X.T @ inverse @ cross
    weights = inverse @ cross + inverse @ X @ information_inverse @ effective_target
    adequacy = float(target_variance - cross @ inverse @ cross)
    geometry = float(effective_target @ information_inverse @ effective_target)
    risk = adequacy + geometry
    direct_risk = float(weights @ sigma @ weights - 2 * weights @ cross + target_variance)
    return weights, {
        "adequacy_term": adequacy,
        "geometry_term": geometry,
        "formula_risk": risk,
        "direct_quadratic_risk": direct_risk,
    }


def rbf_covariance(features: np.ndarray, amplitude: float, lengthscale: float,
                   nugget: float) -> np.ndarray:
    sq = np.sum((features[:, None, :] - features[None, :, :]) ** 2, axis=2)
    kernel = amplitude ** 2 * np.exp(-sq / (2 * lengthscale ** 2))
    kernel += nugget ** 2 * np.eye(len(features))
    return kernel


def null_basis(matrix: np.ndarray, tolerance: float = 1e-10) -> np.ndarray:
    _u, singular, vh = np.linalg.svd(matrix, full_matrices=True)
    rank = int(np.sum(singular > tolerance))
    return vh[rank:].T


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=100000)
    parser.add_argument("--random-competitors", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--out", type=Path, default=Path("results/exact_risk"))
    args = parser.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    n_nodes = 7
    train_edges = np.asarray([
        (0, 1), (1, 2), (2, 0), (2, 3), (3, 4),
        (4, 5), (5, 6), (0, 3), (1, 5), (3, 6),
    ], dtype=int)
    target_edge = np.asarray([[0, 6]], dtype=int)
    X = incidence(train_edges, n_nodes)
    target = incidence(target_edge, n_nodes)[0]
    assert np.linalg.matrix_rank(X) == n_nodes
    node_embedding = rng.normal(size=(n_nodes, 3))
    all_edges = np.vstack([train_edges, target_edge])
    pair_features = node_embedding[all_edges[:, 0]] + node_embedding[all_edges[:, 1]]
    noise_variance = .20 ** 2

    regimes = {}
    regimes["factorized"] = np.zeros((len(all_edges), len(all_edges)))
    regimes["independent_mismatch"] = .70 ** 2 * np.eye(len(all_edges))
    regimes["correlated_mismatch"] = rbf_covariance(
        pair_features, amplitude=.70, lengthscale=2.0, nugget=.08
    )

    rows = []
    all_pass = True
    for regime, latent_covariance in regimes.items():
        train_covariance = latent_covariance[:-1, :-1]
        cross = latent_covariance[:-1, -1]
        target_variance = float(latent_covariance[-1, -1])
        sigma = train_covariance + noise_variance * np.eye(len(train_edges))
        weights, terms = exact_blup(X, target, sigma, cross, target_variance)
        joint = np.block([
            [sigma, cross[:, None]],
            [cross[None, :], np.asarray([[target_variance]])],
        ])
        draws = rng.multivariate_normal(
            np.zeros(len(train_edges) + 1), joint, size=args.replicates,
            check_valid="raise", tol=1e-9,
        )
        empirical_error = draws[:, :-1] @ weights - draws[:, -1]
        empirical_mse = float(np.mean(empirical_error ** 2))
        relative_formula_error = abs(empirical_mse - terms["formula_risk"]) / terms["formula_risk"]

        geometry_weights = X @ np.linalg.inv(X.T @ X) @ target
        geometry_risk = float(
            geometry_weights @ sigma @ geometry_weights
            - 2 * geometry_weights @ cross + target_variance
        )
        basis = null_basis(X.T)
        random_risks = []
        for _ in range(args.random_competitors):
            competitor = weights + basis @ rng.normal(scale=.5, size=basis.shape[1])
            random_risks.append(float(
                competitor @ sigma @ competitor - 2 * competitor @ cross + target_variance
            ))
        minimum_random_excess = float(min(random_risks) - terms["formula_risk"])
        unbiasedness_error = float(np.max(np.abs(X.T @ weights - target)))
        regime_pass = bool(
            relative_formula_error <= .03
            and abs(terms["formula_risk"] - terms["direct_quadratic_risk"]) <= 1e-9
            and terms["formula_risk"] <= geometry_risk + 1e-10
            and minimum_random_excess >= -1e-9
            and unbiasedness_error <= 1e-9
        )
        all_pass &= regime_pass
        rows.append({
            "regime": regime, **terms,
            "empirical_mse": empirical_mse,
            "relative_formula_error": relative_formula_error,
            "pure_geometry_predictor_risk": geometry_risk,
            "blup_gain_over_geometry": geometry_risk - terms["formula_risk"],
            "minimum_random_unbiased_excess": minimum_random_excess,
            "unbiasedness_error": unbiasedness_error,
            "gate_pass": regime_pass,
        })
    frame = pd.DataFrame(rows); frame.to_csv(args.out / "risk_verification.csv", index=False)
    verdict = {
        "theorem_checked": "exact BLUP risk and minimum risk among linear unbiased predictors",
        "replicates_per_regime": args.replicates,
        "random_unbiased_competitors_per_regime": args.random_competitors,
        "gate": "Monte Carlo relative error <=3%, algebra equality and unbiasedness <=1e-9, exact predictor no worse than pure geometry or any sampled unbiased competitor",
        "gate_pass": bool(all_pass),
    }
    (args.out / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(frame.to_string(index=False))
    print("\n" + json.dumps(verdict, indent=2))
    if not verdict["gate_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
