#!/usr/bin/env python3
"""Known-covariance positive and negative controls for the estimator."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from estimated_witness import (
    EPS,
    fit_predict,
    incidence,
    make_kernel_matrices,
    safe_edge_split,
    weighted_mse,
)
from run_norman_estimated_witness import inner_splits, select_parameters


def pair_features(node_features: np.ndarray, edges: np.ndarray) -> np.ndarray:
    upper = np.triu_indices(node_features.shape[1])
    rows = []
    for left, right in edges:
        outer = 0.5 * (
            np.outer(node_features[left], node_features[right])
            + np.outer(node_features[right], node_features[left])
        )
        rows.append(np.r_[
            node_features[left] + node_features[right],
            np.abs(node_features[left] - node_features[right]),
            outer[upper],
        ])
    return np.asarray(rows)


def true_conditional(
    covariance: np.ndarray,
    response: np.ndarray,
    train: np.ndarray,
    test: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    train_covariance = covariance[np.ix_(train, train)]
    cross = covariance[np.ix_(train, test)]
    solved = np.linalg.solve(train_covariance, cross)
    risk = np.diag(covariance[np.ix_(test, test)]) - np.sum(cross * solved, axis=0)
    return solved.T @ response[train], risk


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=12)
    parser.add_argument("--outputs", type=int, default=512)
    parser.add_argument("--out", type=Path, default=Path("results/estimator_controls"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    n_nodes = 12
    all_edges = np.asarray([
        (left, right) for left in range(n_nodes) for right in range(left + 1, n_nodes)
    ], dtype=int)
    design = incidence(all_edges, n_nodes)
    reliability = np.ones(len(all_edges))
    rows = []
    for regime, true_rho in (("positive", 0.70), ("negative", 0.0)):
        for seed in range(args.seeds):
            rng = np.random.default_rng(20261100 + 97 * seed + int(1000 * true_rho))
            node_features = rng.normal(size=(n_nodes, 6))
            features = pair_features(node_features, all_edges)
            all_kernel = make_kernel_matrices(
                np.arange(len(all_edges)), np.arange(len(all_edges)),
                design, features, length_factor=1.0,
            )
            noise_ratio = 0.30
            covariance = (
                (1.0 - true_rho) * all_kernel.geometry_train
                + true_rho * all_kernel.discrepancy_train
                + noise_ratio * np.eye(len(all_edges))
            )
            response = rng.multivariate_normal(
                np.zeros(len(all_edges)), covariance, size=args.outputs
            ).T
            train, test = safe_edge_split(
                np.arange(len(all_edges)), all_edges, n_nodes, 0.20,
                np.random.default_rng(31001 + seed),
            )
            nested = inner_splits(
                train, all_edges, n_nodes, 3, 0.20, 45001 + seed
            )
            parameters, _ = select_parameters(
                nested, response, reliability, design, features,
                [0.5, 1.0, 2.0], [0.0, 0.25, 0.50, 0.75, 1.0],
                [0.10, 0.30, 1.0],
            )
            baseline_parameters, _ = select_parameters(
                nested, response, reliability, design, features,
                [1.0], [0.0], [0.10, 0.30, 1.0],
            )
            kernels = make_kernel_matrices(
                train, test, design, features, float(parameters["length_factor"])
            )
            estimated = fit_predict(
                response[train], kernels, float(parameters["rho"]),
                float(parameters["noise_ratio"]),
            )
            baseline_kernels = make_kernel_matrices(
                train, test, design, features,
                float(baseline_parameters["length_factor"]),
            )
            baseline = fit_predict(
                response[train], baseline_kernels, 0.0,
                float(baseline_parameters["noise_ratio"]),
            )
            oracle_prediction, true_risk = true_conditional(
                covariance, response, train, test
            )
            estimated_mse = weighted_mse(
                estimated.prediction, response[test], reliability[test]
            )
            baseline_mse = weighted_mse(
                baseline.prediction, response[test], reliability[test]
            )
            oracle_mse = weighted_mse(
                oracle_prediction, response[test], reliability[test]
            )
            rows.append({
                "regime": regime,
                "seed": seed,
                "true_rho": true_rho,
                "selected_rho": float(parameters["rho"]),
                "selected_length_factor": float(parameters["length_factor"]),
                "selected_noise_ratio": float(parameters["noise_ratio"]),
                "estimated_mse": estimated_mse,
                "geometry_mse": baseline_mse,
                "oracle_mse": oracle_mse,
                "relative_gain_over_geometry": (baseline_mse - estimated_mse)
                / max(baseline_mse, EPS),
                "risk_spearman_true": float(spearmanr(estimated.risk, true_risk).statistic),
            })
    frame = pd.DataFrame(rows)
    frame.to_csv(args.out / "control_rows.csv", index=False)
    summary = frame.groupby("regime").agg(
        seeds=("seed", "nunique"),
        selected_rho_median=("selected_rho", "median"),
        selected_rho_mean=("selected_rho", "mean"),
        relative_gain_mean=("relative_gain_over_geometry", "mean"),
        risk_spearman_mean=("risk_spearman_true", "mean"),
        estimated_mse_mean=("estimated_mse", "mean"),
        geometry_mse_mean=("geometry_mse", "mean"),
        oracle_mse_mean=("oracle_mse", "mean"),
    ).reset_index()
    summary.to_csv(args.out / "summary.csv", index=False)
    positive = summary[summary.regime == "positive"].iloc[0]
    negative = summary[summary.regime == "negative"].iloc[0]
    verdict = {
        "positive_control": "known discrepancy kernel",
        "negative_control": "factorized geometry plus noise only",
        "gate": "positive median rho >=0.50, gain >=5%, risk Spearman >=0.50; negative median rho <=0.10 and no >5% spurious gain",
        "positive": positive.to_dict(),
        "negative": negative.to_dict(),
    }
    verdict["gate_pass"] = bool(
        positive.selected_rho_median >= 0.50
        and positive.relative_gain_mean >= 0.05
        and positive.risk_spearman_mean >= 0.50
        and negative.selected_rho_median <= 0.10
        and negative.relative_gain_mean <= 0.05
    )
    (args.out / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(summary.to_string(index=False))
    print("\n" + json.dumps(verdict, indent=2))
    if not verdict["gate_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

