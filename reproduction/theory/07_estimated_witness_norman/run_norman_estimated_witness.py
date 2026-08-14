#!/usr/bin/env python3
"""Nested-CV test of estimated Witness Risk on Norman double perturbations."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from estimated_witness import (
    EPS,
    fit_predict,
    incidence,
    make_kernel_matrices,
    oracle_conditional,
    row_cosine,
    safe_edge_split,
    weighted_cosine,
    weighted_mean,
    weighted_mse,
)


def safe_spearman(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 4 or np.ptp(left) <= EPS or np.ptp(right) <= EPS:
        return float("nan")
    return float(spearmanr(left, right).statistic)


def load_single_profiles(data_root: Path, output_genes: np.ndarray, nodes: list[str]) -> np.ndarray:
    mean = np.load(data_root / "pseudobulk_mean.npy")
    perturbations = json.loads((data_root / "perts.json").read_text())
    perturbation_id = {name: index for index, name in enumerate(perturbations)}
    gene_names = np.asarray([
        (line.split("\t") + [""])[1] or line.split("\t")[0]
        for line in gzip.open(data_root / "GSE133344_filtered_genes.tsv.gz", "rt")
        .read().splitlines()
    ])
    gene_id = {name: index for index, name in enumerate(gene_names)}
    output_index = np.asarray([gene_id[name] for name in output_genes], dtype=int)
    effect = mean[:, output_index] - mean[perturbation_id["CTRL"], output_index]
    return np.vstack([effect[perturbation_id[node]] for node in nodes])


def symmetric_pair_features(single_profiles: np.ndarray, edges: np.ndarray, components: int) -> np.ndarray:
    standardized = StandardScaler().fit_transform(single_profiles)
    latent = PCA(n_components=min(components, len(single_profiles) - 1), random_state=0).fit_transform(
        standardized
    )
    upper = np.triu_indices(latent.shape[1])
    rows = []
    for left, right in edges:
        outer = 0.5 * (
            np.outer(latent[left], latent[right])
            + np.outer(latent[right], latent[left])
        )
        rows.append(np.r_[
            latent[left] + latent[right],
            np.abs(latent[left] - latent[right]),
            outer[upper],
        ])
    return np.asarray(rows, dtype=float)


def inner_splits(
    train: np.ndarray,
    edges: np.ndarray,
    n_nodes: int,
    count: int,
    fraction: float,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    return [
        safe_edge_split(
            train, edges, n_nodes, fraction,
            np.random.default_rng(seed + 1009 * repeat),
        )
        for repeat in range(count)
    ]


def score_candidate(
    splits: list[tuple[np.ndarray, np.ndarray]],
    response: np.ndarray,
    reliability: np.ndarray,
    design: np.ndarray,
    features: np.ndarray,
    length_factor: float,
    rho: float,
    noise_ratio: float,
) -> float:
    scores = []
    for fit, validation in splits:
        kernels = make_kernel_matrices(
            fit, validation, design, features, length_factor
        )
        result = fit_predict(response[fit], kernels, rho, noise_ratio)
        scores.append(weighted_mse(
            result.prediction, response[validation], reliability[validation]
        ))
    return float(np.mean(scores))


def select_parameters(
    splits: list[tuple[np.ndarray, np.ndarray]],
    response: np.ndarray,
    reliability: np.ndarray,
    design: np.ndarray,
    features: np.ndarray,
    length_grid: list[float],
    rho_grid: list[float],
    noise_grid: list[float],
) -> tuple[dict, pd.DataFrame]:
    rows = []
    for rho in rho_grid:
        local_lengths = [length_grid[0]] if rho == 0.0 else length_grid
        for length_factor in local_lengths:
            for noise_ratio in noise_grid:
                score = score_candidate(
                    splits, response, reliability, design, features,
                    length_factor, rho, noise_ratio,
                )
                rows.append({
                    "length_factor": length_factor,
                    "rho": rho,
                    "noise_ratio": noise_ratio,
                    "inner_weighted_mse": score,
                })
    frame = pd.DataFrame(rows).sort_values(
        ["inner_weighted_mse", "rho", "noise_ratio", "length_factor"]
    ).reset_index(drop=True)
    return frame.iloc[0].to_dict(), frame


def evaluate_prediction(
    prediction: np.ndarray,
    residual_truth: np.ndarray,
    additive: np.ndarray,
    double_effect: np.ndarray,
    reliability: np.ndarray,
) -> dict:
    residual_energy = weighted_mean(np.mean(residual_truth ** 2, axis=1), reliability)
    full_truth_energy = weighted_mean(np.mean(double_effect ** 2, axis=1), reliability)
    return {
        "residual_mse": weighted_mse(prediction, residual_truth, reliability),
        "residual_nmse": weighted_mse(prediction, residual_truth, reliability)
        / (residual_energy + EPS),
        "residual_cosine": weighted_cosine(prediction, residual_truth, reliability),
        "full_effect_mse": weighted_mse(additive + prediction, double_effect, reliability),
        "full_effect_nmse": weighted_mse(additive + prediction, double_effect, reliability)
        / (full_truth_energy + EPS),
        "full_effect_cosine": weighted_cosine(additive + prediction, double_effect, reliability),
    }


def paired_test(frame: pd.DataFrame, metric: str, higher_better: bool) -> dict:
    wide = frame.pivot(index="seed", columns="strategy", values=metric)
    if higher_better:
        difference = wide["estimated_witness"] - wide["geometry_only"]
        relative = difference / np.maximum(np.abs(wide["geometry_only"]), EPS)
    else:
        difference = wide["geometry_only"] - wide["estimated_witness"]
        relative = difference / np.maximum(np.abs(wide["geometry_only"]), EPS)
    return {
        "metric": metric,
        "direction": "higher" if higher_better else "lower",
        "difference_mean": float(difference.mean()),
        "relative_improvement_mean": float(relative.mean()),
        "wins": int((difference > 0).sum()),
        "win_rate": float((difference > 0).mean()),
        "p_one_sided_wilcoxon": float(wilcoxon(difference, alternative="greater").pvalue)
        if not np.allclose(difference, 0) else 1.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--residuals", type=Path, default=Path("../04_formal_gate/results/pseudobulk_residuals.npz"))
    parser.add_argument("--data-root", type=Path, default=Path("../../data/norman"))
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--seed-offset", type=int, default=0)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--inner-splits", type=int, default=3)
    parser.add_argument("--inner-fraction", type=float, default=0.18)
    parser.add_argument("--pca-components", type=int, default=12)
    parser.add_argument("--out", type=Path, default=Path("results/smoke"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    data = np.load(args.residuals, allow_pickle=True)
    pairs = data["pairs"].astype(str)
    response = data["residual"].astype(float)
    additive = data["additive"].astype(float)
    double_effect = data["double_effect"].astype(float)
    reliability = data["weights"].astype(float)
    output_genes = data["genes"].astype(str)
    nodes = sorted({node for pair in pairs for node in pair.split("+")})
    node_id = {node: index for index, node in enumerate(nodes)}
    edges = np.asarray([
        (node_id[pair.split("+")[0]], node_id[pair.split("+")[1]]) for pair in pairs
    ], dtype=int)
    design = incidence(edges, len(nodes))
    single_profiles = load_single_profiles(args.data_root, output_genes, nodes)
    pair_features = symmetric_pair_features(single_profiles, edges, args.pca_components)

    length_grid = [0.25, 0.5, 1.0, 2.0, 4.0]
    rho_grid = [0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0]
    noise_grid = [0.03, 0.10, 0.30, 1.0, 3.0, 10.0]
    all_indices = np.arange(len(edges))
    metric_rows: list[dict] = []
    target_rows: list[pd.DataFrame] = []
    tuning_rows: list[pd.DataFrame] = []

    for local_seed in range(args.seeds):
        seed = args.seed_offset + local_seed
        train, test = safe_edge_split(
            all_indices, edges, len(nodes), args.test_fraction,
            np.random.default_rng(7411 + seed),
        )
        nested = inner_splits(
            train, edges, len(nodes), args.inner_splits, args.inner_fraction,
            99173 + seed,
        )
        witness_parameters, witness_tuning = select_parameters(
            nested, response, reliability, design, pair_features,
            length_grid, rho_grid, noise_grid,
        )
        geometry_parameters, geometry_tuning = select_parameters(
            nested, response, reliability, design, pair_features,
            [length_grid[0]], [0.0], noise_grid,
        )
        for name, tuning in (("estimated_witness", witness_tuning), ("geometry_only", geometry_tuning)):
            local = tuning.copy()
            local.insert(0, "seed", seed)
            local.insert(1, "strategy", name)
            tuning_rows.append(local)

        fitted = {}
        fitted_kernels = {}
        for strategy, parameters in (
            ("estimated_witness", witness_parameters),
            ("geometry_only", geometry_parameters),
        ):
            kernels = make_kernel_matrices(
                train, test, design, pair_features,
                float(parameters["length_factor"]),
            )
            result = fit_predict(
                response[train], kernels,
                float(parameters["rho"]), float(parameters["noise_ratio"]),
            )
            fitted[strategy] = result
            fitted_kernels[strategy] = kernels
            metrics = evaluate_prediction(
                result.prediction, response[test], additive[test],
                double_effect[test], reliability[test],
            )
            metric_rows.append({
                "seed": seed, "strategy": strategy, **metrics,
                "selected_length_factor": float(parameters["length_factor"]),
                "selected_rho": float(parameters["rho"]),
                "selected_noise_ratio": float(parameters["noise_ratio"]),
                "covariance_scale": result.scale,
                "estimated_K_diagonal": result.discrepancy_variance,
                "estimated_noise_variance": result.noise_variance,
            })

        oracle = oracle_conditional(response, train, test)
        oracle_metrics = evaluate_prediction(
            oracle.prediction, response[test], additive[test],
            double_effect[test], reliability[test],
        )
        metric_rows.append({
            "seed": seed, "strategy": "empirical_oracle", **oracle_metrics,
            "selected_length_factor": np.nan, "selected_rho": np.nan,
            "selected_noise_ratio": np.nan, "covariance_scale": 1.0,
            "estimated_K_diagonal": np.nan,
            "estimated_noise_variance": oracle.noise_variance,
        })

        estimated = fitted["estimated_witness"]
        geometry = fitted["geometry_only"]
        estimated_kernels = fitted_kernels["estimated_witness"]
        covariance_dir = args.out / "estimated_covariances"
        covariance_dir.mkdir(exist_ok=True)
        estimated_rho = float(witness_parameters["rho"])
        np.savez_compressed(
            covariance_dir / f"seed_{seed:03d}.npz",
            train_indices=train,
            test_indices=test,
            train_pairs=pairs[train],
            test_pairs=pairs[test],
            K_hat=(
                estimated.scale * estimated_rho
                * estimated_kernels.discrepancy_train
            ).astype(np.float32),
            k_t_hat=(
                estimated.scale * estimated_rho
                * estimated_kernels.discrepancy_cross
            ).astype(np.float32),
            k_tt_hat=np.full(
                len(test), estimated.scale * estimated_rho, dtype=np.float32
            ),
            noise_variance=np.asarray(estimated.noise_variance),
            witness_risk=estimated.risk.astype(np.float32),
            length_factor=np.asarray(float(witness_parameters["length_factor"])),
            rho=np.asarray(estimated_rho),
            noise_ratio=np.asarray(float(witness_parameters["noise_ratio"])),
            covariance_scale=np.asarray(estimated.scale),
        )
        realized_estimated = np.mean((estimated.prediction - response[test]) ** 2, axis=1)
        realized_geometry = np.mean((geometry.prediction - response[test]) ** 2, axis=1)
        realized_oracle = np.mean((oracle.prediction - response[test]) ** 2, axis=1)
        target_rows.append(pd.DataFrame({
            "seed": seed,
            "pair": pairs[test],
            "reliability_weight": reliability[test],
            "estimated_witness_risk": estimated.risk,
            "geometry_risk": geometry.risk,
            "oracle_risk": oracle.risk,
            "estimated_realized_mse": realized_estimated,
            "geometry_realized_mse": realized_geometry,
            "oracle_realized_mse": realized_oracle,
            "estimated_residual_cosine": row_cosine(estimated.prediction, response[test]),
            "geometry_residual_cosine": row_cosine(geometry.prediction, response[test]),
            "selected_rho": float(witness_parameters["rho"]),
            "selected_length_factor": float(witness_parameters["length_factor"]),
            "selected_noise_ratio": float(witness_parameters["noise_ratio"]),
        }))

    metrics = pd.DataFrame(metric_rows)
    targets = pd.concat(target_rows, ignore_index=True)
    tuning = pd.concat(tuning_rows, ignore_index=True)
    metrics.to_csv(args.out / "per_seed_metrics.csv", index=False)
    targets.to_csv(args.out / "target_rows.csv", index=False)
    tuning.to_csv(args.out / "inner_cv_grid.csv", index=False)

    metric_columns = [
        "residual_mse", "residual_nmse", "residual_cosine",
        "full_effect_mse", "full_effect_nmse", "full_effect_cosine",
        "selected_rho", "selected_noise_ratio",
    ]
    summary = metrics.groupby("strategy")[metric_columns].agg(["mean", "std"])
    summary.columns = ["_".join(column) for column in summary.columns]
    summary.reset_index().to_csv(args.out / "summary.csv", index=False)

    tests = [
        paired_test(metrics, "residual_mse", False),
        paired_test(metrics, "residual_nmse", False),
        paired_test(metrics, "residual_cosine", True),
        paired_test(metrics, "full_effect_mse", False),
        paired_test(metrics, "full_effect_nmse", False),
        paired_test(metrics, "full_effect_cosine", True),
    ]
    pd.DataFrame(tests).to_csv(args.out / "paired_tests.csv", index=False)

    risk_rows = []
    for seed, local in targets.groupby("seed"):
        risk_rows.append({
            "seed": int(seed),
            "spearman_estimated_risk_vs_oracle_risk": safe_spearman(
                local.estimated_witness_risk.to_numpy(), local.oracle_risk.to_numpy()
            ),
            "spearman_estimated_risk_vs_realized_mse": safe_spearman(
                local.estimated_witness_risk.to_numpy(), local.estimated_realized_mse.to_numpy()
            ),
            "spearman_geometry_risk_vs_realized_mse": safe_spearman(
                local.geometry_risk.to_numpy(), local.geometry_realized_mse.to_numpy()
            ),
            "oracle_formula_relative_error": float(np.max(
                np.abs(local.oracle_risk - local.oracle_realized_mse)
                / np.maximum(local.oracle_realized_mse, EPS)
            )),
            "oracle_efficiency_estimated": weighted_mean(
                local.oracle_risk.to_numpy() / np.maximum(local.estimated_realized_mse.to_numpy(), EPS),
                local.reliability_weight.to_numpy(),
            ),
            "oracle_efficiency_geometry": weighted_mean(
                local.oracle_risk.to_numpy() / np.maximum(local.geometry_realized_mse.to_numpy(), EPS),
                local.reliability_weight.to_numpy(),
            ),
        })
    risk_frame = pd.DataFrame(risk_rows)
    risk_frame.to_csv(args.out / "risk_alignment.csv", index=False)
    risk_summary = risk_frame.agg(["mean", "std", "median"]).T.reset_index()
    risk_summary.columns = ["metric", "mean", "std", "median"]
    risk_summary.to_csv(args.out / "risk_alignment_summary.csv", index=False)

    test_map = {row["metric"]: row for row in tests}
    residual_test = test_map["residual_mse"]
    full_test = test_map["full_effect_mse"]
    risk_oracle_mean = float(risk_frame.spearman_estimated_risk_vs_oracle_risk.mean())
    risk_realized_mean = float(risk_frame.spearman_estimated_risk_vs_realized_mse.mean())
    gate = {
        "claim": "K, k_t and k_tt estimated only from training doubles plus target-safe single-perturbation descriptors yield a deployable Witness Risk and improve held-out Norman prediction.",
        "scope": "Norman development-panel pseudobulk additive residual; nested outer-edge holdout; empirical oracle is diagnostic and never used for fitting or selection",
        "leakage_audit": {
            "outer_test_double_outcomes_used_for_hyperparameter_selection": False,
            "outer_test_double_outcomes_used_for_estimated_K_kt_ktt": False,
            "test_single_perturbation_profiles_used": True,
            "test_reliability_weights_used_only_for_evaluation": True,
        },
        "frozen_formal_gate": {
            "risk": "mean within-split Spearman(estimated risk, oracle risk) >=0.20 and Spearman(estimated risk, realized estimated-predictor MSE) >=0.25",
            "prediction": "mean residual-MSE improvement >=3%, win rate >=70%, one-sided paired Wilcoxon p<0.05; mean full-effect MSE noninferior",
        },
        "risk_spearman_oracle_mean": risk_oracle_mean,
        "risk_spearman_realized_mean": risk_realized_mean,
        "prediction_test": residual_test,
        "full_effect_safety_test": full_test,
        "selected_rho_mean": float(
            metrics.loc[metrics.strategy == "estimated_witness", "selected_rho"].mean()
        ),
    }
    gate["risk_gate_pass"] = bool(
        risk_oracle_mean >= 0.20 and risk_realized_mean >= 0.25
    )
    gate["prediction_gate_pass"] = bool(
        residual_test["relative_improvement_mean"] >= 0.03
        and residual_test["win_rate"] >= 0.70
        and residual_test["p_one_sided_wilcoxon"] < 0.05
        and full_test["difference_mean"] >= 0.0
    )
    gate["gate_pass"] = bool(gate["risk_gate_pass"] and gate["prediction_gate_pass"])
    (args.out / "verdict.json").write_text(json.dumps(gate, indent=2))
    print(summary.to_string())
    print("\n" + pd.DataFrame(tests).to_string(index=False))
    print("\n" + risk_summary.to_string(index=False))
    print("\n" + json.dumps(gate, indent=2))
    if args.seeds >= 30 and not gate["gate_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
