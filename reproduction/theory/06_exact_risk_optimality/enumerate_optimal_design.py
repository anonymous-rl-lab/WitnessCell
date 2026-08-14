#!/usr/bin/env python3
"""Enumerate candidate interventions and certify exact target-risk optimality."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from verify_exact_risk import exact_blup, incidence


def build_instance(seed: int, n_nodes: int, base_size: int,
                   candidate_size: int, target_size: int):
    rng = np.random.default_rng(20261101 + seed)
    edges = np.asarray(list(itertools.combinations(range(n_nodes), 2)), dtype=int)
    for _ in range(1000):
        permutation = rng.permutation(len(edges))
        base = permutation[:base_size]
        if np.linalg.matrix_rank(incidence(edges[base], n_nodes)) == n_nodes:
            break
    else:
        raise RuntimeError("could not construct full-rank base design")
    candidates = permutation[base_size:base_size + candidate_size]
    targets = permutation[base_size + candidate_size:base_size + candidate_size + target_size]
    node_embedding = rng.normal(size=(n_nodes, 3))
    node_scale = rng.lognormal(mean=-.1, sigma=.8, size=n_nodes)
    pair_features = node_embedding[edges[:, 0]] + node_embedding[edges[:, 1]]
    pair_amplitude = .55 * np.sqrt(node_scale[edges[:, 0]] * node_scale[edges[:, 1]])
    sq = np.sum((pair_features[:, None, :] - pair_features[None, :, :]) ** 2, axis=2)
    correlation = np.exp(-sq / (2 * 1.5 ** 2))
    covariance = pair_amplitude[:, None] * correlation * pair_amplitude[None, :]
    covariance += .10 ** 2 * np.eye(len(edges))
    return rng, edges, base, candidates, targets, covariance


def design_risks(edges: np.ndarray, observed: np.ndarray, targets: np.ndarray,
                 covariance: np.ndarray, n_nodes: int, noise_variance: float):
    X = incidence(edges[observed], n_nodes)
    sigma = covariance[np.ix_(observed, observed)] + noise_variance * np.eye(len(observed))
    pure_information_inverse = np.linalg.inv(X.T @ X)
    exact_total = 0.0; geometry_total = 0.0; adequacy_total = 0.0
    weights = []
    for target_index in targets:
        target = incidence(edges[[target_index]], n_nodes)[0]
        cross = covariance[observed, target_index]
        target_variance = float(covariance[target_index, target_index])
        a, terms = exact_blup(X, target, sigma, cross, target_variance)
        exact_total += terms["formula_risk"]
        geometry_total += noise_variance * float(target @ pure_information_inverse @ target)
        adequacy_total += terms["adequacy_term"]
        weights.append(a)
    return {
        "exact": exact_total / len(targets),
        "pure_geometry": geometry_total / len(targets),
        "adequacy_only": adequacy_total / len(targets),
    }, np.asarray(weights), sigma


def empirical_mse(rng: np.random.Generator, edges: np.ndarray, observed: np.ndarray,
                  targets: np.ndarray, covariance: np.ndarray,
                  noise_variance: float, weights: np.ndarray,
                  replicates: int) -> float:
    sigma = covariance[np.ix_(observed, observed)] + noise_variance * np.eye(len(observed))
    cross = covariance[np.ix_(observed, targets)]
    target_covariance = covariance[np.ix_(targets, targets)]
    joint = np.block([[sigma, cross], [cross.T, target_covariance]])
    draws = rng.multivariate_normal(
        np.zeros(len(observed) + len(targets)), joint,
        size=replicates, check_valid="raise", tol=1e-9,
    )
    prediction = draws[:, :len(observed)] @ weights.T
    error = prediction - draws[:, len(observed):]
    return float(np.mean(error ** 2))


def one_instance(seed: int, args):
    rng, edges, base, candidates, targets, covariance = build_instance(
        seed, args.nodes, args.base_size, args.candidates, args.targets
    )
    subsets = list(itertools.combinations(range(args.candidates), args.budget))
    subset_rows = []
    cache = {}
    for subset_id, local_subset in enumerate(subsets):
        selected = candidates[np.asarray(local_subset, dtype=int)]
        observed = np.r_[base, selected]
        risks, weights, _sigma = design_risks(
            edges, observed, targets, covariance, args.nodes, args.noise ** 2
        )
        cache[subset_id] = (observed, weights)
        subset_rows.append({
            "seed": seed, "subset_id": subset_id,
            "selected_local": ";".join(map(str, local_subset)), **risks,
        })
    table = pd.DataFrame(subset_rows)
    selected_ids = {
        "exact": int(table.exact.idxmin()),
        "pure_geometry": int(table.pure_geometry.idxmin()),
        "adequacy_only": int(table.adequacy_only.idxmin()),
    }
    optimum = float(table.exact.min())
    random_average = float(table.exact.mean())
    strategy_rows = []
    for strategy, subset_id in selected_ids.items():
        observed, weights = cache[subset_id]
        exact_risk = float(table.loc[subset_id, "exact"])
        empirical = empirical_mse(
            rng, edges, observed, targets, covariance, args.noise ** 2,
            weights, args.mc_replicates
        )
        strategy_rows.append({
            "seed": seed, "strategy": strategy, "subset_id": subset_id,
            "exact_target_risk": exact_risk,
            "empirical_target_mse": empirical,
            "mc_relative_error": abs(empirical - exact_risk) / exact_risk,
            "regret_to_global_optimum": exact_risk - optimum,
            "relative_regret": exact_risk / optimum - 1.0,
            "gain_vs_random_average": random_average - exact_risk,
        })
    return table, pd.DataFrame(strategy_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--instances", type=int, default=30)
    parser.add_argument("--mc-replicates", type=int, default=10000)
    parser.add_argument("--nodes", type=int, default=8)
    parser.add_argument("--base-size", type=int, default=9)
    parser.add_argument("--candidates", type=int, default=10)
    parser.add_argument("--targets", type=int, default=4)
    parser.add_argument("--budget", type=int, default=3)
    parser.add_argument("--noise", type=float, default=.20)
    parser.add_argument("--out", type=Path, default=Path("results/optimal_design"))
    args = parser.parse_args(); args.out.mkdir(parents=True, exist_ok=True)

    subset_tables, strategy_tables = [], []
    for seed in range(args.instances):
        subset, strategy = one_instance(seed, args)
        subset_tables.append(subset); strategy_tables.append(strategy)
    subsets = pd.concat(subset_tables, ignore_index=True)
    strategies = pd.concat(strategy_tables, ignore_index=True)
    subsets.to_csv(args.out / "all_subsets.csv", index=False)
    strategies.to_csv(args.out / "selected_designs.csv", index=False)
    summary = strategies.groupby("strategy").agg(
        instances=("seed", "nunique"), exact_risk_mean=("exact_target_risk", "mean"),
        empirical_mse_mean=("empirical_target_mse", "mean"),
        mc_relative_error_max=("mc_relative_error", "max"),
        relative_regret_mean=("relative_regret", "mean"),
        relative_regret_median=("relative_regret", "median"),
        strict_optimum_rate=("regret_to_global_optimum", lambda x: float(np.mean(np.asarray(x) <= 1e-12))),
    ).reset_index()
    summary.to_csv(args.out / "summary.csv", index=False)
    wide = strategies.pivot(index="seed", columns="strategy", values="exact_target_risk")
    improvements = []
    for comparator in ("pure_geometry", "adequacy_only"):
        relative_gain = (wide[comparator] - wide["exact"]) / wide[comparator]
        improvements.append({
            "comparison": f"exact_vs_{comparator}",
            "strict_win_rate": float((relative_gain > 1e-10).mean()),
            "noninferiority_rate": float((relative_gain >= -1e-10).mean()),
            "relative_gain_mean": float(relative_gain.mean()),
            "relative_gain_median": float(relative_gain.median()),
        })
    improvement_frame = pd.DataFrame(improvements)
    improvement_frame.to_csv(args.out / "comparisons.csv", index=False)
    exact_rows = strategies[strategies.strategy == "exact"]
    verdict = {
        "certificate": "exact criterion was exhaustively minimized over every candidate subset",
        "subsets_per_instance": int(len(list(itertools.combinations(range(args.candidates), args.budget)))),
        "instances": args.instances,
        "gate": "zero enumerated regret for exact design; Monte Carlo error <=5%; exact noninferior in every instance and strictly better than each reduced criterion in at least one instance",
        "maximum_exact_enumeration_regret": float(exact_rows.regret_to_global_optimum.max()),
        "maximum_exact_mc_relative_error": float(exact_rows.mc_relative_error.max()),
        "comparisons": improvements,
    }
    verdict["gate_pass"] = bool(
        verdict["maximum_exact_enumeration_regret"] <= 1e-12
        and verdict["maximum_exact_mc_relative_error"] <= .05
        and all(x["noninferiority_rate"] == 1.0 and x["strict_win_rate"] > 0 for x in improvements)
    )
    (args.out / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(summary.to_string(index=False))
    print("\n" + improvement_frame.to_string(index=False))
    print("\n" + json.dumps(verdict, indent=2))
    if not verdict["gate_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
