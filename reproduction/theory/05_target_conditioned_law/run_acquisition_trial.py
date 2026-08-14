#!/usr/bin/env python3
"""Outcome-blind acquisition trial for target-conditioned witness design.

Random, global D-optimal, and target-conditioned V-optimal policies receive
the same number of candidate double perturbations.  Candidate outcomes are
revealed only after selection.  D-optimal and V-optimal achieve the same global
rank/nullity improvement at a fixed budget; they differ only in whether the
new directions support the predeclared target set.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.linear_model import Ridge

from run_target_geometry_diagnosis import incidence, make_truth, row_cosine


EPS = 1e-12


def build_acquisition_design(n_components: int, side: int, target_components: set[int]):
    base, candidates, targets = [], [], []
    for component in range(n_components):
        offset = component * 2 * side
        left = np.arange(offset, offset + side)
        right = np.arange(offset + side, offset + 2 * side)
        held_cross = {(int(left[2]), int(right[2])), (int(left[3]), int(right[3]))}
        for i in left:
            for j in right:
                edge = (int(i), int(j))
                if edge not in held_cross:
                    base.append(edge)
        candidates.append((int(left[0]), int(left[1]), component))
        if component in target_components:
            targets.extend([
                (int(left[2]), int(left[3]), component),
                (int(right[2]), int(right[3]), component),
            ])
    return (
        np.asarray(base, dtype=int),
        pd.DataFrame(candidates, columns=["i", "j", "component"]),
        pd.DataFrame(targets, columns=["i", "j", "component"]),
    )


def greedy_order(base_edges: np.ndarray, candidate_edges: np.ndarray,
                 target_edges: np.ndarray, n_nodes: int, strategy: str,
                 ridge: float) -> np.ndarray:
    q = incidence(base_edges, n_nodes).T @ incidence(base_edges, n_nodes)
    q += ridge * np.eye(n_nodes)
    target_design = incidence(target_edges, n_nodes)
    remaining = list(range(len(candidate_edges)))
    order = []
    while remaining:
        inverse = np.linalg.inv(q)
        scores = []
        for position in remaining:
            x = incidence(candidate_edges[[position]], n_nodes)[0]
            denominator = 1.0 + float(x @ inverse @ x)
            if strategy == "d_optimal":
                score = float(x @ inverse @ x)
            elif strategy == "v_optimal":
                score = float(np.sum((target_design @ inverse @ x) ** 2) / denominator)
            else:
                raise ValueError(strategy)
            scores.append(score)
        local = int(np.argmax(scores))
        selected = remaining.pop(local)
        order.append(selected)
        x = incidence(candidate_edges[[selected]], n_nodes)[0]
        q += np.outer(x, x)
    return np.asarray(order, dtype=int)


def evaluate_selection(base_edges: np.ndarray, candidate_edges: np.ndarray,
                       selected: np.ndarray, target_edges: np.ndarray,
                       base_y: np.ndarray, candidate_y: np.ndarray,
                       target_truth: np.ndarray, n_nodes: int, fit_ridge: float):
    train_edges = np.vstack([base_edges, candidate_edges[selected]]) if len(selected) else base_edges
    train_y = np.vstack([base_y, candidate_y[selected]]) if len(selected) else base_y
    model = Ridge(alpha=fit_ridge, fit_intercept=False).fit(incidence(train_edges, n_nodes), train_y)
    prediction = model.predict(incidence(target_edges, n_nodes))
    cosine = float(np.mean(row_cosine(prediction, target_truth)))
    nmse = float(np.mean((prediction - target_truth) ** 2) / (np.mean(target_truth ** 2) + EPS))
    rank = int(np.linalg.matrix_rank(incidence(train_edges, n_nodes), tol=1e-9))
    return cosine, nmse, n_nodes - rank


def one_seed(seed: int, args) -> list[dict]:
    rng = np.random.default_rng(args.seed_base + seed)
    target_components = set(rng.choice(
        args.components, args.target_components, replace=False
    ).tolist())
    base, candidate_frame, target_frame = build_acquisition_design(
        args.components, args.side, target_components
    )
    candidate_frame = candidate_frame.iloc[rng.permutation(len(candidate_frame))].reset_index(drop=True)
    candidate_edges = candidate_frame[["i", "j"]].to_numpy(dtype=int)
    target_edges = target_frame[["i", "j"]].to_numpy(dtype=int)
    n_nodes = args.components * 2 * args.side
    theta = make_truth(rng, args.components, args.side, args.output_dim, args.contrast_scale)
    base_truth = theta[base[:, 0]] + theta[base[:, 1]]
    candidate_truth = theta[candidate_edges[:, 0]] + theta[candidate_edges[:, 1]]
    target_truth = theta[target_edges[:, 0]] + theta[target_edges[:, 1]]
    signal_rms = float(np.sqrt(np.mean(base_truth ** 2)))
    base_y = base_truth + rng.normal(scale=args.noise * signal_rms, size=base_truth.shape)
    candidate_y = candidate_truth + rng.normal(
        scale=args.noise * signal_rms, size=candidate_truth.shape
    )

    orders = {
        strategy: greedy_order(
            base, candidate_edges, target_edges, n_nodes, strategy, args.risk_ridge
        ) for strategy in ("d_optimal", "v_optimal")
    }
    rows = []
    for budget in args.budgets:
        if budget > len(candidate_edges):
            continue
        for strategy, order in orders.items():
            selected = order[:budget]
            cosine, nmse, nullity = evaluate_selection(
                base, candidate_edges, selected, target_edges, base_y, candidate_y,
                target_truth, n_nodes, args.fit_ridge
            )
            selected_components = set(candidate_frame.iloc[selected].component.astype(int))
            rows.append({
                "seed": seed, "strategy": strategy, "rep": 0, "budget": budget,
                "cosine": cosine, "nmse": nmse, "global_nullity": nullity,
                "target_components_witnessed": len(selected_components & target_components),
            })
        for rep in range(args.random_reps):
            random_rng = np.random.default_rng(args.random_seed_base + 1009 * seed + rep)
            selected = random_rng.choice(len(candidate_edges), size=budget, replace=False)
            cosine, nmse, nullity = evaluate_selection(
                base, candidate_edges, selected, target_edges, base_y, candidate_y,
                target_truth, n_nodes, args.fit_ridge
            )
            selected_components = set(candidate_frame.iloc[selected].component.astype(int))
            rows.append({
                "seed": seed, "strategy": "random", "rep": rep, "budget": budget,
                "cosine": cosine, "nmse": nmse, "global_nullity": nullity,
                "target_components_witnessed": len(selected_components & target_components),
            })
    return rows


def paired_test(frame: pd.DataFrame, budget: int, comparator: str) -> dict:
    vopt = frame[(frame.strategy == "v_optimal") & (frame.budget == budget)].set_index("seed")
    if comparator == "random":
        other = frame[(frame.strategy == comparator) & (frame.budget == budget)].groupby("seed").agg(
            cosine=("cosine", "mean"), nmse=("nmse", "mean"),
            global_nullity=("global_nullity", "mean"),
        )
    else:
        other = frame[(frame.strategy == comparator) & (frame.budget == budget)].set_index("seed")
    delta_cos = vopt.cosine - other.cosine
    delta_nmse = vopt.nmse - other.nmse
    if np.allclose(delta_cos, 0):
        p_value = 1.0
    else:
        p_value = float(wilcoxon(delta_cos, alternative="greater").pvalue)
    return {
        "budget": budget,
        "comparison": f"v_optimal_minus_{comparator}",
        "delta_cosine_mean": float(delta_cos.mean()),
        "delta_nmse_mean": float(delta_nmse.mean()),
        "cosine_win_rate": float((delta_cos > 0).mean()),
        "p_one_sided_wilcoxon": p_value,
        "mean_global_nullity_difference": float((vopt.global_nullity - other.global_nullity).mean()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=50)
    parser.add_argument("--seed-base", type=int, default=20260817)
    parser.add_argument("--random-seed-base", type=int, default=88001)
    parser.add_argument("--random-reps", type=int, default=20)
    parser.add_argument("--components", type=int, default=12)
    parser.add_argument("--target-components", type=int, default=4)
    parser.add_argument("--side", type=int, default=6)
    parser.add_argument("--output-dim", type=int, default=64)
    parser.add_argument("--contrast-scale", type=float, default=1.5)
    parser.add_argument("--noise", type=float, default=0.10)
    parser.add_argument("--fit-ridge", type=float, default=1e-5)
    parser.add_argument("--risk-ridge", type=float, default=1e-3)
    parser.add_argument("--budgets", type=int, nargs="+", default=[0, 1, 2, 4, 8, 12])
    parser.add_argument("--out", type=Path, default=Path("results/acquisition"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    frame = pd.DataFrame([row for seed in range(args.seeds) for row in one_seed(seed, args)])
    frame.to_csv(args.out / "acquisition_rows.csv", index=False)
    summary = frame.groupby(["strategy", "budget"]).agg(
        n=("cosine", "size"), cosine_mean=("cosine", "mean"), cosine_sd=("cosine", "std"),
        nmse_mean=("nmse", "mean"), nullity_mean=("global_nullity", "mean"),
        target_witnesses_mean=("target_components_witnessed", "mean"),
    ).reset_index()
    summary.to_csv(args.out / "acquisition_summary.csv", index=False)

    tests = [
        paired_test(frame, budget, comparator)
        for budget in args.budgets if budget > 0 and budget <= args.components
        for comparator in ("random", "d_optimal")
    ]
    pd.DataFrame(tests).to_csv(args.out / "paired_tests.csv", index=False)
    primary_budget = args.target_components
    primary = {x["comparison"]: x for x in tests if x["budget"] == primary_budget}
    random_test = primary["v_optimal_minus_random"]
    dopt_test = primary["v_optimal_minus_d_optimal"]
    audit = {
        "candidate_outcomes_used_for_selection": False,
        "target_outcomes_used_for_selection": False,
        "target_identities_used_by_v_optimal": True,
        "all_targets_seen2": True,
        "same_budget_all_strategies": True,
        "primary_budget": primary_budget,
    }
    verdict = {
        "claim": "At matched budget and matched global nullity improvement, target-conditioned geometry yields higher target prediction utility than random or global D-optimal acquisition.",
        "audit": audit,
        "primary_vopt_minus_random": random_test,
        "primary_vopt_minus_dopt": dopt_test,
        "pass_rule": "at the primary budget, V-optimal must exceed both comparators by >=0.15 cosine, win >=90% of seeds, p<0.01, and have zero mean global-nullity difference",
    }
    verdict["gate_pass"] = bool(all(
        test["delta_cosine_mean"] >= 0.15
        and test["cosine_win_rate"] >= 0.90
        and test["p_one_sided_wilcoxon"] < 0.01
        and abs(test["mean_global_nullity_difference"]) < 1e-9
        for test in (random_test, dopt_test)
    ))
    (args.out / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(summary.to_string(index=False))
    print("\n" + pd.DataFrame(tests).to_string(index=False))
    print("\n" + json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()

