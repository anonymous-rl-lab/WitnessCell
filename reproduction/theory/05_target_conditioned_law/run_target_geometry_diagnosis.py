#!/usr/bin/env python3
"""Held-out target test for target-conditioned intervention geometry.

Every target has both endpoints represented in the training graph (seen2).
The number of training edges and the global nullity are fixed within a seed.
Only the target direction relative to the observed unsigned-incidence row
space changes.  This separates target estimability from exposure and global
identifiability.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score


EPS = 1e-12


def incidence(edges: np.ndarray, n_nodes: int) -> np.ndarray:
    matrix = np.zeros((len(edges), n_nodes), dtype=float)
    matrix[np.arange(len(edges)), edges[:, 0]] = 1.0
    matrix[np.arange(len(edges)), edges[:, 1]] = 1.0
    return matrix


def row_cosine(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.sum(a * b, axis=1) / (
        np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + EPS
    )


def build_design(n_components: int, side: int, witnessed_components: set[int]):
    """Construct equal-size bipartite components and strictly held-out targets."""
    train, targets, candidate_witnesses = [], [], []
    for component in range(n_components):
        offset = component * 2 * side
        left = np.arange(offset, offset + side)
        right = np.arange(offset + side, offset + 2 * side)
        held_cross = {(int(left[2]), int(right[2])), (int(left[3]), int(right[3]))}
        for i in left:
            for j in right:
                edge = (int(i), int(j))
                if edge not in held_cross:
                    train.append(edge)
        witness = (int(left[0]), int(left[1]))
        candidate_witnesses.append(witness)
        if component in witnessed_components:
            train.append(witness)
        status = "witnessed" if component in witnessed_components else "unwitnessed"
        target_specs = [
            (int(left[2]), int(right[2]), "cross"),
            (int(left[3]), int(right[3]), "cross"),
            (int(left[2]), int(left[3]), "within"),
            (int(right[2]), int(right[3]), "within"),
        ]
        for i, j, pair_type in target_specs:
            targets.append({
                "i": i, "j": j, "component": component,
                "pair_type": pair_type, "component_status": status,
                "unsupported_truth": int(status == "unwitnessed" and pair_type == "within"),
            })
    return np.asarray(train, dtype=int), pd.DataFrame(targets), np.asarray(candidate_witnesses, dtype=int)


def make_truth(rng: np.random.Generator, n_components: int, side: int,
               output_dim: int, contrast_scale: float) -> np.ndarray:
    n_nodes = n_components * 2 * side
    theta = rng.normal(scale=0.45, size=(n_nodes, output_dim))
    for component in range(n_components):
        offset = component * 2 * side
        contrast = rng.normal(scale=contrast_scale, size=output_dim)
        theta[offset:offset + side] += contrast
        theta[offset + side:offset + 2 * side] -= contrast
    return theta


def structural_scores(train_edges: np.ndarray, target_edges: np.ndarray,
                      n_nodes: int, ridge: float):
    design = incidence(train_edges, n_nodes)
    q = design.T @ design
    eigvals, eigvecs = np.linalg.eigh(q)
    null = eigvecs[:, eigvals <= 1e-9]
    inverse = np.linalg.inv(q + ridge * np.eye(n_nodes))
    target_design = incidence(target_edges, n_nodes)
    null_projection = np.sum((target_design @ null) ** 2, axis=1) / 2.0
    regularized_variance = np.einsum("ij,jk,ik->i", target_design, inverse, target_design)
    return null_projection, regularized_variance, int(null.shape[1]), float(eigvals.max())


def safe_auc(y: np.ndarray, score: np.ndarray) -> float:
    return float(roc_auc_score(y, score)) if len(np.unique(y)) == 2 else float("nan")


def one_seed(seed: int, args) -> tuple[pd.DataFrame, dict]:
    rng = np.random.default_rng(args.seed_base + seed)
    witnessed = set(rng.choice(args.components, args.components // 2, replace=False).tolist())
    train, target_frame, _ = build_design(args.components, args.side, witnessed)
    n_nodes = args.components * 2 * args.side
    target_edges = target_frame[["i", "j"]].to_numpy(dtype=int)
    theta = make_truth(rng, args.components, args.side, args.output_dim, args.contrast_scale)
    train_truth = theta[train[:, 0]] + theta[train[:, 1]]
    signal_rms = float(np.sqrt(np.mean(train_truth ** 2)))
    train_y = train_truth + rng.normal(scale=args.noise * signal_rms, size=train_truth.shape)
    target_truth = theta[target_edges[:, 0]] + theta[target_edges[:, 1]]

    model = Ridge(alpha=args.fit_ridge, fit_intercept=False).fit(
        incidence(train, n_nodes), train_y
    )
    prediction = model.predict(incidence(target_edges, n_nodes))
    null_risk, variance_risk, nullity, max_eigenvalue = structural_scores(
        train, target_edges, n_nodes, args.risk_ridge
    )
    degree = np.bincount(train.ravel(), minlength=n_nodes)
    nmse = np.mean((prediction - target_truth) ** 2, axis=1) / (
        np.mean(target_truth ** 2, axis=1) + EPS
    )
    target_frame = target_frame.assign(
        seed=seed,
        seen_category="seen2",
        train_edges=len(train),
        global_nullity=nullity,
        endpoint_degree_sum=degree[target_edges[:, 0]] + degree[target_edges[:, 1]],
        null_projection=null_risk,
        regularized_variance=variance_risk,
        cosine=row_cosine(prediction, target_truth),
        nmse=nmse,
        catastrophic=(nmse >= args.failure_nmse).astype(int),
    )
    overlap = set(map(tuple, train)) & set(map(tuple, target_edges))
    audit = {
        "seed": seed,
        "train_target_overlap": len(overlap),
        "all_target_endpoints_seen": bool(np.all(degree[np.unique(target_edges)] > 0)),
        "train_edges": len(train),
        "global_nullity": nullity,
        "max_information_eigenvalue": max_eigenvalue,
    }
    return target_frame, audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=100)
    parser.add_argument("--seed-base", type=int, default=20260807)
    parser.add_argument("--components", type=int, default=8)
    parser.add_argument("--side", type=int, default=7)
    parser.add_argument("--output-dim", type=int, default=64)
    parser.add_argument("--contrast-scale", type=float, default=1.5)
    parser.add_argument("--noise", type=float, default=0.10)
    parser.add_argument("--fit-ridge", type=float, default=1e-5)
    parser.add_argument("--risk-ridge", type=float, default=1e-3)
    parser.add_argument("--failure-nmse", type=float, default=0.25)
    parser.add_argument("--out", type=Path, default=Path("results/diagnosis"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    frames, audits = [], []
    for seed in range(args.seeds):
        frame, audit = one_seed(seed, args)
        frames.append(frame); audits.append(audit)
    rows = pd.concat(frames, ignore_index=True)
    rows.to_csv(args.out / "target_rows.csv", index=False)

    group = rows.groupby(["component_status", "pair_type"]).agg(
        n=("nmse", "size"), cosine_mean=("cosine", "mean"), cosine_sd=("cosine", "std"),
        nmse_mean=("nmse", "mean"), catastrophic_rate=("catastrophic", "mean"),
        null_risk_mean=("null_projection", "mean"), degree_sum_mean=("endpoint_degree_sum", "mean"),
    ).reset_index()
    group.to_csv(args.out / "group_summary.csv", index=False)

    y_error = rows["nmse"].to_numpy()
    y_fail = rows["catastrophic"].to_numpy()
    metrics = {
        "null_projection": rows["null_projection"].to_numpy(),
        "regularized_variance": rows["regularized_variance"].to_numpy(),
        "within_pair_indicator": (rows["pair_type"] == "within").astype(float).to_numpy(),
        "inverse_endpoint_degree": 1.0 / rows["endpoint_degree_sum"].to_numpy(),
        "seen_category": np.full(len(rows), 2.0),
        "global_nullity": rows["global_nullity"].to_numpy(dtype=float),
        "train_edge_count": rows["train_edges"].to_numpy(dtype=float),
    }
    metric_rows = []
    for name, score in metrics.items():
        rho = spearmanr(score, y_error).statistic if np.ptp(score) > 1e-10 else np.nan
        metric_rows.append({
            "metric": name,
            "spearman_with_nmse": float(rho) if np.isfinite(rho) else np.nan,
            "failure_auroc": safe_auc(y_fail, score) if np.ptp(score) > 1e-10 else 0.5,
        })
    metric_frame = pd.DataFrame(metric_rows)
    metric_frame.to_csv(args.out / "metric_comparison.csv", index=False)

    audit_ok = (
        all(a["train_target_overlap"] == 0 for a in audits)
        and all(a["all_target_endpoints_seen"] for a in audits)
        and len({a["train_edges"] for a in audits}) == 1
        and len({a["global_nullity"] for a in audits}) == 1
        and rows["seen_category"].nunique() == 1
    )
    geometry_auc = float(metric_frame.loc[
        metric_frame.metric == "null_projection", "failure_auroc"
    ].iloc[0])
    baseline_auc = float(metric_frame.loc[
        metric_frame.metric == "within_pair_indicator", "failure_auroc"
    ].iloc[0])
    unwitnessed_within = group[
        (group.component_status == "unwitnessed") & (group.pair_type == "within")
    ].iloc[0]
    other = rows[~((rows.component_status == "unwitnessed") & (rows.pair_type == "within"))]
    verdict = {
        "claim": "Target estimability is governed by target direction relative to the observed intervention row space, not by seen2 exposure or global nullity alone.",
        "audit_pass": bool(audit_ok),
        "all_targets_strictly_held_out": True,
        "all_targets_seen2": True,
        "fixed_train_edges": int(rows.train_edges.iloc[0]),
        "fixed_global_nullity": int(rows.global_nullity.iloc[0]),
        "geometry_failure_auroc": geometry_auc,
        "within_only_baseline_auroc": baseline_auc,
        "unsupported_within_nmse": float(unwitnessed_within.nmse_mean),
        "other_targets_nmse": float(other.nmse.mean()),
        "pass_rule": "audit_pass and geometry AUROC >= 0.90 and geometry exceeds within-only baseline by >= 0.10 and unsupported NMSE >= 3x other NMSE",
    }
    verdict["gate_pass"] = bool(
        audit_ok and geometry_auc >= 0.90 and geometry_auc - baseline_auc >= 0.10
        and verdict["unsupported_within_nmse"] >= 3.0 * verdict["other_targets_nmse"]
    )
    (args.out / "audit.json").write_text(json.dumps(audits, indent=2))
    (args.out / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(group.to_string(index=False))
    print("\n" + metric_frame.to_string(index=False))
    print("\n" + json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
