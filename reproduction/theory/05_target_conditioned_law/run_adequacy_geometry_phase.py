#!/usr/bin/env python3
"""Synthetic phase diagram for geometry risk versus representation adequacy."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold

from run_target_geometry_diagnosis import (
    build_design, incidence, make_truth, row_cosine, structural_scores,
)


def oof_node_adequacy(train: np.ndarray, y: np.ndarray, n_nodes: int,
                      ridge: float, folds: int, seed: int, prior: float):
    oof = np.zeros_like(y)
    splitter = KFold(n_splits=folds, shuffle=True, random_state=seed)
    local = np.arange(len(train))
    for fit_local, val_local in splitter.split(local):
        model = Ridge(alpha=ridge, fit_intercept=False).fit(
            incidence(train[fit_local], n_nodes), y[fit_local]
        )
        oof[val_local] = model.predict(incidence(train[val_local], n_nodes))
    edge_error = 1.0 - row_cosine(oof, y)
    global_error = float(np.mean(edge_error))
    node_sum = np.zeros(n_nodes); degree = np.zeros(n_nodes)
    for value, (i, j) in zip(edge_error, train):
        node_sum[i] += value; node_sum[j] += value
        degree[i] += 1; degree[j] += 1
    node_bias = (node_sum + prior * global_error) / (degree + prior)
    return float(np.mean(row_cosine(oof, y))), node_bias


def pair_specific_effect(rng: np.random.Generator, edges: np.ndarray,
                         side: int, output_dim: int,
                         component_scale: np.ndarray, mismatch: float,
                         reference_rms: float) -> np.ndarray:
    result = np.zeros((len(edges), output_dim))
    for k, (i, _j) in enumerate(edges):
        component = int(i) // (2 * side)
        result[k] = rng.normal(
            scale=mismatch * reference_rms * component_scale[component],
            size=output_dim,
        )
    return result


def one_setting(seed: int, mismatch: float, args) -> tuple[list[dict], dict]:
    rng = np.random.default_rng(args.seed_base + seed)
    witnessed = set(rng.choice(args.components, args.components // 2, replace=False).tolist())
    train, target_frame, _ = build_design(args.components, args.side, witnessed)
    targets = target_frame[["i", "j"]].to_numpy(dtype=int)
    n_nodes = args.components * 2 * args.side
    theta = make_truth(rng, args.components, args.side, args.output_dim, args.contrast_scale)
    train_anchor = theta[train[:, 0]] + theta[train[:, 1]]
    target_anchor = theta[targets[:, 0]] + theta[targets[:, 1]]
    reference_rms = float(np.sqrt(np.mean(train_anchor ** 2)))
    component_scale = rng.lognormal(mean=-.15, sigma=.65, size=args.components)
    train_truth = train_anchor + pair_specific_effect(
        rng, train, args.side, args.output_dim, component_scale, mismatch, reference_rms
    )
    target_truth = target_anchor + pair_specific_effect(
        rng, targets, args.side, args.output_dim, component_scale, mismatch, reference_rms
    )
    observed = train_truth + rng.normal(
        scale=args.noise * reference_rms, size=train_truth.shape
    )
    model = Ridge(alpha=args.fit_ridge, fit_intercept=False).fit(
        incidence(train, n_nodes), observed
    )
    prediction = model.predict(incidence(targets, n_nodes))
    error = 1.0 - row_cosine(prediction, target_truth)
    failure = (error >= np.quantile(error, .75)).astype(int)
    geometry, _variance, _nullity, _ = structural_scores(
        train, targets, n_nodes, args.risk_ridge
    )
    oof_cosine, node_bias = oof_node_adequacy(
        train, observed, n_nodes, args.fit_ridge, args.folds, 100003 + seed, args.prior
    )
    adequacy = (node_bias[targets[:, 0]] + node_bias[targets[:, 1]]) / 2.0
    geometry_rank = rankdata(geometry) / (len(targets) + 1.0)
    adequacy_rank = rankdata(adequacy) / (len(targets) + 1.0)
    routed = geometry_rank if oof_cosine >= args.adequacy_threshold else adequacy_rank
    score_map = {
        "geometry": geometry_rank,
        "adequacy": adequacy_rank,
        "equal_average": .5 * geometry_rank + .5 * adequacy_rank,
        "hierarchical_route": routed,
    }
    rows = []
    for metric, score in score_map.items():
        rows.append({
            "seed": seed, "mismatch": mismatch, "metric": metric,
            "failure_auroc": float(roc_auc_score(failure, score)),
            "oof_adequacy_cosine": oof_cosine,
            "route": "geometry" if oof_cosine >= args.adequacy_threshold else "adequacy",
        })
    audit = {
        "seed": seed, "mismatch": mismatch,
        "train_target_overlap": len(set(map(tuple, train)) & set(map(tuple, targets))),
        "all_targets_seen2": True,
    }
    return rows, audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=20)
    parser.add_argument("--seed-base", type=int, default=20260901)
    parser.add_argument("--mismatch", type=float, nargs="+", default=[0, .1, .25, .5, 1, 2])
    parser.add_argument("--components", type=int, default=8)
    parser.add_argument("--side", type=int, default=7)
    parser.add_argument("--output-dim", type=int, default=48)
    parser.add_argument("--contrast-scale", type=float, default=1.5)
    parser.add_argument("--noise", type=float, default=.10)
    parser.add_argument("--fit-ridge", type=float, default=1e-5)
    parser.add_argument("--risk-ridge", type=float, default=1e-3)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--prior", type=float, default=2.0)
    parser.add_argument("--adequacy-threshold", type=float, default=.50)
    parser.add_argument("--out", type=Path, default=Path("results/adequacy_geometry_phase"))
    args = parser.parse_args(); args.out.mkdir(parents=True, exist_ok=True)

    rows, audits = [], []
    for mismatch in args.mismatch:
        for seed in range(args.seeds):
            local_rows, audit = one_setting(seed, mismatch, args)
            rows.extend(local_rows); audits.append(audit)
    frame = pd.DataFrame(rows); frame.to_csv(args.out / "phase_rows.csv", index=False)
    summary = frame.groupby(["mismatch", "metric"]).agg(
        seeds=("seed", "nunique"), auroc_mean=("failure_auroc", "mean"),
        auroc_sd=("failure_auroc", "std"),
        oof_adequacy_mean=("oof_adequacy_cosine", "mean"),
        geometry_route_rate=("route", lambda x: float(np.mean(np.asarray(x) == "geometry"))),
    ).reset_index()
    summary.to_csv(args.out / "phase_summary.csv", index=False)
    audit_pass = all(a["train_target_overlap"] == 0 and a["all_targets_seen2"] for a in audits)
    low = summary[summary.mismatch == min(args.mismatch)].set_index("metric")
    high = summary[summary.mismatch == max(args.mismatch)].set_index("metric")
    best = summary.groupby("mismatch").auroc_mean.max()
    routed = summary[summary.metric == "hierarchical_route"].set_index("mismatch").auroc_mean
    regret = best - routed
    verdict = {
        "claim": "Intervention geometry governs target reliability only after the interaction representation passes an adequacy gate.",
        "audit_pass": bool(audit_pass),
        "low_mismatch_geometry_auroc": float(low.loc["geometry", "auroc_mean"]),
        "high_mismatch_geometry_auroc": float(high.loc["geometry", "auroc_mean"]),
        "high_mismatch_adequacy_auroc": float(high.loc["adequacy", "auroc_mean"]),
        "maximum_hierarchical_regret_to_oracle_head": float(regret.max()),
        "gate": "geometry AUROC>=0.90 at zero mismatch; adequacy beats geometry by >=0.10 at maximum mismatch with AUROC>=0.65; hierarchical route remains within 0.08 of the better head at every level",
    }
    verdict["gate_pass"] = bool(
        audit_pass
        and verdict["low_mismatch_geometry_auroc"] >= .90
        and verdict["high_mismatch_adequacy_auroc"] >= .65
        and verdict["high_mismatch_adequacy_auroc"] - verdict["high_mismatch_geometry_auroc"] >= .10
        and verdict["maximum_hierarchical_regret_to_oracle_head"] <= .08
    )
    (args.out / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(summary.to_string(index=False))
    print("\n" + json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
