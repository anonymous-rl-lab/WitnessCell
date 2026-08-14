#!/usr/bin/env python3
"""Four-metric CPU pilot using the old scPerturBench/pertpy definitions.

This covers Pearson, MSE, Euclidean energy distance, and symmetric Gaussian KL.
Wasserstein and the final aggregate transformation remain reserved for the
official full scorer, so the output cannot be used as a formal SOTA claim.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from score_directional import pearson, welch_scores


EPS = 1e-8


def dense(value):
    return value.toarray() if hasattr(value, "toarray") else np.asarray(value)


def edistance(left: np.ndarray, right: np.ndarray) -> float:
    within_left = float(cdist(left, left, metric="euclidean").mean())
    within_right = float(cdist(right, right, metric="euclidean").mean())
    between = float(cdist(left, right, metric="euclidean").mean())
    return 2.0 * between - within_left - within_right


def symmetric_kl(left: np.ndarray, right: np.ndarray) -> float:
    left_mean = left.mean(axis=0)
    right_mean = right.mean(axis=0)
    left_std = left.std(axis=0) + EPS
    right_std = right.std(axis=0) + EPS
    difference = np.square(left_mean - right_mean)
    forward = (
        np.log(right_std / left_std)
        + (np.square(left_std) + difference) / (2.0 * np.square(right_std))
        - 0.5
    )
    reverse = (
        np.log(left_std / right_std)
        + (np.square(right_std) + difference) / (2.0 * np.square(left_std))
        - 0.5
    )
    return float(np.log2(1.0 + np.mean(forward + reverse)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--prediction", type=Path, action="append", required=True)
    parser.add_argument("--published", type=Path, required=True)
    parser.add_argument("--panel", type=int, default=100)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    source = ad.read_h5ad(args.data, backed="r")
    labels = source.obs.perturbation.astype(str).to_numpy()
    published = pd.read_csv(args.published)
    rows = []
    for path in args.prediction:
        archive = np.load(path, allow_pickle=False)
        seed = int(archive["seed"])
        dataset = str(archive["dataset"])
        control = archive["control"].astype(float)
        control_var = archive["control_variance"].astype(float)
        control_n = int(archive["control_count"])
        std = np.sqrt(np.maximum(control_var, 0.0))
        rng = np.random.default_rng(3000017 + seed)
        for index, condition in enumerate(archive["conditions"].astype(str)):
            truth_mean = archive["truth"][index].astype(float)
            score = welch_scores(
                truth_mean,
                archive["truth_variance"][index].astype(float),
                int(archive["test_cell_counts"][index]),
                control,
                control_var,
                control_n,
            )
            genes = np.argsort(-np.abs(score))[: args.panel]
            observed = dense(source.X[np.flatnonzero(labels == condition)])[:, genes].astype(float)
            predicted = rng.normal(
                loc=archive["prediction"][index, genes].astype(float),
                scale=std[genes],
                size=(len(observed), len(genes)),
            )
            rows.append({
                "DataSet": dataset,
                "seed": seed,
                "perturb": condition,
                "DEG": args.panel,
                "method": "WitnessCell",
                "pearson_distance": pearson(
                    predicted.mean(axis=0) - control[genes],
                    observed.mean(axis=0) - control[genes],
                ),
                "mse": float(np.mean(np.square(
                    predicted.mean(axis=0) - observed.mean(axis=0)
                ))),
                "edistance": edistance(predicted, observed),
                "sym_kldiv": symmetric_kl(predicted, observed),
            })
    frame = pd.DataFrame(rows)
    frame.to_csv(args.out / "four_metric_raw.csv", index=False)

    metrics = ["pearson_distance", "mse", "edistance", "sym_kldiv"]
    assembled = []
    for seed, own in frame.groupby("seed"):
        baseline = published[
            (published.DataSet == own.DataSet.iloc[0])
            & (published.seed == seed)
            & (published.DEG == args.panel)
            & (published.metric.isin(metrics))
            & (published.perturb.isin(set(own.perturb)))
        ].pivot_table(
            index=["perturb", "method"],
            columns="metric",
            values="performance",
        ).reset_index()
        joined = pd.concat([
            baseline,
            own[["perturb", "method", *metrics]],
        ], ignore_index=True)
        joined["rank_cor"] = joined.groupby("perturb").pearson_distance.rank(
            ascending=False, method="min"
        )
        for metric in ("mse", "edistance", "sym_kldiv"):
            joined[f"rank_{metric}"] = joined.groupby("perturb")[metric].rank(
                ascending=True, method="min"
            )
        rank_columns = ["rank_cor", "rank_mse", "rank_edistance", "rank_sym_kldiv"]
        joined["rank_four_metric"] = joined[rank_columns].mean(axis=1)
        assembled.append(joined.assign(seed=int(seed)))
    combined = pd.concat(assembled, ignore_index=True)
    leaderboard = combined.groupby("method", as_index=False).agg(
        mean_rank_four_metric=("rank_four_metric", "mean"),
        mean_rank_cor=("rank_cor", "mean"),
        mean_rank_mse=("rank_mse", "mean"),
        mean_rank_edistance=("rank_edistance", "mean"),
        mean_rank_sym_kldiv=("rank_sym_kldiv", "mean"),
    ).sort_values(["mean_rank_four_metric", "method"])
    leaderboard["position"] = np.arange(1, len(leaderboard) + 1)
    leaderboard.to_csv(args.out / "four_metric_leaderboard.csv", index=False)
    witness = leaderboard[leaderboard.method == "WitnessCell"].iloc[0]
    verdict = {
        "status": "FOUR_METRIC_CPU_PILOT",
        "not_a_formal_sota_claim": True,
        "panel": args.panel,
        "conditions_scored": len(frame),
        "witness_position": int(witness.position),
        "witness_mean_rank": float(witness.mean_rank_four_metric),
        "missing_formal_components": [
            "official frozen DEG_hvg5000.pkl identity check",
            "official OTT-JAX Wasserstein",
            "DEG recovery score",
            "published aggregate-score transformation",
        ],
    }
    (args.out / "verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(leaderboard.to_string(index=False))
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
