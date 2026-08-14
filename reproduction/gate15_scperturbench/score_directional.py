#!/usr/bin/env python3
"""Cheap, exact-mean directional gate against published scPerturBench rows.

This is deliberately not the final six-metric leaderboard.  It evaluates the
two mean-response quantities underlying the published Pearson and MSE columns,
for top-100 and top-5000 condition-specific DE genes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


EPS = 1e-12


def pearson(left: np.ndarray, right: np.ndarray) -> float:
    left = left - left.mean()
    right = right - right.mean()
    return float(np.sum(left * right) / max(np.linalg.norm(left) * np.linalg.norm(right), EPS))


def welch_scores(
    truth_mean: np.ndarray,
    truth_variance: np.ndarray,
    truth_n: int,
    control_mean: np.ndarray,
    control_variance: np.ndarray,
    control_n: int,
) -> np.ndarray:
    truth_unbiased = truth_variance * truth_n / max(truth_n - 1, 1)
    control_unbiased = control_variance * control_n / max(control_n - 1, 1)
    denominator = np.sqrt(
        truth_unbiased / max(truth_n, 1)
        + control_unbiased / max(control_n, 1)
    )
    return (truth_mean - control_mean) / np.maximum(denominator, EPS)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction", type=Path, action="append", required=True)
    parser.add_argument("--published", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    published = pd.read_csv(args.published)
    rows = []
    for path in args.prediction:
        data = np.load(path, allow_pickle=False)
        dataset = str(data["dataset"])
        seed = int(data["seed"])
        control = data["control"].astype(float)
        control_variance = data["control_variance"].astype(float)
        control_n = int(data["control_count"])
        for index, condition in enumerate(data["conditions"].astype(str)):
            truth = data["truth"][index].astype(float)
            score = welch_scores(
                truth,
                data["truth_variance"][index].astype(float),
                int(data["test_cell_counts"][index]),
                control,
                control_variance,
                control_n,
            )
            # scPerturBench's calDEG explicitly re-sorts scanpy t-test scores
            # by absolute magnitude before taking the first N genes.
            order = np.argsort(-np.abs(score))
            for panel in (100, 5000):
                genes = order[: min(panel, len(order))]
                for method, key in (
                    ("WitnessCell", "prediction"),
                    ("factorized_only", "factorized_prediction"),
                ):
                    prediction = data[key][index].astype(float)
                    rows.append({
                        "DataSet": dataset,
                        "seed": seed,
                        "perturb": condition,
                        "subgroup": str(data["subgroups"][index]),
                        "DEG": panel,
                        "method": method,
                        "cor": pearson(
                            prediction[genes] - control[genes],
                            truth[genes] - control[genes],
                        ),
                        "mse": float(np.mean((prediction[genes] - truth[genes]) ** 2)),
                    })
    frame = pd.DataFrame(rows)
    frame.to_csv(args.out / "directional_metrics.csv", index=False)

    # Compare only against the raw published mean metrics.  The full benchmark
    # gate must still run the official distributional scorer.
    comparison_rows = []
    for (dataset, seed, panel), own in frame.groupby(["DataSet", "seed", "DEG"]):
        base = published[
            (published.DataSet == dataset)
            & (published.seed == seed)
            & (published.DEG == panel)
            & (published.metric.isin(["mse", "pearson_distance"]))
        ].copy()
        for method, current in own.groupby("method"):
            test = set(current.perturb)
            selected = base[base.perturb.isin(test)]
            pcc = selected[selected.metric == "pearson_distance"]
            mse = selected[selected.metric == "mse"]
            comparison_rows.append({
                "DataSet": dataset,
                "seed": int(seed),
                "DEG": int(panel),
                "method": method,
                "conditions": len(current),
                "mean_cor": float(current.cor.mean()),
                "published_best_mean_cor": float(pcc.groupby("method").performance.mean().max()),
                "mean_mse": float(current.mse.mean()),
                "published_best_mean_mse": float(mse.groupby("method").performance.mean().min()),
                "witness_vs_factorized_cor": float("nan"),
                "witness_vs_factorized_mse_gain": float("nan"),
            })
    comparison = pd.DataFrame(comparison_rows)
    for (dataset, seed, panel), block in comparison.groupby(["DataSet", "seed", "DEG"]):
        values = block.set_index("method")
        if {"WitnessCell", "factorized_only"}.issubset(values.index):
            mask = (
                (comparison.DataSet == dataset)
                & (comparison.seed == seed)
                & (comparison.DEG == panel)
            )
            comparison.loc[mask, "witness_vs_factorized_cor"] = (
                values.loc["WitnessCell", "mean_cor"]
                - values.loc["factorized_only", "mean_cor"]
            )
            comparison.loc[mask, "witness_vs_factorized_mse_gain"] = (
                values.loc["factorized_only", "mean_mse"]
                - values.loc["WitnessCell", "mean_mse"]
            ) / max(values.loc["factorized_only", "mean_mse"], EPS)
    comparison.to_csv(args.out / "directional_summary.csv", index=False)

    leaderboard_rows = []
    for (dataset, panel), own_panel in frame.groupby(["DataSet", "DEG"]):
        assembled = []
        for seed, own_seed in own_panel.groupby("seed"):
            baseline = published[
                (published.DataSet == dataset)
                & (published.seed == seed)
                & (published.DEG == panel)
                & (published.metric.isin(["mse", "pearson_distance"]))
                & (published.perturb.isin(set(own_seed.perturb)))
            ].pivot_table(
                index=["perturb", "method"],
                columns="metric",
                values="performance",
            ).reset_index().rename(columns={"pearson_distance": "cor"})
            # The ablation is diagnostic only and is not submitted as a 29th
            # benchmark method, so it must not alter official-method ranks.
            current = own_seed[
                own_seed.method == "WitnessCell"
            ][["perturb", "method", "cor", "mse"]]
            joined = pd.concat([baseline, current], ignore_index=True)
            joined["rank_cor"] = joined.groupby("perturb").cor.rank(
                ascending=False, method="min"
            )
            joined["rank_mse"] = joined.groupby("perturb").mse.rank(
                ascending=True, method="min"
            )
            joined["rank_two_metric"] = (
                joined.rank_cor + joined.rank_mse
            ) / 2.0
            assembled.append(joined.assign(seed=int(seed)))
        ranked = pd.concat(assembled, ignore_index=True)
        summary = ranked.groupby("method", as_index=False).agg(
            mean_rank_two_metric=("rank_two_metric", "mean"),
            mean_rank_cor=("rank_cor", "mean"),
            mean_rank_mse=("rank_mse", "mean"),
        ).sort_values(["mean_rank_two_metric", "method"])
        summary.insert(0, "DEG", int(panel))
        summary.insert(0, "DataSet", dataset)
        summary["rank_position_two_metric"] = np.arange(1, len(summary) + 1)
        leaderboard_rows.append(summary)
    leaderboard = pd.concat(leaderboard_rows, ignore_index=True)
    leaderboard.to_csv(args.out / "two_metric_leaderboard.csv", index=False)
    witness = comparison[comparison.method == "WitnessCell"]
    verdict = {
        "status": "DIRECTIONAL_GATE_ONLY",
        "not_a_full_sota_claim": True,
        "seeds": sorted(frame.seed.unique().astype(int).tolist()),
        "conditions": int(frame[frame.method == "WitnessCell"].perturb.nunique()),
        "mean_witness_vs_factorized_mse_gain": float(witness.witness_vs_factorized_mse_gain.mean()),
        "mean_witness_vs_factorized_cor_delta": float(witness.witness_vs_factorized_cor.mean()),
        "beats_published_best_mean_mse_cells": int((witness.mean_mse < witness.published_best_mean_mse).sum()),
        "beats_published_best_mean_cor_cells": int((witness.mean_cor > witness.published_best_mean_cor).sum()),
        "comparison_cells": int(len(witness)),
        "two_metric_rank_by_panel": {
            str(int(panel)): int(
                leaderboard[
                    (leaderboard.method == "WitnessCell")
                    & (leaderboard.DEG == panel)
                ].rank_position_two_metric.iloc[0]
            )
            for panel in sorted(leaderboard.DEG.unique())
        },
        "two_metric_mean_rank_by_panel": {
            str(int(panel)): float(
                leaderboard[
                    (leaderboard.method == "WitnessCell")
                    & (leaderboard.DEG == panel)
                ].mean_rank_two_metric.iloc[0]
            )
            for panel in sorted(leaderboard.DEG.unique())
        },
        "claim_boundary": "Recomputes scanpy-style Welch-t DEG panels from the official data; final SOTA still requires the official full distributional scorer and its frozen DEG artifact.",
    }
    (args.out / "verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(comparison.to_string(index=False))
    print(leaderboard.to_string(index=False))
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
