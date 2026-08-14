#!/usr/bin/env python3
"""Seed-level bootstrap confidence intervals for the formal Norman gate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def bootstrap_mean(values: np.ndarray, replicates: int, rng: np.random.Generator) -> tuple[float, float]:
    draw = rng.integers(0, len(values), size=(replicates, len(values)))
    means = values[draw].mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=Path("results/formal_30split"))
    parser.add_argument("--bootstrap", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20261207)
    args = parser.parse_args()
    metrics = pd.read_csv(args.results / "per_seed_metrics.csv")
    alignment = pd.read_csv(args.results / "risk_alignment.csv").set_index("seed")
    wide_mse = metrics.pivot(index="seed", columns="strategy", values="residual_mse")
    wide_cosine = metrics.pivot(index="seed", columns="strategy", values="residual_cosine")
    wide_full_cosine = metrics.pivot(index="seed", columns="strategy", values="full_effect_cosine")
    series = {
        "relative_residual_mse_improvement": (
            (wide_mse.geometry_only - wide_mse.estimated_witness)
            / wide_mse.geometry_only
        ).to_numpy(),
        "residual_cosine_gain": (
            wide_cosine.estimated_witness - wide_cosine.geometry_only
        ).to_numpy(),
        "full_effect_cosine_gain": (
            wide_full_cosine.estimated_witness - wide_full_cosine.geometry_only
        ).to_numpy(),
        "spearman_estimated_risk_vs_oracle_risk": alignment[
            "spearman_estimated_risk_vs_oracle_risk"
        ].to_numpy(),
        "spearman_estimated_risk_vs_realized_mse": alignment[
            "spearman_estimated_risk_vs_realized_mse"
        ].to_numpy(),
        "oracle_efficiency_estimated": alignment["oracle_efficiency_estimated"].to_numpy(),
        "oracle_efficiency_geometry": alignment["oracle_efficiency_geometry"].to_numpy(),
    }
    rng = np.random.default_rng(args.seed)
    rows = []
    for metric, values in series.items():
        low, high = bootstrap_mean(values, args.bootstrap, rng)
        rows.append({
            "metric": metric,
            "n_seeds": len(values),
            "mean": float(values.mean()),
            "sd": float(values.std(ddof=1)),
            "bootstrap_95ci_low": low,
            "bootstrap_95ci_high": high,
            "bootstrap_replicates": args.bootstrap,
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(args.results / "confidence_intervals.csv", index=False)
    summary = {
        "resampling_unit": "outer split seed",
        "bootstrap_replicates": args.bootstrap,
        "seed": args.seed,
        "interval": "percentile 95% interval of the split-level mean",
        "metrics": rows,
    }
    (args.results / "confidence_intervals.json").write_text(json.dumps(summary, indent=2))
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()

