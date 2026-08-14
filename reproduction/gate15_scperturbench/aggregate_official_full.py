#!/usr/bin/env python3
"""Insert WitnessCell raw metrics into the published combo leaderboard."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


NAME_MAP = {"bioLord": "biolord"}


def names(series: pd.Series) -> pd.Series:
    return series.replace(NAME_MAP)


def inverse_minmax(values: pd.Series) -> pd.Series:
    span = float(values.max() - values.min())
    if span == 0.0:
        return pd.Series(np.ones(len(values)), index=values.index)
    return (values.max() - values) / span


def assemble_panel(
    panel: int,
    published_raw: pd.DataFrame,
    published_summary: pd.DataFrame,
    witness_raw: pd.DataFrame | None,
) -> pd.DataFrame:
    metrics = ["pearson_distance", "mse", "edistance", "sym_kldiv"]
    if panel == 100:
        metrics += ["wasserstein"]
    base = published_raw[
        (published_raw.DEG == panel) & published_raw.metric.isin(metrics)
    ].copy()
    base["method"] = names(base.method)
    base["op"] = base.perturb.astype(str) + "_" + base.seed.astype(str)
    base = base.pivot_table(
        index=["DataSet", "op", "method"],
        columns="metric",
        values="performance",
    ).reset_index()
    # The downloadable published raw table floors negative Pearson values at
    # zero, whereas the released aggregate table retains the signed values.
    # Consequently the official aggregate ``cor`` column is the authoritative
    # baseline source for PCC-delta.  All other baseline metrics remain rebuilt
    # from the raw table below.
    official_cor = published_summary[["DataSet", "op", "method", "cor"]].copy()
    official_cor["method"] = names(official_cor.method)
    base = base.merge(
        official_cor,
        on=["DataSet", "op", "method"],
        validate="one_to_one",
    ).drop(columns="pearson_distance")
    if panel == 100:
        deg = published_summary[["DataSet", "op", "method", "DEG_score"]].copy()
        deg["method"] = names(deg.method)
        base = base.merge(deg, on=["DataSet", "op", "method"], validate="one_to_one")
        base = base.rename(columns={"DEG_score": "common_degs"})

    if witness_raw is not None:
        own = witness_raw[
            (witness_raw.DEG == panel)
            & witness_raw.metric.isin(metrics + (["common_degs"] if panel == 100 else []))
        ].copy()
        own["op"] = own.perturb.astype(str) + "_" + own.seed.astype(str)
        own = own.pivot_table(
            index=["DataSet", "op", "method"],
            columns="metric",
            values="performance",
        ).reset_index()
        own = own.rename(columns={"pearson_distance": "cor"})
        base = pd.concat([base, own], ignore_index=True)

    blocks = []
    for (dataset, op), block in base.groupby(["DataSet", "op"], sort=False):
        block = block.copy()
        block["edistance_score"] = inverse_minmax(block["edistance"])
        block["mse_score"] = inverse_minmax(block["mse"])
        block["sym_score"] = inverse_minmax(block["sym_kldiv"])
        score_columns = ["cor", "edistance_score", "mse_score", "sym_score"]
        rank_columns = []
        block["Rank_pcc"] = block.cor.rank(ascending=False, method="min")
        block["Rank_edistance"] = block.edistance.rank(ascending=True, method="min")
        block["Rank_mse"] = block.mse.rank(ascending=True, method="min")
        block["Rank_sym"] = block.sym_kldiv.rank(ascending=True, method="min")
        rank_columns.extend(["Rank_pcc", "Rank_edistance", "Rank_mse", "Rank_sym"])
        if panel == 100:
            block["DEG_score"] = block["common_degs"]
            block["was_score"] = inverse_minmax(block["wasserstein"])
            block["Rank_deg"] = block.common_degs.rank(ascending=False, method="min")
            block["Rank_was"] = block.wasserstein.rank(ascending=True, method="min")
            score_columns.extend(["DEG_score", "was_score"])
            rank_columns.extend(["Rank_deg", "Rank_was"])
        block["ac_score"] = block[score_columns].mean(axis=1)
        block["ac_rank"] = block[rank_columns].mean(axis=1)
        block["Rank"] = block.ac_rank.rank(ascending=True, method="min")
        block["DEG"] = panel
        blocks.append(block)
    return pd.concat(blocks, ignore_index=True)


def validate_reconstruction(rebuilt: pd.DataFrame, official: pd.DataFrame, panel: int) -> dict:
    columns = [
        "cor", "edistance_score", "mse_score", "sym_score", "ac_score",
        "ac_rank", "Rank_pcc", "Rank_edistance", "Rank_mse", "Rank_sym", "Rank",
    ]
    if panel == 100:
        columns += ["DEG_score", "was_score", "Rank_deg", "Rank_was"]
    left = rebuilt.copy()
    right = official.copy()
    left["method"] = names(left.method)
    right["method"] = names(right.method)
    joined = left.merge(
        right[["DataSet", "op", "method", *columns]],
        on=["DataSet", "op", "method"],
        suffixes=("_new", "_official"),
        validate="one_to_one",
    )
    errors = {}
    for column in columns:
        errors[column] = float(np.nanmax(np.abs(
            joined[f"{column}_new"].astype(float)
            - joined[f"{column}_official"].astype(float)
        )))
    rank_columns = [column for column in columns if column.startswith("Rank") or column == "ac_rank"]
    exact_ranks = all(errors[column] < 1e-9 for column in rank_columns)
    # Published display scores are rounded, while min-max scores reconstructed
    # from four-decimal raw distances retain more precision.
    score_close = all(errors[column] < 5e-4 for column in columns if column not in rank_columns)
    return {
        "panel": panel,
        "rows": len(joined),
        "max_abs_errors": errors,
        "exact_published_ranks": exact_ranks,
        "published_scores_within_rounding_tolerance": score_close,
        "status": "PASS" if exact_ranks and score_close else "FAIL",
    }


def write_metric_leaderboard(full: pd.DataFrame, out: Path) -> None:
    columns = {
        "Rank_pcc": "PCC-delta",
        "Rank_mse": "MSE",
        "Rank_edistance": "E-distance",
        "Rank_sym": "symmetric KL",
        "Rank_was": "Wasserstein",
        "Rank_deg": "Common-DEGs",
    }
    top = full[full.DEG == 100]
    rows = []
    for column, label in columns.items():
        values = top.groupby("method")[column].mean().sort_values()
        for position, (method, value) in enumerate(values.items(), start=1):
            rows.append({
                "metric": label,
                "method": method,
                "mean_metric_rank": float(value),
                "position": position,
            })
    pd.DataFrame(rows).to_csv(out / "official_top100_per_metric_leaderboard.csv", index=False)


def write_cluster_bootstrap(
    full: pd.DataFrame,
    out: Path,
    replicates: int = 20000,
    seed: int = 20260809,
) -> None:
    top = full[full.DEG == 100].copy()
    top["perturb"] = top.op.str.rsplit("_", n=1).str[0]
    wide = top.pivot(
        index=["DataSet", "op", "perturb"], columns="method", values="Rank"
    ).reset_index()
    rng = np.random.default_rng(seed)
    comparisons = {}
    for competitor in ("scouter", "linearModel", "baseReg"):
        base = wide[["DataSet", "op", "perturb", "WitnessCell", competitor]].dropna().copy()
        base["difference"] = base.WitnessCell - base[competitor]
        groups = []
        for _, group in base.groupby("DataSet", sort=True):
            groups.append([
                cluster.difference.to_numpy()
                for _, cluster in group.groupby("perturb", sort=True)
            ])
        samples = np.empty(replicates)
        for repeat in range(replicates):
            values = []
            for clusters in groups:
                selected = rng.integers(0, len(clusters), size=len(clusters))
                values.extend(np.concatenate([clusters[index] for index in selected]))
            samples[repeat] = np.mean(values)
        comparisons[competitor] = {
            "rank_difference_witness_minus_competitor": float(base.difference.mean()),
            "cluster_bootstrap_95_ci": [
                float(np.quantile(samples, 0.025)),
                float(np.quantile(samples, 0.975)),
            ],
            "bootstrap_probability_witness_better": float(np.mean(samples < 0)),
            "operation_seed_units": int(len(base)),
            "unique_dataset_perturbation_clusters": int(
                base.groupby(["DataSet", "perturb"]).ngroups
            ),
        }
    report = {
        "status": "PASS_PAIRED_CLUSTER_BOOTSTRAP",
        "seed": seed,
        "replicates": replicates,
        "cluster": "dataset x perturbation; split-seed rows remain inside cluster",
        "comparisons": comparisons,
        "claim_boundary": "Uncertainty is over benchmark perturbation clusters, not cells; the four benchmark datasets themselves are fixed, not sampled from a population.",
    }
    (out / "paired_cluster_bootstrap.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    pd.DataFrame([
        {"competitor": competitor, **values}
        for competitor, values in comparisons.items()
    ]).to_csv(out / "paired_cluster_bootstrap.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--published-raw", type=Path, required=True)
    parser.add_argument("--published-top100", type=Path, required=True)
    parser.add_argument("--published-top5000", type=Path, required=True)
    parser.add_argument("--witness-raw", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    published_raw = pd.read_csv(args.published_raw)
    top100 = pd.read_csv(args.published_top100)
    top5000 = pd.read_csv(args.published_top5000)
    combo = {"Norman", "Wessels", "Schmidt", "Replogle_exp6"}
    published_raw = published_raw[published_raw.DataSet.isin(combo)].copy()
    top100 = top100[top100.DataSet.isin(combo)].copy()
    top5000 = top5000[top5000.DataSet.isin(combo)].copy()

    validation100 = validate_reconstruction(
        assemble_panel(100, published_raw, top100, None), top100, 100
    )
    validation5000 = validate_reconstruction(
        assemble_panel(5000, published_raw, top5000, None), top5000, 5000
    )
    if validation100["status"] != "PASS" or validation5000["status"] != "PASS":
        raise RuntimeError(json.dumps([validation100, validation5000], indent=2))

    witness = pd.concat([pd.read_csv(path) for path in args.witness_raw], ignore_index=True)
    full = pd.concat([
        assemble_panel(100, published_raw, top100, witness),
        assemble_panel(5000, published_raw, top5000, witness),
    ], ignore_index=True)
    full.to_csv(args.out / "published_plus_witness_per_operation.csv", index=False)
    write_metric_leaderboard(full, args.out)
    write_cluster_bootstrap(full, args.out)

    per_cell = full.groupby(
        ["DataSet", "DEG", "method"], as_index=False
    ).agg(
        mean_operation_rank=("Rank", "mean"),
        mean_metric_rank=("ac_rank", "mean"),
        mean_composite_score=("ac_score", "mean"),
        operations=("op", "nunique"),
    )
    per_cell["dataset_panel_position"] = per_cell.groupby(
        ["DataSet", "DEG"]
    ).mean_operation_rank.rank(ascending=True, method="min")
    per_cell.to_csv(args.out / "dataset_panel_leaderboard.csv", index=False)

    # Sensitivity analysis: each dataset and panel receives equal weight.  This
    # is useful for checking dataset dominance, but it is not the aggregation
    # used by the official Fig. 4 plotting code.
    equal_weight = per_cell.groupby("method", as_index=False).agg(
        mean_rank_4datasets_x_2panels=("mean_operation_rank", "mean"),
        mean_dataset_panel_position=("dataset_panel_position", "mean"),
        cells_won=("dataset_panel_position", lambda value: int(np.sum(value == 1))),
        cells=("dataset_panel_position", "size"),
    ).sort_values(["mean_rank_4datasets_x_2panels", "method"])
    equal_weight["overall_position"] = equal_weight.mean_rank_4datasets_x_2panels.rank(
        ascending=True, method="min"
    )
    equal_weight.to_csv(args.out / "formal_equal_weight_leaderboard.csv", index=False)

    dataset_equal_top100 = per_cell[per_cell.DEG == 100].groupby("method", as_index=False).agg(
        mean_rank_4datasets_top100=("mean_operation_rank", "mean"),
        mean_dataset_position_top100=("dataset_panel_position", "mean"),
        datasets_won_top100=("dataset_panel_position", lambda value: int(np.sum(value == 1))),
    ).sort_values(["mean_rank_4datasets_top100", "method"])
    dataset_equal_top100["top100_position"] = dataset_equal_top100.mean_rank_4datasets_top100.rank(
        ascending=True, method="min"
    )
    dataset_equal_top100.to_csv(
        args.out / "dataset_equal_weight_top100_sensitivity.csv", index=False
    )

    # Official Fig. 4 aggregation: average each method's per-operation Rank
    # over all perturbation/seed operations.  The top-100 panel is the paper's
    # primary six-metric accuracy evaluation.
    official_top100 = full[full.DEG == 100].groupby("method", as_index=False).agg(
        mean_operation_rank=("Rank", "mean"),
        mean_metric_rank=("ac_rank", "mean"),
        mean_composite_score=("ac_score", "mean"),
        operation_seed_units=("op", "size"),
    ).sort_values(["mean_operation_rank", "method"])
    official_top100["official_position"] = official_top100.mean_operation_rank.rank(
        ascending=True, method="min"
    )
    official_top100.to_csv(args.out / "formal_top100_leaderboard.csv", index=False)

    both_panels = full.groupby("method", as_index=False).agg(
        mean_operation_rank=("Rank", "mean"),
        mean_metric_rank=("ac_rank", "mean"),
        mean_composite_score=("ac_score", "mean"),
        operation_panel_units=("op", "size"),
    ).sort_values(["mean_operation_rank", "method"])
    both_panels["robustness_position"] = both_panels.mean_operation_rank.rank(
        ascending=True, method="min"
    )
    both_panels.to_csv(args.out / "formal_both_panels_leaderboard.csv", index=False)

    witness_equal = equal_weight[equal_weight.method == "WitnessCell"].iloc[0]
    witness_dataset_equal = dataset_equal_top100[
        dataset_equal_top100.method == "WitnessCell"
    ].iloc[0]
    witness_top100 = official_top100[official_top100.method == "WitnessCell"].iloc[0]
    witness_both = both_panels[both_panels.method == "WitnessCell"].iloc[0]
    eligible_published = sorted(set(published_raw.method.replace(NAME_MAP)))
    verdict = {
        "status": "PASS_FORMAL_AGGREGATION",
        "official_aggregation_reconstruction": [validation100, validation5000],
        "benchmark_total_methods": 27,
        "published_combo_track_methods": len(eligible_published),
        "published_combo_track_method_names": eligible_published,
        "new_combo_track_entrants_after_witness": len(eligible_published) + 1,
        "official_primary": "Top-100 perturbation-affected genes; six metrics; mean per-operation Rank over all four combination datasets, matching the official Fig. 4 aggregation code.",
        "official_primary_top100_position": int(witness_top100.official_position),
        "official_primary_top100_mean_rank": float(witness_top100.mean_operation_rank),
        "official_primary_operation_seed_units": int(witness_top100.operation_seed_units),
        "sota_official_primary": bool(witness_top100.official_position == 1),
        "both_panels_robustness_position": int(witness_both.robustness_position),
        "both_panels_robustness_mean_rank": float(witness_both.mean_operation_rank),
        "dataset_equal_weight_top100_position": int(witness_dataset_equal.top100_position),
        "dataset_equal_weight_top100_mean_rank": float(witness_dataset_equal.mean_rank_4datasets_top100),
        "dataset_equal_weight_both_panels_position": int(witness_equal.overall_position),
        "dataset_equal_weight_both_panels_mean_rank": float(witness_equal.mean_rank_4datasets_x_2panels),
        "dataset_panel_cells_won": int(witness_equal.cells_won),
        "sota_primary": bool(witness_top100.official_position == 1),
        "sota_top100": bool(witness_top100.official_position == 1),
        "claim_boundary": "Formal SOTA claim is limited to the genetic-combination track: 15 published eligible methods plus WitnessCell, four official combination datasets, three official split seeds, and the official top-100 six-metric aggregation. The benchmark contains 27 methods overall across incompatible tasks.",
    }
    (args.out / "formal_verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(official_top100.head(16).to_string(index=False))
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
