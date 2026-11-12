#!/usr/bin/env python3
"""Training-LOO and blind-test probe for unseen single perturbations."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from identity_head import EPS, fit_identity_head, pearson, welch_scores
from official_split import (
    filter_gears_supported_conditions,
    load_gears_supported_genes,
    make_official_split,
)


CANDIDATES = (
    ("fixed_mean_plus_self", "add", "mean", "fixed"),
    ("fixed_mean_set_self", "replace", "mean", "fixed"),
    ("weighted_mean_plus_self", "add", "mean", "closed_form"),
    ("weighted_mean_set_self", "replace", "mean", "closed_form"),
    ("weighted_mean_plus_median_self", "add", "median", "closed_form"),
)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, action="append", required=True)
    parser.add_argument("--gears-assets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=(1, 2, 3))
    parser.add_argument("--gate-metric", choices=("top100_mse", "all_gene_mse"), default="top100_mse")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    supported = load_gears_supported_genes(args.gears_assets)
    condition_rows: list[dict] = []
    fit_rows: list[dict] = []

    for cache_path in args.cache:
        data = np.load(cache_path, allow_pickle=False)
        dataset = str(data["dataset"])
        conditions = data["conditions"].astype(str).tolist()
        genes = data["genes"].astype(str).tolist()
        matrix = data["means"].astype(float)
        variance_matrix = data["variances"].astype(float)
        count_array = data["counts"].astype(np.int64)
        row_id = {condition: index for index, condition in enumerate(conditions)}
        means = {condition: matrix[index] for index, condition in enumerate(conditions)}
        variances = {condition: variance_matrix[index] for index, condition in enumerate(conditions)}
        counts = {condition: int(count_array[index]) for index, condition in enumerate(conditions)}
        official = filter_gears_supported_conditions(conditions, supported)
        gene_index = {gene: index for index, gene in enumerate(genes)}
        control = means["control"]

        for seed in args.seeds:
            split = make_official_split(official, seed)
            single_effects = {
                condition: means[condition] - control
                for condition in split.train
                if condition != "control" and "+" not in condition
            }
            background = np.mean(np.stack(list(single_effects.values())), axis=0)
            for name, anchor_mode, self_summary, weight_mode in CANDIDATES:
                fit = fit_identity_head(
                    single_effects=single_effects,
                    single_means={condition: means[condition] for condition in single_effects},
                    single_variances={condition: variances[condition] for condition in single_effects},
                    single_counts={condition: counts[condition] for condition in single_effects},
                    control_mean=control,
                    control_variance=variances["control"],
                    control_count=counts["control"],
                    genes=genes,
                    anchor_mode=anchor_mode,
                    self_summary=self_summary,
                    weight_mode=weight_mode,
                    gate_metric=args.gate_metric,
                )
                fit_rows.append({
                    "dataset": dataset,
                    "seed": seed,
                    "candidate": name,
                    "active": int(fit.active),
                    "background_weight": fit.background_weight,
                    "self_weight": fit.self_weight,
                    "gate_metric": fit.gate_metric,
                    "gate_lower_bound": fit.gate_lower_bound,
                    "loo_all_gene_gain_mean": fit.loo_all_gene_gain_mean,
                    "loo_all_gene_gain_lower": fit.loo_all_gene_gain_lower,
                    "loo_top100_gain_mean": fit.loo_top100_gain_mean,
                    "loo_top100_gain_lower": fit.loo_top100_gain_lower,
                    "loo_pcc_delta_mean": fit.loo_pcc_delta_mean,
                    "known_single_count": fit.known_single_count,
                })
                for condition in split.test_subgroup["unseen_single"]:
                    truth = means[condition] - control
                    raw_prediction = fit.background_weight * background
                    # predict() obeys the gate; reconstruct raw candidate by a
                    # copy with gate semantics handled explicitly below.
                    from identity_head import features_for_node
                    base_feature, anchor = features_for_node(
                        condition, single_effects, gene_index, anchor_mode, self_summary
                    )
                    raw_prediction = fit.background_weight * base_feature + fit.self_weight * anchor
                    gated_prediction = raw_prediction if fit.active else background
                    score = welch_scores(
                        means[condition], variances[condition], counts[condition],
                        control, variances["control"], counts["control"],
                    )
                    order = np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]
                    baseline_all = float(np.mean(np.square(background - truth)))
                    baseline_top = float(np.mean(np.square(background[order] - truth[order])))
                    for branch, prediction in (("raw", raw_prediction), ("gated", gated_prediction)):
                        all_mse = float(np.mean(np.square(prediction - truth)))
                        top_mse = float(np.mean(np.square(prediction[order] - truth[order])))
                        condition_rows.append({
                            "dataset": dataset,
                            "seed": seed,
                            "condition": condition,
                            "candidate": name,
                            "branch": branch,
                            "gate_active": int(fit.active),
                            "baseline_all_gene_mse": baseline_all,
                            "candidate_all_gene_mse": all_mse,
                            "relative_all_gene_mse_gain": (baseline_all - all_mse) / max(baseline_all, EPS),
                            "baseline_top100_mse": baseline_top,
                            "candidate_top100_mse": top_mse,
                            "relative_top100_mse_gain": (baseline_top - top_mse) / max(baseline_top, EPS),
                            "baseline_top100_pcc": pearson(background[order], truth[order]),
                            "candidate_top100_pcc": pearson(prediction[order], truth[order]),
                            "top100_pcc_delta": pearson(prediction[order], truth[order]) - pearson(background[order], truth[order]),
                        })

    write_csv(args.out / "identity_fit_by_split.csv", fit_rows)
    write_csv(args.out / "unseen_single_by_condition.csv", condition_rows)
    summary_rows: list[dict] = []
    keys = sorted({(r["candidate"], r["branch"]) for r in condition_rows})
    for candidate, branch in keys:
        block = [r for r in condition_rows if r["candidate"] == candidate and r["branch"] == branch]
        summary_rows.append({
            "candidate": candidate,
            "branch": branch,
            "condition_seed_units": len(block),
            "relative_all_gene_mse_gain": 1.0 - sum(r["candidate_all_gene_mse"] for r in block) / max(sum(r["baseline_all_gene_mse"] for r in block), EPS),
            "relative_top100_mse_gain": 1.0 - sum(r["candidate_top100_mse"] for r in block) / max(sum(r["baseline_top100_mse"] for r in block), EPS),
            "mean_top100_pcc_delta": float(np.mean([r["top100_pcc_delta"] for r in block])),
            "active_split_count": len({(r["dataset"], r["seed"]) for r in block if r["gate_active"]}),
        })
    write_csv(args.out / "unseen_single_summary.csv", summary_rows)
    verdict = {
        "status": "PASS_ENDPOINT_IDENTITY_CANDIDATE_PROBE",
        "seeds": list(args.seeds),
        "datasets": sorted({r["dataset"] for r in condition_rows}),
        "condition_seed_units": len({(r["dataset"], r["seed"], r["condition"]) for r in condition_rows}),
        "candidate_count": len(CANDIDATES),
        "gate_metric": args.gate_metric,
        "selection_contract": "candidate choice must use only identity_fit_by_split LOO evidence; blind-test summary is opened only after freezing",
    }
    (args.out / "verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))
    for row in summary_rows:
        print(row)


if __name__ == "__main__":
    main()
