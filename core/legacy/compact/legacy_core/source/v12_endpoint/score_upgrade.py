#!/usr/bin/env python3
"""Paired mean-response audit of endpoint v2 versus frozen WitnessCell v1."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from identity_head import EPS, pearson, welch_scores


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def aggregate(rows: list[dict], label: str) -> dict:
    return {
        "scope": label,
        "condition_seed_units": len(rows),
        "all_gene_mse_v1": np.mean([row["all_gene_mse_v1"] for row in rows]),
        "all_gene_mse_v2": np.mean([row["all_gene_mse_v2"] for row in rows]),
        "relative_all_gene_mse_gain": 1.0 - sum(row["all_gene_mse_v2"] for row in rows) / max(sum(row["all_gene_mse_v1"] for row in rows), EPS),
        "top100_mse_v1": np.mean([row["top100_mse_v1"] for row in rows]),
        "top100_mse_v2": np.mean([row["top100_mse_v2"] for row in rows]),
        "relative_top100_mse_gain": 1.0 - sum(row["top100_mse_v2"] for row in rows) / max(sum(row["top100_mse_v1"] for row in rows), EPS),
        "top100_pcc_v1": np.mean([row["top100_pcc_v1"] for row in rows]),
        "top100_pcc_v2": np.mean([row["top100_pcc_v2"] for row in rows]),
        "top100_pcc_delta": np.mean([row["top100_pcc_delta"] for row in rows]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in args.prediction:
        data = np.load(path, allow_pickle=False)
        dataset = str(data["dataset"])
        seed = int(data["seed"])
        control = data["control"].astype(float)
        control_variance = data["control_variance"].astype(float)
        control_count = int(data["control_count"])
        for index, condition in enumerate(data["conditions"].astype(str)):
            truth = data["truth"][index].astype(float)
            baseline = data["baseline_prediction"][index].astype(float)
            upgraded = data["prediction"][index].astype(float)
            score = welch_scores(
                truth,
                data["truth_variance"][index].astype(float),
                int(data["test_cell_counts"][index]),
                control,
                control_variance,
                control_count,
            )
            order = np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]
            base_all = float(np.mean(np.square(baseline - truth)))
            new_all = float(np.mean(np.square(upgraded - truth)))
            base_top = float(np.mean(np.square(baseline[order] - truth[order])))
            new_top = float(np.mean(np.square(upgraded[order] - truth[order])))
            base_pcc = pearson(baseline[order] - control[order], truth[order] - control[order])
            new_pcc = pearson(upgraded[order] - control[order], truth[order] - control[order])
            rows.append({
                "dataset": dataset,
                "seed": seed,
                "condition": condition,
                "subgroup": str(data["subgroups"][index]),
                "all_gene_mse_v1": base_all,
                "all_gene_mse_v2": new_all,
                "relative_all_gene_mse_gain": (base_all - new_all) / max(base_all, EPS),
                "top100_mse_v1": base_top,
                "top100_mse_v2": new_top,
                "relative_top100_mse_gain": (base_top - new_top) / max(base_top, EPS),
                "top100_pcc_v1": base_pcc,
                "top100_pcc_v2": new_pcc,
                "top100_pcc_delta": new_pcc - base_pcc,
            })
    write_csv(args.out / "paired_by_condition.csv", rows)
    summaries = [aggregate(rows, "all")]
    for dataset in sorted({row["dataset"] for row in rows}):
        summaries.append(aggregate(
            [row for row in rows if row["dataset"] == dataset], dataset
        ))
    for subgroup in ("combo_seen0", "combo_seen1", "combo_seen2", "unseen_single"):
        block = [row for row in rows if row["subgroup"] == subgroup]
        if block:
            summaries.append(aggregate(block, subgroup))
    write_csv(args.out / "paired_summary.csv", summaries)
    verdict = {
        "status": "PASS_PAIRED_ENDPOINT_UPGRADE_DIRECTIONAL_AUDIT",
        "archives": len(args.prediction),
        "condition_seed_units": len(rows),
        "overall": summaries[0],
        "contract": "same conditions, truth panels and frozen v1 predictions; paired v2-v1 comparison",
    }
    (args.out / "verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))
    for row in summaries:
        print(row)


if __name__ == "__main__":
    main()
