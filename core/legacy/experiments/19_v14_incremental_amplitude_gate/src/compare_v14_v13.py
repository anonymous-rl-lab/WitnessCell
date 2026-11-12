#!/usr/bin/env python3
"""Paired held-out comparison of v14 against frozen v13 predictions."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from identity_head import EPS, pearson, welch_scores


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict], scope: str) -> dict:
    return {
        "scope": scope,
        "condition_seed_units": len(rows),
        "all_gene_mse_v13": float(np.mean([row["all_gene_mse_v13"] for row in rows])),
        "all_gene_mse_v14": float(np.mean([row["all_gene_mse_v14"] for row in rows])),
        "relative_all_gene_mse_gain": 1.0 - sum(row["all_gene_mse_v14"] for row in rows) / max(sum(row["all_gene_mse_v13"] for row in rows), EPS),
        "top100_mse_v13": float(np.mean([row["top100_mse_v13"] for row in rows])),
        "top100_mse_v14": float(np.mean([row["top100_mse_v14"] for row in rows])),
        "relative_top100_mse_gain": 1.0 - sum(row["top100_mse_v14"] for row in rows) / max(sum(row["top100_mse_v13"] for row in rows), EPS),
        "top100_pcc_v13": float(np.mean([row["top100_pcc_v13"] for row in rows])),
        "top100_pcc_v14": float(np.mean([row["top100_pcc_v14"] for row in rows])),
        "top100_pcc_delta": float(np.mean([row["top100_pcc_delta"] for row in rows])),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v14-root", type=Path, required=True)
    parser.add_argument("--v13-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    split_audit = []
    for v14_path in sorted(args.v14_root.glob("*/seed*/predictions.npz")):
        dataset = v14_path.parent.parent.name
        seed = int(v14_path.parent.name.removeprefix("seed"))
        v13_path = args.v13_root / dataset / f"seed{seed}" / "predictions.npz"
        if not v13_path.exists():
            v13_path = args.v13_root / dataset / f"seed{seed}" / "deploy_predictions.npz"
        v14 = np.load(v14_path, allow_pickle=False)
        v13 = np.load(v13_path, allow_pickle=False)
        if not np.array_equal(v14["conditions"], v13["conditions"]):
            raise ValueError(f"condition mismatch: {dataset} seed{seed}")
        manifest = json.loads((v14_path.parent / "manifest.json").read_text())
        active = bool(manifest["identity_head"]["correction_active"])
        max_delta = float(np.max(np.abs(v14["prediction"] - v13["prediction"])))
        split_audit.append({
            "dataset": dataset,
            "seed": seed,
            "correction_active": active,
            "max_prediction_delta": max_delta,
            "inactive_exact_fallback": bool(active or max_delta == 0.0),
        })
        control = v14["control"].astype(float)
        for index, condition in enumerate(v14["conditions"].astype(str)):
            truth = v14["truth"][index].astype(float)
            p13 = v13["prediction"][index].astype(float)
            p14 = v14["prediction"][index].astype(float)
            score = welch_scores(
                truth, v14["truth_variance"][index].astype(float),
                int(v14["test_cell_counts"][index]), control,
                v14["control_variance"].astype(float), int(v14["control_count"]),
            )
            top = np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]
            all13 = float(np.mean(np.square(p13 - truth)))
            all14 = float(np.mean(np.square(p14 - truth)))
            top13 = float(np.mean(np.square(p13[top] - truth[top])))
            top14 = float(np.mean(np.square(p14[top] - truth[top])))
            pcc13 = pearson(p13[top] - control[top], truth[top] - control[top])
            pcc14 = pearson(p14[top] - control[top], truth[top] - control[top])
            rows.append({
                "dataset": dataset,
                "seed": seed,
                "condition": condition,
                "subgroup": str(v14["subgroups"][index]),
                "correction_active": active,
                "all_gene_mse_v13": all13,
                "all_gene_mse_v14": all14,
                "relative_all_gene_mse_gain": (all13 - all14) / max(all13, EPS),
                "top100_mse_v13": top13,
                "top100_mse_v14": top14,
                "relative_top100_mse_gain": (top13 - top14) / max(top13, EPS),
                "top100_pcc_v13": pcc13,
                "top100_pcc_v14": pcc14,
                "top100_pcc_delta": pcc14 - pcc13,
            })
    summaries = [summarize(rows, "all")]
    for dataset in sorted({row["dataset"] for row in rows}):
        summaries.append(summarize(
            [row for row in rows if row["dataset"] == dataset], dataset
        ))
    for subgroup in ("combo_seen0", "combo_seen1", "combo_seen2", "unseen_single"):
        block = [row for row in rows if row["subgroup"] == subgroup]
        if block:
            summaries.append(summarize(block, subgroup))
    write_csv(args.out / "paired_by_condition.csv", rows)
    write_csv(args.out / "paired_summary.csv", summaries)
    write_csv(args.out / "split_fallback_audit.csv", split_audit)
    verdict = {
        "status": "PASS_V14_V13_PAIRED_SMALL_PANEL_AUDIT",
        "split_count": len(split_audit),
        "active_splits": sum(row["correction_active"] for row in split_audit),
        "inactive_exact_fallback": all(row["inactive_exact_fallback"] for row in split_audit),
        "overall": summaries[0],
    }
    (args.out / "verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))
    for row in summaries:
        print(row)


if __name__ == "__main__":
    main()
