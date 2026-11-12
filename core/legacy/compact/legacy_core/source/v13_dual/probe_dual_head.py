#!/usr/bin/env python3
"""Training-only candidate screen for the WitnessCell dual endpoint head.

This script never opens test-condition means.  It ranks predeclared recipes
only from known-single pseudo-LOO rows and writes enough per-split evidence to
freeze one recipe before a blind test audit.
"""
from __future__ import annotations

import argparse
import csv
import json
import pickle
from pathlib import Path

import numpy as np

from dual_head import fit_dual_head
from identity_head import EPS
from official_split import (
    filter_gears_supported_conditions,
    load_gears_supported_genes,
    make_official_split,
)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def candidates(stage: str) -> list[dict]:
    rows = []
    if stage == "background":
        # Keep the directional object fixed to the one-coordinate self anchor.
        # This screen asks only whether a GO-neighbour residual belongs in the
        # dense background and, if so, how many neighbours are stable under
        # training-single LOO.  No held-out test mean is opened here.
        rows.append({
            "dense_mode": "mean_only",
            "sparse_mode": "self_only",
            "support_mode": "go_residual",
            "support_k": 1,
            "go_top_k": 40,
        })
        for go_top_k in (5, 10, 20, 40, 80):
            rows.append({
                "dense_mode": "mean_go",
                "sparse_mode": "self_only",
                "support_mode": "go_residual",
                "support_k": 1,
                "go_top_k": go_top_k,
            })
        return rows
    if stage == "restored":
        # Restore the v12 mean+GO+self dense endpoint program exactly, then
        # probe only an additional training-side sparse directional residual.
        # k=40 is the frozen v12 GO neighbourhood; only sparse support varies.
        rows.append({
            "dense_mode": "mean_go_self",
            "sparse_mode": "self_only",
            "support_mode": "go_residual",
            "support_k": 1,
            "go_top_k": 40,
        })
        for sparse_mode in ("go_only", "go_self", "scale_dense"):
            for support_mode in (
                "go_residual", "global_activity",
                "go_de_vote", "global_de_frequency",
            ):
                for support_k in (25, 50, 100, 200, 400):
                    rows.append({
                        "dense_mode": "mean_go_self",
                        "sparse_mode": sparse_mode,
                        "support_mode": support_mode,
                        "support_k": support_k,
                        "go_top_k": 40,
                    })
        return rows
    if stage == "joint_restore":
        return [{
            "dense_mode": "mean_go",
            "sparse_mode": "joint_self",
            "support_mode": "go_residual",
            "support_k": 1,
            "go_top_k": 40,
        }]
    if stage == "joint_k":
        return [{
            "dense_mode": "mean_go",
            "sparse_mode": "joint_self",
            "support_mode": "go_residual",
            "support_k": 1,
            "go_top_k": go_top_k,
        } for go_top_k in (3, 10, 20, 40, 80)]
    if stage == "coarse":
        for go_top_k in (10, 20, 40):
            for support_mode in ("go_residual", "go_program", "hybrid"):
                for support_k in (50, 100, 200):
                    rows.append({
                        "dense_mode": "mean_go",
                        "sparse_mode": "go_self",
                        "support_mode": support_mode,
                        "support_k": support_k,
                        "go_top_k": go_top_k,
                    })
        # Minimal ablations: no dense GO and no sparse GO.
        for support_k in (50, 100, 200):
            rows.append({
                "dense_mode": "mean_only",
                "sparse_mode": "go_self",
                "support_mode": "hybrid",
                "support_k": support_k,
                "go_top_k": 20,
            })
        rows.append({
            "dense_mode": "mean_go",
            "sparse_mode": "self_only",
            "support_mode": "hybrid",
            "support_k": 100,
            "go_top_k": 20,
        })
        return rows
    if stage == "focused":
        # mean+self is the stable first-stage identity head.  The second stage
        # adds only a masked GO residual, so the direction gate cannot silently
        # switch the generic self anchor on and off.
        for go_top_k in (10, 20, 40):
            for dense_mode in ("mean_self", "mean_go_self"):
                for support_mode in (
                    "go_residual", "go_program", "hybrid",
                    "go_de_vote", "global_de_frequency", "de_hybrid",
                ):
                    for support_k in (50, 100, 150, 200, 300):
                        rows.append({
                            "dense_mode": dense_mode,
                            "sparse_mode": "go_only",
                            "support_mode": support_mode,
                            "support_k": support_k,
                            "go_top_k": go_top_k,
                        })
        return rows
    if stage == "de_vote":
        for go_top_k in (10, 20, 40):
            for dense_mode in ("mean_self", "mean_go_self"):
                for support_mode in (
                    "go_de_vote", "global_de_frequency", "de_hybrid"
                ):
                    for support_k in (50, 100, 150, 200, 300):
                        rows.append({
                            "dense_mode": dense_mode,
                            "sparse_mode": "go_only",
                            "support_mode": support_mode,
                            "support_k": support_k,
                            "go_top_k": go_top_k,
                        })
        return rows
    if stage == "calibration":
        for go_top_k in (10, 20, 40):
            for support_mode in (
                "go_program", "global_activity",
                "go_de_vote", "global_de_frequency",
            ):
                for support_k in (50, 100, 200, 400):
                    rows.append({
                        "dense_mode": "mean_go_self",
                        "sparse_mode": "scale_dense",
                        "support_mode": support_mode,
                        "support_k": support_k,
                        "go_top_k": go_top_k,
                    })
        return rows
    raise ValueError(stage)


def deployed_aggregate(head) -> dict:
    records = list(head.records)
    baseline_all = sum(row["baseline_all_gene_mse"] for row in records)
    baseline_top = sum(row["baseline_top100_mse"] for row in records)
    if not head.dense_active:
        deployed_all = baseline_all
        deployed_top = baseline_top
        pcc_delta = 0.0
    elif head.sparse_active:
        deployed_all = sum(row["final_all_gene_mse"] for row in records)
        deployed_top = sum(row["final_top100_mse"] for row in records)
        pcc_delta = float(np.mean([
            row["final_top100_pcc"] - row["baseline_top100_pcc"]
            for row in records
        ]))
    else:
        deployed_all = sum(row["dense_all_gene_mse"] for row in records)
        deployed_top = sum(row["dense_top100_mse"] for row in records)
        pcc_delta = float(np.mean([
            row["dense_top100_pcc"] - row["baseline_top100_pcc"]
            for row in records
        ]))
    return {
        "baseline_all": baseline_all,
        "deployed_all": deployed_all,
        "baseline_top": baseline_top,
        "deployed_top": deployed_top,
        "pcc_delta": pcc_delta,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, action="append", required=True)
    parser.add_argument("--gears-assets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--stage",
        choices=(
            "background", "restored", "joint_restore", "joint_k", "coarse", "focused", "de_vote",
            "calibration",
        ),
        default="coarse",
    )
    parser.add_argument("--nested-gate", action="store_true")
    parser.add_argument("--seeds", type=int, nargs="+", default=(1, 2, 3))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    with (args.gears_assets / "gene2go_all.pkl").open("rb") as stream:
        gene2go = pickle.load(stream)
    supported = load_gears_supported_genes(args.gears_assets)
    cache_payloads = []
    for cache_path in args.cache:
        data = np.load(cache_path, allow_pickle=False)
        conditions = data["conditions"].astype(str).tolist()
        cache_payloads.append({
            "dataset": str(data["dataset"]),
            "conditions": conditions,
            "genes": data["genes"].astype(str).tolist(),
            "means": dict(zip(conditions, data["means"].astype(float))),
            "variances": dict(zip(conditions, data["variances"].astype(float))),
            "counts": dict(zip(conditions, map(int, data["counts"].astype(np.int64)))),
            "official": filter_gears_supported_conditions(conditions, supported),
        })

    split_rows: list[dict] = []
    summary_rows: list[dict] = []
    for candidate_id, recipe in enumerate(candidates(args.stage), start=1):
        blocks = []
        dense_active = 0
        sparse_active = 0
        wessels_dense_active = 0
        wessels_sparse_active = 0
        for payload in cache_payloads:
            control = payload["means"]["control"]
            for seed in args.seeds:
                split = make_official_split(payload["official"], seed)
                effects = {
                    condition: payload["means"][condition] - control
                    for condition in split.train
                    if condition != "control" and "+" not in condition
                }
                head = fit_dual_head(
                    single_effects=effects,
                    single_means={node: payload["means"][node] for node in effects},
                    single_variances={node: payload["variances"][node] for node in effects},
                    single_counts={node: payload["counts"][node] for node in effects},
                    control_mean=control,
                    control_variance=payload["variances"]["control"],
                    control_count=payload["counts"]["control"],
                    genes=payload["genes"],
                    gene2go=gene2go,
                    nested_gate=args.nested_gate,
                    **recipe,
                )
                block = deployed_aggregate(head)
                blocks.append(block)
                dense_active += int(head.dense_active)
                sparse_active += int(head.sparse_active)
                if payload["dataset"] == "Wessels":
                    wessels_dense_active += int(head.dense_active)
                    wessels_sparse_active += int(head.sparse_active)
                split_rows.append({
                    "candidate_id": candidate_id,
                    **recipe,
                    "nested_gate": int(args.nested_gate),
                    "dataset": payload["dataset"],
                    "seed": seed,
                    "known_single_count": head.known_single_count,
                    "dense_active": int(head.dense_active),
                    "sparse_active": int(head.sparse_active),
                    "dense_weights": ";".join(f"{value:.12g}" for value in head.dense_weights),
                    "sparse_weights": ";".join(f"{value:.12g}" for value in head.sparse_weights),
                    "dense_all_gene_gain_mean": head.dense_all_gene_gain_mean,
                    "dense_all_gene_gain_lower": head.dense_all_gene_gain_lower,
                    "final_all_gene_gain_mean": head.final_all_gene_gain_mean,
                    "final_all_gene_gain_lower": head.final_all_gene_gain_lower,
                    "final_top100_gain_mean": head.final_top100_gain_mean,
                    "final_top100_gain_lower": head.final_top100_gain_lower,
                    "sparse_incremental_top100_gain_mean": head.sparse_incremental_top100_gain_mean,
                    "sparse_incremental_top100_gain_lower": head.sparse_incremental_top100_gain_lower,
                    "final_top100_pcc_delta_mean": head.final_top100_pcc_delta_mean,
                    "final_top100_pcc_delta_lower": head.final_top100_pcc_delta_lower,
                    "sparse_incremental_top100_pcc_delta_mean": head.sparse_incremental_top100_pcc_delta_mean,
                    "sparse_incremental_top100_pcc_delta_lower": head.sparse_incremental_top100_pcc_delta_lower,
                    "deployed_all_gene_gain": 1.0 - block["deployed_all"] / max(block["baseline_all"], EPS),
                    "deployed_top100_gain": 1.0 - block["deployed_top"] / max(block["baseline_top"], EPS),
                    "deployed_top100_pcc_delta": block["pcc_delta"],
                })
        baseline_all = sum(block["baseline_all"] for block in blocks)
        deployed_all = sum(block["deployed_all"] for block in blocks)
        baseline_top = sum(block["baseline_top"] for block in blocks)
        deployed_top = sum(block["deployed_top"] for block in blocks)
        summary_rows.append({
            "candidate_id": candidate_id,
            **recipe,
            "nested_gate": int(args.nested_gate),
            "split_count": len(blocks),
            "dense_active_splits": dense_active,
            "sparse_active_splits": sparse_active,
            "wessels_dense_active_splits": wessels_dense_active,
            "wessels_sparse_active_splits": wessels_sparse_active,
            "pooled_deployed_all_gene_gain": 1.0 - deployed_all / max(baseline_all, EPS),
            "pooled_deployed_top100_gain": 1.0 - deployed_top / max(baseline_top, EPS),
            "mean_deployed_top100_pcc_delta": float(np.mean([block["pcc_delta"] for block in blocks])),
        })

    # Training-only lexicographic selection: preserve Wessels fallback, then
    # maximize top-100 gain, using all-gene gain and recipe id as tie breakers.
    eligible = [row for row in summary_rows if row["wessels_dense_active_splits"] == 0]
    ranked = sorted(
        eligible,
        key=lambda row: (
            -row["pooled_deployed_top100_gain"],
            -row["pooled_deployed_all_gene_gain"],
            row["candidate_id"],
        ),
    )
    selected = ranked[0] if ranked else None
    write_csv(args.out / "dual_head_training_by_split.csv", split_rows)
    write_csv(args.out / "dual_head_training_summary.csv", summary_rows)
    verdict = {
        "status": "PASS_DUAL_HEAD_TRAINING_ONLY_SCREEN",
        "stage": args.stage,
        "nested_gate": args.nested_gate,
        "candidate_count": len(summary_rows),
        "selection_uses_test": False,
        "selection_contract": "Wessels dense fallback first; pooled deployed training-LOO top100 gain; all-gene tie break",
        "selected": selected,
    }
    (args.out / "verdict.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))
    print("TOP TRAINING-ONLY CANDIDATES")
    for row in ranked[:10]:
        print(row)


if __name__ == "__main__":
    main()
