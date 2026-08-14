#!/usr/bin/env python3
"""Formal Phase D: regrade fixed Gate 21 identities with canonical WMSE."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from metric_core import WeightVector, wmse


THRESHOLD = 0.0923227147328771
REPLICATES = 20_000
BOOTSTRAP_SEED = 20260811
PERMUTATION_SEED = 20260812


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--gene-contract", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    args.out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(repo / "experiments/21_frozen_selective_prediction"))
    from selective_core import (  # type: ignore
        cluster_bootstrap_ratio,
        evaluate_threshold,
        pair_balanced_weights,
        weighted_quantile,
        within_seed_permutation,
    )

    with np.load(args.gene_contract, allow_pickle=False) as data:
        query_ids = data["query_ids"].astype(str)
        seeds = data["seeds"].astype(np.int64)
        pairs = data["pairs"].astype(str)
        genes = data["genes"].astype(str)
        estimated_prediction = data["estimated_prediction"].astype(np.float64)
        geometry_prediction = data["geometry_prediction"].astype(np.float64)
        truth = data["truth"].astype(np.float64)
        weights = data["canonical_weights"].astype(np.float64)
        scoring_evaluable = data["scoring_evaluable"].astype(bool)
        risk = data["estimated_witness_risk"].astype(np.float64)
        geometry_risk = data["geometry_risk"].astype(np.float64)
        stored_accepted = data["accepted_primary"].astype(bool)
        stored_threshold = float(data["threshold"])
    if np.float64(stored_threshold).tobytes() != np.float64(THRESHOLD).tobytes():
        raise RuntimeError("Gate 21 threshold bit pattern drift")
    accepted = risk <= THRESHOLD
    if not np.array_equal(accepted, stored_accepted):
        raise RuntimeError("Gate 21 accepted identities drifted")

    estimated_wmse = np.full(len(query_ids), np.nan, dtype=np.float64)
    geometry_wmse = np.full(len(query_ids), np.nan, dtype=np.float64)
    for index in np.flatnonzero(scoring_evaluable):
        weight = WeightVector(genes, weights[index])
        estimated_wmse[index] = wmse(estimated_prediction[index], truth[index], weight)
        geometry_wmse[index] = wmse(geometry_prediction[index], truth[index], weight)
    full_frame = pd.DataFrame(
        {
            "query_id": query_ids,
            "seed": seeds,
            "pair": pairs,
            "estimated_witness_risk": risk,
            "geometry_risk": geometry_risk,
            "estimated_wmse": estimated_wmse,
            "geometry_wmse": geometry_wmse,
            "scoring_evaluable": scoring_evaluable,
            "accepted_primary": accepted,
        }
    )
    frame = full_frame.loc[full_frame["scoring_evaluable"]].copy()
    if len(frame) / len(full_frame) < 0.90:
        raise RuntimeError("Gate 21 complete-case coverage fell below frozen 90% floor")
    primary = evaluate_threshold(
        frame, "estimated_witness_risk", "estimated_wmse", THRESHOLD
    )
    bootstrap = cluster_bootstrap_ratio(
        frame,
        "estimated_witness_risk",
        "estimated_wmse",
        THRESHOLD,
        replicates=REPLICATES,
        seed=BOOTSTRAP_SEED,
    )
    permuted = within_seed_permutation(
        frame,
        "estimated_witness_risk",
        "estimated_wmse",
        THRESHOLD,
        replicates=REPLICATES,
        seed=PERMUTATION_SEED,
    )
    permutation_p = float(
        (1 + np.sum(permuted <= float(primary["accepted_mse"]))) / (len(permuted) + 1)
    )
    criteria = {
        "coverage_non_degenerate": 0.35 <= primary["pair_balanced_coverage"] <= 0.65,
        "practical_effect": primary["accepted_over_all_mse"] <= 0.80,
        "cluster_uncertainty": float(np.quantile(bootstrap, 0.975)) < 1.0,
        "random_selection_control": permutation_p < 0.05,
        "rejected_separation": primary["accepted_mse"] < primary["rejected_mse"],
    }

    pair_weights = pair_balanced_weights(frame)
    geometry_threshold = weighted_quantile(
        geometry_risk,
        float(primary["pair_balanced_coverage"]),
        pair_weights,
    )
    geometry = evaluate_threshold(
        frame, "geometry_risk", "estimated_wmse", geometry_threshold
    )
    geometry_bootstrap = cluster_bootstrap_ratio(
        frame,
        "geometry_risk",
        "estimated_wmse",
        geometry_threshold,
        replicates=REPLICATES,
        seed=BOOTSTRAP_SEED,
    )
    geometry_permuted = within_seed_permutation(
        frame,
        "geometry_risk",
        "estimated_wmse",
        geometry_threshold,
        replicates=REPLICATES,
        seed=PERMUTATION_SEED,
    )
    geometry_p = float(
        (1 + np.sum(geometry_permuted <= float(geometry["accepted_mse"])))
        / (len(geometry_permuted) + 1)
    )
    geometry_criteria = {
        "coverage_non_degenerate": 0.35 <= geometry["pair_balanced_coverage"] <= 0.65,
        "practical_effect": geometry["accepted_over_all_mse"] <= 0.80,
        "cluster_uncertainty": float(np.quantile(geometry_bootstrap, 0.975)) < 1.0,
        "random_selection_control": geometry_p < 0.05,
        "rejected_separation": geometry["accepted_mse"] < geometry["rejected_mse"],
    }
    result = {
        "status": "PASS_PHASE_D_EXECUTION",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "full_rows": int(len(full_frame)),
        "full_pairs": int(full_frame["pair"].nunique()),
        "evaluable_rows": int(len(frame)),
        "evaluable_pairs": int(frame["pair"].nunique()),
        "evaluable_fraction": float(len(frame) / len(full_frame)),
        "accepted_full": int(accepted.sum()),
        "rejected_full": int((~accepted).sum()),
        "threshold_repr": repr(float(THRESHOLD)),
        "primary": primary,
        "bootstrap_accepted_over_all_ci95": [
            float(np.quantile(bootstrap, 0.025)),
            float(np.quantile(bootstrap, 0.975)),
        ],
        "within_seed_random_selection_p_one_sided": permutation_p,
        "risk_wmse_spearman": float(spearmanr(risk, estimated_wmse).statistic),
        "criteria": criteria,
        "GATE21_WMSE": (
            "PASS" if bool(np.all(scoring_evaluable)) and all(criteria.values())
            else "FAIL" if bool(np.all(scoring_evaluable))
            else "NOT_ADJUDICATED_FULL_213"
        ),
        "GATE21_WMSE_COMPLETE_CASE": "PASS" if all(criteria.values()) else "FAIL",
        "geometry_matched_coverage": {
            "threshold": geometry_threshold,
            "result": geometry,
            "bootstrap_accepted_over_all_ci95": [
                float(np.quantile(geometry_bootstrap, 0.025)),
                float(np.quantile(geometry_bootstrap, 0.975)),
            ],
            "within_seed_random_selection_p_one_sided": geometry_p,
            "criteria": geometry_criteria,
            "all_five_pass": all(geometry_criteria.values()),
            "witness_minus_geometry_selective_wmse": float(
                primary["accepted_mse"] - geometry["accepted_mse"]
            ),
        },
        "source_gene_contract_sha256": sha256(args.gene_contract),
    }
    table_path = args.out / "gate21_wmse_rows.csv"
    array_path = args.out / "gate21_inference_arrays.npz"
    full_frame.to_csv(table_path, index=False)
    np.savez_compressed(
        array_path,
        witness_bootstrap_ratio=bootstrap,
        witness_random_selection_wmse=permuted,
        geometry_bootstrap_ratio=geometry_bootstrap,
        geometry_random_selection_wmse=geometry_permuted,
    )
    result["artifacts"] = {
        "rows_sha256": sha256(table_path),
        "inference_arrays_sha256": sha256(array_path),
    }
    (args.out / "decision_verdict.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
