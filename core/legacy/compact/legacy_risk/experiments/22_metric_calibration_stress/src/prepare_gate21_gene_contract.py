#!/usr/bin/env python3
"""Regenerate and seal Gate 21 gene-level arrays without scoring WMSE.

This is a pre-freeze asset-reconstruction step.  It reproduces the immutable
Gate 07 predictions used by Gate 21, verifies every previously revealed scalar
loss/risk value, aligns the candidate-blind Norman truth-side weight contract,
and writes gene-level vectors.  It deliberately does not evaluate the new
weighted loss; that operation is reserved for formal Phase D after freeze.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


THRESHOLD = 0.0923227147328771
SEEDS = range(100, 130)
# Gate 07 was created under an older BLAS stack and did not archive prediction
# vectors.  Re-solving the same float32 residual system is required to recover
# those vectors; a 5-ppm parity bound is strict enough to detect algorithm or
# split drift while allowing documented cross-BLAS solve variation.
RTOL = 5e-6
ATOL = 5e-9


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity_hash(values: list[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--norman-contract", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    gate07 = repo / "experiments/07_estimated_witness_norman"
    sys.path.insert(0, str(gate07))
    from estimated_witness import (  # type: ignore
        fit_predict,
        incidence,
        make_kernel_matrices,
        row_cosine,
        safe_edge_split,
    )
    from run_norman_estimated_witness import (  # type: ignore
        load_single_profiles,
        symmetric_pair_features,
    )

    residual_path = repo / "experiments/04_formal_gate/results/pseudobulk_residuals.npz"
    formal_dir = gate07 / "results/formal_30split"
    revealed_path = (
        repo
        / "experiments/21_frozen_selective_prediction/results/formal_reveal/revealed_target_rows.csv"
    )
    calibration_path = gate07 / "results/smoke/target_rows.csv"
    metric_path = formal_dir / "per_seed_metrics.csv"

    data = np.load(residual_path, allow_pickle=True)
    pairs = np.asarray(data["pairs"], dtype=str)
    response = np.asarray(data["residual"], dtype=np.float64)
    output_genes = np.asarray(data["genes"], dtype=str)
    nodes = sorted({node for pair in pairs for node in pair.split("+")})
    node_id = {node: index for index, node in enumerate(nodes)}
    edges = np.asarray(
        [(node_id[pair.split("+")[0]], node_id[pair.split("+")[1]]) for pair in pairs],
        dtype=int,
    )
    design = incidence(edges, len(nodes))
    singles = load_single_profiles(repo / "data/norman", output_genes, nodes)
    features = symmetric_pair_features(singles, edges, components=12)

    revealed = pd.read_csv(revealed_path)
    revealed_ids = revealed["query_id"].astype(str).tolist()
    if not revealed["query_id"].is_unique or len(revealed) != 213:
        raise RuntimeError("Gate 21 revealed identity contract is not 213 unique rows")
    calibration_pairs = set(pd.read_csv(calibration_path)["pair"].astype(str))
    metrics = pd.read_csv(metric_path)

    regenerated: dict[str, dict] = {}
    for seed in SEEDS:
        all_indices = np.arange(len(edges))
        train, test = safe_edge_split(
            all_indices,
            edges,
            len(nodes),
            0.20,
            np.random.default_rng(7411 + seed),
        )
        local_metrics = metrics.loc[metrics["seed"] == seed]
        if set(local_metrics["strategy"]) != {
            "estimated_witness",
            "geometry_only",
            "empirical_oracle",
        }:
            raise RuntimeError(f"unexpected Gate 07 strategies for seed {seed}")

        fitted = {}
        for strategy in ("estimated_witness", "geometry_only"):
            row = local_metrics.loc[local_metrics["strategy"] == strategy].iloc[0]
            kernels = make_kernel_matrices(
                train,
                test,
                design,
                features,
                float(row["selected_length_factor"]),
            )
            fitted[strategy] = fit_predict(
                response[train],
                kernels,
                float(row["selected_rho"]),
                float(row["selected_noise_ratio"]),
            )

        for position, edge_index in enumerate(test):
            pair = str(pairs[edge_index])
            if pair in calibration_pairs:
                continue
            query_id = f"{seed:03d}::{pair}"
            regenerated[query_id] = {
                "seed": seed,
                "pair": pair,
                "truth": response[edge_index],
                "estimated_prediction": fitted["estimated_witness"].prediction[position],
                "geometry_prediction": fitted["geometry_only"].prediction[position],
                "estimated_risk": float(fitted["estimated_witness"].risk[position]),
                "geometry_risk": float(fitted["geometry_only"].risk[position]),
            }

    if set(regenerated) != set(revealed_ids):
        missing = sorted(set(revealed_ids) - set(regenerated))
        extra = sorted(set(regenerated) - set(revealed_ids))
        raise RuntimeError(f"Gate 21 identity mismatch: missing={missing[:3]} extra={extra[:3]}")

    estimated_prediction = np.stack([regenerated[q]["estimated_prediction"] for q in revealed_ids])
    geometry_prediction = np.stack([regenerated[q]["geometry_prediction"] for q in revealed_ids])
    truth = np.stack([regenerated[q]["truth"] for q in revealed_ids])
    estimated_risk = np.asarray([regenerated[q]["estimated_risk"] for q in revealed_ids])
    geometry_risk = np.asarray([regenerated[q]["geometry_risk"] for q in revealed_ids])
    estimated_mse = np.mean((estimated_prediction - truth) ** 2, axis=1)
    geometry_mse = np.mean((geometry_prediction - truth) ** 2, axis=1)
    estimated_cosine = row_cosine(estimated_prediction, truth)
    geometry_cosine = row_cosine(geometry_prediction, truth)

    comparisons = {
        "estimated_risk": (estimated_risk, revealed["estimated_witness_risk"].to_numpy(float)),
        "geometry_risk": (geometry_risk, revealed["geometry_risk"].to_numpy(float)),
        "estimated_mse": (estimated_mse, revealed["estimated_realized_mse"].to_numpy(float)),
        "geometry_mse": (geometry_mse, revealed["geometry_realized_mse"].to_numpy(float)),
        "estimated_cosine": (
            estimated_cosine,
            revealed["estimated_residual_cosine"].to_numpy(float),
        ),
        "geometry_cosine": (
            geometry_cosine,
            revealed["geometry_residual_cosine"].to_numpy(float),
        ),
    }
    verification = {}
    for name, (actual, expected) in comparisons.items():
        absolute = np.abs(actual - expected)
        local_atol = 5e-6 if name.endswith("cosine") else ATOL
        verification[name] = {
            "max_abs_error": float(np.max(absolute)),
            "max_rel_error": float(
                np.max(absolute / np.maximum(np.abs(expected), np.finfo(float).eps))
            ),
            "rtol": RTOL,
            "atol": local_atol,
            "allclose": bool(np.allclose(actual, expected, rtol=RTOL, atol=local_atol)),
        }
        if not verification[name]["allclose"]:
            raise RuntimeError(
                f"Gate 21 regeneration failed scalar parity for {name}: "
                f"{verification[name]}"
            )

    # The previously revealed risk column—not a cross-BLAS reconstruction—is
    # the immutable Gate 21 deployment score used by the formal stress test.
    frozen_estimated_risk = revealed["estimated_witness_risk"].to_numpy(float)
    frozen_geometry_risk = revealed["geometry_risk"].to_numpy(float)
    accepted = frozen_estimated_risk <= THRESHOLD
    expected_accepted = revealed["accepted_primary"].astype(bool).to_numpy()
    if not np.array_equal(accepted, expected_accepted):
        raise RuntimeError("Gate 21 accepted set changed during regeneration")

    # Align candidate-blind first-half truth weights.  No weighted score is
    # evaluated here; this only seals the arrays Phase D is allowed to read.
    contract = np.load(args.norman_contract, allow_pickle=False)
    contract_genes = np.asarray(contract["genes"], dtype=str)
    contract_conditions = np.asarray(contract["conditions"], dtype=str)
    contract_weights = np.asarray(contract["weights"], dtype=np.float64)
    contract_evaluable = np.asarray(contract["weight_evaluable"], dtype=bool)
    gene_index = {gene: index for index, gene in enumerate(contract_genes)}
    missing_genes = [gene for gene in output_genes if gene not in gene_index]
    overlap_mask = np.asarray([gene in gene_index for gene in output_genes], dtype=bool)
    output_index = np.asarray([gene_index[gene] for gene in output_genes[overlap_mask]], dtype=int)
    condition_index = {condition: index for index, condition in enumerate(contract_conditions)}
    pair_values = revealed["pair"].astype(str).to_numpy()
    missing_pairs = sorted(set(pair_values) - set(condition_index))
    condition_rows = np.asarray([condition_index.get(pair, -1) for pair in pair_values], dtype=int)
    scoring_evaluable = condition_rows >= 0
    scoring_evaluable[scoring_evaluable] &= contract_evaluable[condition_rows[scoring_evaluable]]
    full_condition_weights = np.zeros((len(revealed), len(contract_genes)), dtype=np.float64)
    full_condition_weights[scoring_evaluable] = contract_weights[condition_rows[scoring_evaluable]]
    canonical_weights = np.zeros((len(revealed), len(output_genes)), dtype=np.float64)
    canonical_weights[:, overlap_mask] = full_condition_weights[:, output_index]
    scoring_evaluable &= canonical_weights.sum(axis=1) > 0
    if not np.all(np.isfinite(canonical_weights)):
        raise RuntimeError("Gate 21 canonical weights contain non-finite values")
    overlap_weight_fraction = np.full(len(revealed), np.nan, dtype=np.float64)
    overlap_weight_fraction[scoring_evaluable] = (
        canonical_weights[scoring_evaluable].sum(axis=1)
        / full_condition_weights[scoring_evaluable].sum(axis=1)
    )

    np.savez_compressed(
        args.out,
        query_ids=np.asarray(revealed_ids, dtype=str),
        seeds=revealed["seed"].to_numpy(dtype=np.int64),
        pairs=np.asarray(pair_values, dtype=str),
        genes=np.asarray(output_genes, dtype=str),
        estimated_prediction=estimated_prediction.astype(np.float64),
        geometry_prediction=geometry_prediction.astype(np.float64),
        truth=truth.astype(np.float64),
        canonical_weights=canonical_weights.astype(np.float64),
        scoring_evaluable=scoring_evaluable.astype(bool),
        estimated_witness_risk=frozen_estimated_risk.astype(np.float64),
        geometry_risk=frozen_geometry_risk.astype(np.float64),
        accepted_primary=accepted.astype(bool),
        threshold=np.asarray(THRESHOLD, dtype=np.float64),
    )
    report = {
        "status": "PASS_GATE21_GENE_LEVEL_REGENERATION",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "rows": len(revealed_ids),
        "pairs": int(revealed["pair"].nunique()),
        "genes": int(output_genes.size),
        "canonical_gene_overlap": int(overlap_mask.sum()),
        "canonical_gene_missing_zero_weight": int(len(missing_genes)),
        "canonical_weight_evaluable_rows": int(scoring_evaluable.sum()),
        "canonical_weight_non_evaluable_rows": int((~scoring_evaluable).sum()),
        "canonical_weight_non_evaluable_pairs": missing_pairs,
        "overlap_weight_fraction": {
            "minimum": float(np.nanmin(overlap_weight_fraction)),
            "mean": float(np.nanmean(overlap_weight_fraction)),
            "median": float(np.nanmedian(overlap_weight_fraction)),
        },
        "accepted": int(accepted.sum()),
        "rejected": int((~accepted).sum()),
        "threshold_repr": repr(float(THRESHOLD)),
        "query_identity_sha256": identity_hash(revealed_ids),
        "accepted_identity_sha256": identity_hash(sorted(np.asarray(revealed_ids)[accepted].tolist())),
        "scalar_parity": verification,
        "regeneration_tolerance": {
            "risk_and_mse": {"rtol": RTOL, "atol": ATOL},
            "cosine": {"rtol": RTOL, "atol": 5e-6}
        },
        "formal_risk_source": "previously frozen Gate 21 revealed rows; regenerated risk is parity diagnostic only",
        "weighted_loss_evaluated_before_freeze": False,
        "gene_contract_sha256": sha256(args.out),
        "source_assets": {
            "residuals_sha256": sha256(residual_path),
            "gate07_metrics_sha256": sha256(metric_path),
            "gate21_revealed_rows_sha256": sha256(revealed_path),
            "norman_phase_m_contract_sha256": sha256(args.norman_contract),
        },
    }
    report_path = args.out.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
