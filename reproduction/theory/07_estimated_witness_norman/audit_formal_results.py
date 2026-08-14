#!/usr/bin/env python3
"""Numerical and split-integrity audit for the formal Norman gate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=Path("results/formal_30split"))
    args = parser.parse_args()
    verdict = json.loads((args.results / "verdict.json").read_text())
    targets = pd.read_csv(args.results / "target_rows.csv")
    covariance_files = sorted((args.results / "estimated_covariances").glob("seed_*.npz"))
    audits = []
    maximum_symmetry_error = 0.0
    minimum_eigenvalue = float("inf")
    minimum_single_target_schur = float("inf")
    for path in covariance_files:
        data = np.load(path, allow_pickle=True)
        train = data["train_indices"].astype(int)
        test = data["test_indices"].astype(int)
        K = data["K_hat"].astype(float)
        cross = data["k_t_hat"].astype(float)
        target_variance = data["k_tt_hat"].astype(float)
        risk = data["witness_risk"].astype(float)
        symmetry_error = float(np.max(np.abs(K - K.T)))
        eigenvalue = float(np.linalg.eigvalsh(K).min())
        maximum_symmetry_error = max(maximum_symmetry_error, symmetry_error)
        minimum_eigenvalue = min(minimum_eigenvalue, eigenvalue)
        regularization = max(float(np.trace(K) / max(len(K), 1)) * 1e-8, 1e-12)
        inverse = np.linalg.inv(K + regularization * np.eye(len(K)))
        schur = target_variance - np.sum(cross * (inverse @ cross), axis=0)
        minimum_single_target_schur = min(minimum_single_target_schur, float(schur.min()))
        seed = int(path.stem.split("_")[-1])
        recorded_pairs = set(targets.loc[targets.seed == seed, "pair"])
        audits.append({
            "seed": seed,
            "train_count": len(train),
            "test_count": len(test),
            "train_test_overlap": int(len(np.intersect1d(train, test))),
            "test_pair_names_match": bool(recorded_pairs == set(data["test_pairs"].astype(str))),
            "K_symmetric": bool(symmetry_error <= 1e-6),
            "K_psd": bool(eigenvalue >= -1e-5),
            "risk_finite_positive": bool(np.all(np.isfinite(risk)) and np.all(risk > 0)),
        })
    all_checks = bool(
        verdict["gate_pass"]
        and len(covariance_files) == 30
        and len(targets.seed.unique()) == 30
        and all(
            row["train_test_overlap"] == 0
            and row["test_pair_names_match"]
            and row["K_symmetric"]
            and row["K_psd"]
            and row["risk_finite_positive"]
            for row in audits
        )
    )
    output = {
        "formal_gate_pass": bool(verdict["gate_pass"]),
        "covariance_files": len(covariance_files),
        "formal_seeds": len(targets.seed.unique()),
        "maximum_K_symmetry_error": maximum_symmetry_error,
        "minimum_K_eigenvalue": minimum_eigenvalue,
        "minimum_individual_target_latent_schur": minimum_single_target_schur,
        "note": "Schur diagnostic uses K alone with a numerical inverse and is not the total Witness Risk, which also includes geometry and noise.",
        "split_audits": audits,
        "audit_pass": all_checks,
    }
    (args.results / "audit.json").write_text(json.dumps(output, indent=2))
    print(json.dumps({key: value for key, value in output.items() if key != "split_audits"}, indent=2))
    if not all_checks:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

