#!/usr/bin/env python3
"""Prepare calibration, risk-only test and sealed formal outcomes."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd

from selective_core import sha256


CALIBRATION_COLUMNS = [
    "seed", "pair", "estimated_witness_risk", "geometry_risk",
    "estimated_realized_mse", "geometry_realized_mse",
    "estimated_residual_cosine", "geometry_residual_cosine",
    "selected_rho", "selected_length_factor", "selected_noise_ratio",
]
RISK_ONLY_COLUMNS = [
    "query_id", "seed", "pair", "estimated_witness_risk", "geometry_risk",
    "selected_rho", "selected_length_factor", "selected_noise_ratio",
]
SEALED_COLUMNS = [
    "query_id", "estimated_realized_mse", "geometry_realized_mse",
    "estimated_residual_cosine", "geometry_residual_cosine",
    "reliability_weight",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.out.resolve()
    calibration_source = (
        args.repo / "experiments/07_estimated_witness_norman/results/smoke/target_rows.csv"
    )
    formal_source = (
        args.repo
        / "experiments/07_estimated_witness_norman/results/formal_30split/target_rows.csv"
    )
    calibration = pd.read_csv(calibration_source)
    formal = pd.read_csv(formal_source)
    assert set(CALIBRATION_COLUMNS).issubset(calibration.columns)
    assert set(RISK_ONLY_COLUMNS[3:] + SEALED_COLUMNS[1:]).issubset(formal.columns)
    assert set(calibration.seed.unique()) == set(range(5))
    assert set(formal.seed.unique()) == set(range(100, 130))

    calibration_pairs = set(calibration.pair.astype(str))
    eligible = formal.loc[~formal.pair.astype(str).isin(calibration_pairs)].copy()
    eligible.insert(0, "query_id", eligible.apply(
        lambda row: f"{int(row.seed):03d}::{row.pair}", axis=1
    ))
    assert eligible.query_id.is_unique
    assert not set(eligible.pair.astype(str)) & calibration_pairs
    assert eligible.pair.nunique() >= 25
    assert len(eligible) >= 150

    assets = root / "assets"
    sealed = root / "sealed"
    assets.mkdir(exist_ok=True)
    sealed.mkdir(exist_ok=True)
    calibration_path = assets / "calibration_rows.csv"
    risk_path = assets / "formal_risk_only.csv"
    outcome_path = sealed / "formal_outcomes.csv"
    calibration[CALIBRATION_COLUMNS].to_csv(calibration_path, index=False)
    eligible[RISK_ONLY_COLUMNS].to_csv(risk_path, index=False)
    eligible[SEALED_COLUMNS].to_csv(outcome_path, index=False)
    os.chmod(outcome_path, 0)

    forbidden_tokens = ("mse", "cosine", "reliability", "oracle", "truth")
    assert not any(
        token in column.lower()
        for column in pd.read_csv(risk_path, nrows=0).columns
        for token in forbidden_tokens
    )
    report = {
        "status": "PASS_PRE_REVEAL_ASSET_ISOLATION",
        "calibration_source": str(calibration_source),
        "calibration_source_sha256": sha256(calibration_source),
        "formal_source": str(formal_source),
        "formal_source_sha256": sha256(formal_source),
        "calibration_rows": int(len(calibration)),
        "calibration_pairs": int(len(calibration_pairs)),
        "formal_eligible_rows": int(len(eligible)),
        "formal_eligible_pairs": int(eligible.pair.nunique()),
        "pair_identity_overlap": 0,
        "risk_only_sha256": sha256(risk_path),
        "sealed_outcomes_sha256": sha256(outcome_path),
        "sealed_outcomes_mode": oct(outcome_path.stat().st_mode & 0o777),
        "risk_only_columns": RISK_ONLY_COLUMNS,
        "forbidden_test_columns_present_in_risk_only": False,
    }
    (assets / "ASSET_AUDIT.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        key: report[key] for key in (
            "status", "calibration_rows", "calibration_pairs",
            "formal_eligible_rows", "formal_eligible_pairs", "pair_identity_overlap",
            "forbidden_test_columns_present_in_risk_only",
        )
    }, indent=2))


if __name__ == "__main__":
    main()
