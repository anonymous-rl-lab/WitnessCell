#!/usr/bin/env python3
"""Freeze calibration mapping, thresholds, code and asset hashes before reveal."""
from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

from selective_core import fit_calibration, pair_balanced_weights, sha256, weighted_mean


ROOT = Path(__file__).resolve().parent
GRID = [0.3, 0.5, 0.7, 0.9]
PRIMARY = 0.5
HASHED_FILES = [
    "README.md",
    "PROTOCOL_DRAFT.json",
    "selective_core.py",
    "prepare_assets.py",
    "run_calibration_smoke.py",
    "freeze_gate.py",
    "reveal_formal_test.py",
    "audit.py",
    "assets/calibration_rows.csv",
    "assets/formal_risk_only.csv",
    "assets/ASSET_AUDIT.json",
    "sealed/formal_outcomes.csv",
    "results/calibration_smoke/verdict.json",
]


def main() -> None:
    outcome_path = ROOT / "sealed/formal_outcomes.csv"
    assert (outcome_path.stat().st_mode & 0o777) == 0
    asset_audit = json.loads((ROOT / "assets/ASSET_AUDIT.json").read_text())
    assert asset_audit["status"] == "PASS_PRE_REVEAL_ASSET_ISOLATION"
    assert sha256(outcome_path) == asset_audit["sealed_outcomes_sha256"]
    smoke = json.loads((ROOT / "results/calibration_smoke/verdict.json").read_text())
    assert smoke["formal_outcomes_read"] is False

    calibration = pd.read_csv(ROOT / "assets/calibration_rows.csv")
    mapping, thresholds = fit_calibration(
        calibration, "estimated_witness_risk", "estimated_realized_mse", GRID
    )
    geometry_mapping, geometry_thresholds = fit_calibration(
        calibration, "geometry_risk", "geometry_realized_mse", GRID
    )
    results_dir = ROOT / "results/frozen_calibration"
    results_dir.mkdir(parents=True, exist_ok=True)
    thresholds.to_csv(results_dir / "witness_thresholds.csv", index=False)
    geometry_thresholds.to_csv(results_dir / "geometry_thresholds.csv", index=False)
    (results_dir / "witness_mapping.json").write_text(json.dumps(mapping, indent=2) + "\n")
    (results_dir / "geometry_mapping.json").write_text(
        json.dumps(geometry_mapping, indent=2) + "\n"
    )

    risk_only = pd.read_csv(ROOT / "assets/formal_risk_only.csv")
    primary_threshold = float(
        thresholds.loc[thresholds.target_coverage == PRIMARY, "risk_threshold"].iloc[0]
    )
    weights = pair_balanced_weights(risk_only)
    planned_coverage = weighted_mean(
        (risk_only.estimated_witness_risk.to_numpy(float) <= primary_threshold).astype(float),
        weights,
    )
    draft = json.loads((ROOT / "PROTOCOL_DRAFT.json").read_text())
    draft["protocol_version"] = "21.0.0"
    draft["freeze_status"] = "FROZEN_BEFORE_FORMAL_OUTCOME_REVEAL"
    draft["asset_audit_sha256"] = sha256(ROOT / "assets/ASSET_AUDIT.json")
    draft["sealed_outcomes_sha256"] = sha256(outcome_path)
    draft["calibration_smoke_sha256"] = sha256(
        ROOT / "results/calibration_smoke/verdict.json"
    )
    draft["frozen_witness_thresholds"] = thresholds.to_dict(orient="records")
    draft["frozen_geometry_control_thresholds"] = geometry_thresholds.to_dict(
        orient="records"
    )
    draft["outcome_blind_primary_test_coverage"] = planned_coverage
    draft["formal_outcomes_read_at_freeze"] = False
    protocol_path = ROOT / "FROZEN_PROTOCOL.json"
    protocol_path.write_text(json.dumps(draft, indent=2) + "\n")
    protocol_hash = sha256(protocol_path)
    (ROOT / "FROZEN_PROTOCOL.sha256").write_text(
        f"{protocol_hash}  FROZEN_PROTOCOL.json\n"
    )

    frozen_files = HASHED_FILES + [
        "results/frozen_calibration/witness_thresholds.csv",
        "results/frozen_calibration/geometry_thresholds.csv",
        "results/frozen_calibration/witness_mapping.json",
        "results/frozen_calibration/geometry_mapping.json",
        "FROZEN_PROTOCOL.json",
        "FROZEN_PROTOCOL.sha256",
    ]
    lines = [f"{sha256(ROOT / relative)}  {relative}" for relative in frozen_files]
    (ROOT / "PRE_REVEAL_MANIFEST.sha256").write_text("\n".join(lines) + "\n")
    os.chmod(outcome_path, 0)
    print(json.dumps({
        "status": "FROZEN_BEFORE_FORMAL_OUTCOME_REVEAL",
        "protocol_sha256": protocol_hash,
        "primary_risk_threshold": primary_threshold,
        "outcome_blind_planned_test_coverage": planned_coverage,
        "formal_outcomes_read": False,
    }, indent=2))


if __name__ == "__main__":
    main()
