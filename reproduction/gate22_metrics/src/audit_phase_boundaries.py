#!/usr/bin/env python3
"""Static and filesystem audit of Experiment 22 reveal boundaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.experiment.resolve()
    phase_p = (root / "src/phase_p_prediction_stress.py").read_text()
    phase_m = (root / "src/phase_m_score_controls.py").read_text()
    formal_runner = (root / "run_formal_experiment22.sh").read_text()
    old_name = "METRIC_CALIBRATION_STRESS_PROTOCOL.md"
    checks = {
        "phase_p_verifies_contract_hash": "contract drift" in phase_p
        and "split ledger drift" in phase_p,
        "phase_p_verifies_phase_m_manifest": "sealed Phase M output drift" in phase_p
        and "phase_m_manifest_sha256" in phase_p,
        "phase_p_has_no_contract_write": "savez(args.contracts" not in phase_p
        and "write_text(args.contracts" not in phase_p,
        "phase_m_has_no_candidate_argument": "prediction" not in " ".join(
            line for line in phase_m.splitlines() if "add_argument" in line
        ),
        "formal_runner_enables_actual_firewall_variable": "E22_FORBIDDEN_READ_ROOTS" in formal_runner,
        "formal_runner_refuses_existing_output": "once-only runner refuses overwrite" in formal_runner,
        "formal_output_absent_before_freeze": not (root / "results/formal_e22").exists(),
        "superseded_draft_absent": not (root / old_name).exists(),
        "superseded_draft_not_in_runner": old_name not in formal_runner,
    }
    report = {
        "status": "PASS_PHASE_BOUNDARY_AUDIT" if all(checks.values()) else "FAIL_PHASE_BOUNDARY_AUDIT",
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
