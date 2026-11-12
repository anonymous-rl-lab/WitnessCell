#!/usr/bin/env python3
"""Resume frozen Phase D with a dense positional index after complete-case filtering."""

from __future__ import annotations

import sys
from pathlib import Path


def argument_value(name: str) -> str:
    try:
        return sys.argv[sys.argv.index(name) + 1]
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"missing required argument {name}") from exc


repo = Path(argument_value("--repo")).resolve()
sys.path.insert(0, str(repo / "experiments/21_frozen_selective_prediction"))
sys.path.insert(0, str(repo / "experiments/22_metric_calibration_stress/src"))

import selective_core  # type: ignore  # noqa: E402


_frozen_within_seed_permutation = selective_core.within_seed_permutation


def _dense_index_within_seed_permutation(frame, *args, **kwargs):
    return _frozen_within_seed_permutation(frame.reset_index(drop=True), *args, **kwargs)


selective_core.within_seed_permutation = _dense_index_within_seed_permutation

import phase_d_decision_stress  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(phase_d_decision_stress.main())
