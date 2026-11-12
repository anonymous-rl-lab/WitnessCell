#!/usr/bin/env python3
"""Resume Phase D with dense indices and exact complete-case geometry alignment."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def argument_value(name: str) -> str:
    try:
        return sys.argv[sys.argv.index(name) + 1]
    except (ValueError, IndexError) as exc:
        raise SystemExit(f"missing required argument {name}") from exc


repo = Path(argument_value("--repo")).resolve()
gene_contract = Path(argument_value("--gene-contract")).resolve()
with np.load(gene_contract, allow_pickle=False) as data:
    scoring_evaluable = data["scoring_evaluable"].astype(bool)

sys.path.insert(0, str(repo / "experiments/21_frozen_selective_prediction"))
sys.path.insert(0, str(repo / "experiments/22_metric_calibration_stress/src"))

import selective_core  # type: ignore  # noqa: E402


_frozen_within_seed_permutation = selective_core.within_seed_permutation
_frozen_weighted_quantile = selective_core.weighted_quantile


def _dense_index_within_seed_permutation(frame, *args, **kwargs):
    return _frozen_within_seed_permutation(frame.reset_index(drop=True), *args, **kwargs)


def _complete_case_weighted_quantile(values, quantile, weights):
    values_array = np.asarray(values)
    weights_array = np.asarray(weights)
    if values_array.size == scoring_evaluable.size and weights_array.size == int(scoring_evaluable.sum()):
        values_array = values_array[scoring_evaluable]
    return _frozen_weighted_quantile(values_array, quantile, weights_array)


selective_core.within_seed_permutation = _dense_index_within_seed_permutation
selective_core.weighted_quantile = _complete_case_weighted_quantile

import phase_d_decision_stress  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(phase_d_decision_stress.main())
