"""Frozen and calibration-derived accept/abstain policies."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .exceptions import ValidationError

FROZEN_NORMAN_THRESHOLD = 0.0923227147328771


def weighted_quantile(
    values: Sequence[float] | np.ndarray,
    quantile: float,
    weights: Sequence[float] | np.ndarray,
) -> float:
    values_array = np.asarray(values, dtype=float)
    weights_array = np.asarray(weights, dtype=float)
    if values_array.ndim != 1 or weights_array.shape != values_array.shape:
        raise ValidationError("values and weights must be equal-length vectors")
    if values_array.size == 0:
        raise ValidationError("values and weights must not be empty")
    if not np.all(np.isfinite(values_array)) or not np.all(np.isfinite(weights_array)):
        raise ValidationError("values and weights must contain only finite values")
    if not np.isfinite(quantile) or not 0.0 <= quantile <= 1.0:
        raise ValidationError("quantile must be in [0, 1]")
    if np.any(weights_array < 0) or weights_array.sum() <= 0:
        raise ValidationError("weights must be non-negative with positive sum")
    order = np.argsort(values_array, kind="mergesort")
    ordered_values = values_array[order]
    cdf = np.cumsum(weights_array[order]) / weights_array.sum()
    index = int(np.searchsorted(cdf, quantile, side="left"))
    return float(ordered_values[min(index, len(ordered_values) - 1)])


@dataclass(frozen=True)
class SelectivePolicy:
    """A risk-threshold policy whose only actions are accept and abstain."""

    threshold: float
    provenance: str = "user-calibrated"

    def __post_init__(self) -> None:
        if not np.isfinite(self.threshold) or self.threshold < 0:
            raise ValidationError("threshold must be finite and non-negative")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValidationError("provenance must be a non-empty string")

    @classmethod
    def frozen_norman(cls) -> "SelectivePolicy":
        """Return the retrospective Gate 21 threshold, without transport claims."""
        return cls(FROZEN_NORMAN_THRESHOLD, "gate21-norman-retrospective")

    @classmethod
    def calibrate(
        cls,
        risks: Sequence[float] | np.ndarray,
        *,
        coverage: float = 0.5,
        pair_weights: Sequence[float] | None = None,
    ) -> "SelectivePolicy":
        risk = np.asarray(risks, dtype=float)
        if not np.isfinite(coverage) or not 0.0 < coverage <= 1.0:
            raise ValidationError("coverage must be finite and in (0, 1]")
        weights = (
            np.ones(len(risk), dtype=float)
            if pair_weights is None
            else np.asarray(pair_weights, dtype=float)
        )
        return cls(
            weighted_quantile(risk, coverage, weights),
            f"calibration-pair-balanced-q{coverage:g}",
        )

    def decide(self, risks: Sequence[float] | np.ndarray) -> tuple[str, ...]:
        risk = np.asarray(risks, dtype=float)
        if risk.ndim != 1 or not np.all(np.isfinite(risk)):
            raise ValidationError("risks must be a finite vector")
        return tuple("accept" if value <= self.threshold else "abstain" for value in risk)
