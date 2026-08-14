"""Validated public data objects used by WitnessCell."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from numbers import Real
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from .core import parse_condition
from .exceptions import ValidationError


def _finite_vector(value: Any, width: int, label: str) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != (width,):
        raise ValidationError(
            f"{label} must have shape ({width},), got {array.shape}"
        )
    if not np.all(np.isfinite(array)):
        raise ValidationError(f"{label} contains NaN or infinite values")
    frozen = array.copy()
    frozen.setflags(write=False)
    return frozen


def _string_keyed_mapping(value: Mapping[Any, Any], label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValidationError(f"{label} must be a mapping")
    normalized = {str(key): item for key, item in value.items()}
    if len(normalized) != len(value):
        raise ValidationError(f"{label} contains keys that collide after string conversion")
    return normalized


def _positive_count(value: object, label: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValidationError(f"{label} must be a positive integer")
    number = float(value)
    if not np.isfinite(number) or not number.is_integer() or number < 1:
        raise ValidationError(f"{label} must be a positive integer")
    return int(number)


@dataclass(frozen=True)
class ConditionMoments:
    """Per-condition means, population variances, and cell counts.

    The class is deliberately condition-level: WitnessCell 0.1 predicts
    condition means and does not treat cells as independent replicates.
    """

    genes: tuple[str, ...]
    means: Mapping[str, np.ndarray]
    variances: Mapping[str, np.ndarray]
    counts: Mapping[str, int]

    def __post_init__(self) -> None:
        normalized_genes = tuple(map(str, self.genes))
        if not normalized_genes:
            raise ValidationError("genes must not be empty")
        if any(not gene.strip() for gene in normalized_genes):
            raise ValidationError("genes must contain non-empty names")
        if len(set(normalized_genes)) != len(normalized_genes):
            raise ValidationError("genes must be unique")

        raw_means = _string_keyed_mapping(self.means, "means")
        raw_variances = _string_keyed_mapping(self.variances, "variances")
        raw_counts = _string_keyed_mapping(self.counts, "counts")
        conditions = tuple(raw_means)
        if not conditions:
            raise ValidationError("means must contain at least one condition")
        if set(conditions) != set(raw_variances):
            raise ValidationError("means and variances must have identical conditions")
        if set(conditions) != set(raw_counts):
            raise ValidationError("means and counts must have identical conditions")

        width = len(normalized_genes)
        normalized_means: dict[str, np.ndarray] = {}
        normalized_variances: dict[str, np.ndarray] = {}
        normalized_counts: dict[str, int] = {}
        for condition in conditions:
            try:
                parse_condition(condition)
            except ValueError as exc:
                raise ValidationError(str(exc)) from exc
            normalized_means[condition] = _finite_vector(
                raw_means[condition], width, f"means[{condition!r}]"
            )
            variance = _finite_vector(
                raw_variances[condition], width, f"variances[{condition!r}]"
            )
            if np.any(variance < 0):
                raise ValidationError(f"variances[{condition!r}] contains negatives")
            normalized_variances[condition] = variance
            normalized_counts[condition] = _positive_count(
                raw_counts[condition], f"counts[{condition!r}]"
            )

        object.__setattr__(self, "genes", normalized_genes)
        object.__setattr__(self, "means", MappingProxyType(normalized_means))
        object.__setattr__(self, "variances", MappingProxyType(normalized_variances))
        object.__setattr__(self, "counts", MappingProxyType(normalized_counts))

    @classmethod
    def from_mappings(
        cls,
        *,
        genes: Sequence[str],
        means: Mapping[str, np.ndarray],
        variances: Mapping[str, np.ndarray],
        counts: Mapping[str, int],
    ) -> "ConditionMoments":
        return cls(tuple(genes), means, variances, counts)

    @classmethod
    def from_npz(cls, path: str | Path) -> "ConditionMoments":
        """Load the documented pickle-free condition-moment interchange file."""
        with np.load(Path(path), allow_pickle=False) as archive:
            required = {"genes", "conditions", "means", "variances", "counts"}
            missing = required - set(archive.files)
            if missing:
                raise ValidationError(f"moment archive missing keys: {sorted(missing)}")
            genes = archive["genes"].astype(str)
            conditions = archive["conditions"].astype(str)
            means = np.asarray(archive["means"], dtype=float)
            variances = np.asarray(archive["variances"], dtype=float)
            counts = np.asarray(archive["counts"])
        if means.shape != (len(conditions), len(genes)):
            raise ValidationError("means matrix does not match conditions × genes")
        if variances.shape != means.shape:
            raise ValidationError("variances matrix must match means")
        if counts.shape != (len(conditions),):
            raise ValidationError("counts must have one value per condition")
        return cls.from_mappings(
            genes=genes,
            means=dict(zip(conditions, means, strict=True)),
            variances=dict(zip(conditions, variances, strict=True)),
            counts=dict(zip(conditions, counts, strict=True)),
        )

    def to_npz(self, path: str | Path) -> None:
        """Write a pickle-free, non-object condition-moment archive."""
        conditions = tuple(self.means)
        np.savez_compressed(
            Path(path),
            genes=np.asarray(self.genes, dtype=str),
            conditions=np.asarray(conditions, dtype=str),
            means=np.stack([self.means[c] for c in conditions]),
            variances=np.stack([self.variances[c] for c in conditions]),
            counts=np.asarray([self.counts[c] for c in conditions], dtype=np.int64),
        )

    def require(self, conditions: Sequence[str]) -> None:
        missing = sorted(set(map(str, conditions)) - set(self.means))
        if missing:
            raise ValidationError(f"conditions are absent from moments: {missing}")


@dataclass(frozen=True)
class SplitSpec:
    """Explicit training and validation roles for leakage-safe fitting."""

    train_conditions: tuple[str, ...]
    validation_conditions: tuple[str, ...]

    def __post_init__(self) -> None:
        train = tuple(map(str, self.train_conditions))
        validation = tuple(map(str, self.validation_conditions))
        if len(train) != len(set(train)):
            raise ValidationError("train_conditions contains duplicates")
        if len(validation) != len(set(validation)):
            raise ValidationError("validation_conditions contains duplicates")
        overlap = sorted(set(train) & set(validation))
        if overlap:
            raise ValidationError(
                f"training and validation conditions overlap: {overlap}"
            )
        for condition in train + validation:
            try:
                parse_condition(condition)
            except ValueError as exc:
                raise ValidationError(str(exc)) from exc
        object.__setattr__(self, "train_conditions", train)
        object.__setattr__(self, "validation_conditions", validation)

    @classmethod
    def create(
        cls,
        train_conditions: Sequence[str],
        validation_conditions: Sequence[str] = (),
    ) -> "SplitSpec":
        return cls(tuple(train_conditions), tuple(validation_conditions))


@dataclass(frozen=True)
class PredictionBatch:
    """Condition-mean predictions with explicit factorized fallback output."""

    conditions: tuple[str, ...]
    genes: tuple[str, ...]
    means: np.ndarray
    effects: np.ndarray
    factorized_means: np.ndarray
    factorized_effects: np.ndarray
    risks: np.ndarray | None = None
    decisions: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        conditions = tuple(map(str, self.conditions))
        genes = tuple(map(str, self.genes))
        if not conditions:
            raise ValidationError("conditions must not be empty")
        if not genes or any(not gene.strip() for gene in genes):
            raise ValidationError("genes must contain non-empty names")
        if len(set(genes)) != len(genes):
            raise ValidationError("genes must be unique")
        for condition in conditions:
            try:
                parse_condition(condition)
            except ValueError as exc:
                raise ValidationError(str(exc)) from exc
        object.__setattr__(self, "conditions", conditions)
        object.__setattr__(self, "genes", genes)

        expected = (len(conditions), len(genes))
        for name in ("means", "effects", "factorized_means", "factorized_effects"):
            array = np.asarray(getattr(self, name), dtype=float)
            if array.shape != expected or not np.all(np.isfinite(array)):
                raise ValidationError(f"{name} must be a finite conditions × genes matrix")
            frozen = array.copy()
            frozen.setflags(write=False)
            object.__setattr__(self, name, frozen)
        if self.risks is not None:
            risk = np.asarray(self.risks, dtype=float)
            if risk.shape != (len(conditions),) or not np.all(np.isfinite(risk)):
                raise ValidationError("risks must be one finite value per condition")
            frozen_risk = risk.copy()
            frozen_risk.setflags(write=False)
            object.__setattr__(self, "risks", frozen_risk)
        if self.decisions is not None:
            decisions = tuple(map(str, self.decisions))
            if len(decisions) != len(conditions) or not set(decisions) <= {"accept", "abstain"}:
                raise ValidationError("decisions must be accept/abstain for every condition")
            object.__setattr__(self, "decisions", decisions)

    def with_decisions(
        self, risks: Sequence[float], decisions: Sequence[str]
    ) -> "PredictionBatch":
        risk_array = np.asarray(risks, dtype=float)
        decision_values = tuple(map(str, decisions))
        if risk_array.shape != (len(self.conditions),):
            raise ValidationError("risks must contain one value per condition")
        if len(decision_values) != len(self.conditions):
            raise ValidationError("decisions must contain one value per condition")
        return replace(self, risks=risk_array, decisions=decision_values)

    def to_npz(self, path: str | Path) -> None:
        """Write a deployment-safe archive with no truth or object arrays."""
        payload: dict[str, np.ndarray] = {
            "conditions": np.asarray(self.conditions, dtype=str),
            "genes": np.asarray(self.genes, dtype=str),
            "means": np.asarray(self.means, dtype=float),
            "effects": np.asarray(self.effects, dtype=float),
            "factorized_means": np.asarray(self.factorized_means, dtype=float),
            "factorized_effects": np.asarray(self.factorized_effects, dtype=float),
        }
        if self.risks is not None:
            payload["risks"] = np.asarray(self.risks, dtype=float)
        if self.decisions is not None:
            payload["decisions"] = np.asarray(self.decisions, dtype=str)
        np.savez_compressed(Path(path), **payload)  # type: ignore[arg-type]
