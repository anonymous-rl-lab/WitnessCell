"""High-level WitnessCell estimator API."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .core import (
    DEFAULT_ALPHA_GRID,
    DEFAULT_NOISE_GRID,
    FittedWitnessCore,
    fit_witness_core,
    parse_condition,
)
from .data import ConditionMoments, PredictionBatch, SplitSpec
from .exceptions import NotFittedError, ValidationError


@dataclass(frozen=True)
class WitnessCellConfig:
    """Configuration for the published v14 condition-mean contract."""

    control_label: str = "control"
    alpha_grid: tuple[float, ...] = DEFAULT_ALPHA_GRID
    noise_grid: tuple[float, ...] = DEFAULT_NOISE_GRID
    gate_confidence: float = 0.95

    def __post_init__(self) -> None:
        if not isinstance(self.control_label, str) or not self.control_label.strip():
            raise ValidationError("control_label must not be empty")
        try:
            parts = parse_condition(self.control_label)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
        if len(parts) != 1:
            raise ValidationError("control_label must be a single condition label")
        try:
            alpha_grid = tuple(map(float, self.alpha_grid))
            noise_grid = tuple(map(float, self.noise_grid))
            confidence = float(self.gate_confidence)
        except (TypeError, ValueError) as exc:
            raise ValidationError("configuration grids and confidence must be numeric") from exc
        if (
            not alpha_grid
            or not np.all(np.isfinite(alpha_grid))
            or any(value < 0 for value in alpha_grid)
        ):
            raise ValidationError("alpha_grid must contain finite non-negative values")
        if len(set(alpha_grid)) != len(alpha_grid):
            raise ValidationError("alpha_grid must not contain duplicates")
        if (
            not noise_grid
            or not np.all(np.isfinite(noise_grid))
            or any(value <= 0 for value in noise_grid)
        ):
            raise ValidationError("noise_grid must contain finite positive values")
        if len(set(noise_grid)) != len(noise_grid):
            raise ValidationError("noise_grid must not contain duplicates")
        if not np.isfinite(confidence) or not 0.5 < confidence < 1.0:
            raise ValidationError("gate_confidence must be between 0.5 and 1")
        object.__setattr__(self, "alpha_grid", alpha_grid)
        object.__setattr__(self, "noise_grid", noise_grid)
        object.__setattr__(self, "gate_confidence", confidence)

    @property
    def contract_id(self) -> str:
        if (
            self.control_label == "control"
            and self.alpha_grid == DEFAULT_ALPHA_GRID
            and self.noise_grid == DEFAULT_NOISE_GRID
            and self.gate_confidence == 0.95
        ):
            return "witnesscell-v14-frozen"
        return "custom-nonfrozen"


class WitnessCell:
    """Evidence-gated predictor for unmeasured genetic combinations.

    ``fit`` consumes condition-level moments under explicit training and
    validation roles. ``predict`` returns both the complete WitnessCell mean
    and the factorized fallback mean.  Final-target outcomes are never inputs.
    """

    def __init__(self, config: WitnessCellConfig | None = None) -> None:
        self.config = config or WitnessCellConfig()
        self._state: FittedWitnessCore | None = None

    @property
    def is_fitted(self) -> bool:
        return self._state is not None

    def __sklearn_is_fitted__(self) -> bool:
        return self.is_fitted

    @property
    def state(self) -> FittedWitnessCore:
        if self._state is None:
            raise NotFittedError("call fit before requesting fitted state")
        return self._state

    def fit(
        self,
        moments: ConditionMoments,
        split: SplitSpec,
        *,
        gene2go: Mapping[str, Sequence[Any]],
    ) -> "WitnessCell":
        """Fit endpoint gates and interaction calibration without test targets."""
        if not isinstance(moments, ConditionMoments):
            raise ValidationError("moments must be a ConditionMoments instance")
        if not isinstance(split, SplitSpec):
            raise ValidationError("split must be a SplitSpec instance")
        if not isinstance(gene2go, Mapping):
            raise ValidationError("gene2go must be a mapping")

        # Revalidate at the trust boundary so legacy or manually reconstructed
        # public objects cannot bypass the dataclass invariants.
        moments = ConditionMoments.from_mappings(
            genes=moments.genes,
            means=moments.means,
            variances=moments.variances,
            counts=moments.counts,
        )
        split = SplitSpec.create(
            split.train_conditions,
            split.validation_conditions,
        )
        required = (
            (self.config.control_label,)
            + split.train_conditions
            + split.validation_conditions
        )
        moments.require(required)
        if self.config.control_label not in split.train_conditions:
            raise ValidationError("the control condition must be in train_conditions")
        for condition in required:
            parts = parse_condition(condition)
            if condition != self.config.control_label and self.config.control_label in parts:
                raise ValidationError("control_label cannot be an endpoint of a combination")
        normalized_go: dict[str, tuple[str, ...]] = {}
        for gene, terms in gene2go.items():
            if not isinstance(gene, str) or not gene:
                raise ValidationError("gene2go keys must be non-empty strings")
            if isinstance(terms, (str, bytes)) or not isinstance(terms, Iterable):
                raise ValidationError("each gene2go value must be an iterable of terms")
            normalized_go[gene] = tuple(map(str, terms))
        self._state = fit_witness_core(
            means=moments.means,
            variances=moments.variances,
            counts=moments.counts,
            genes=moments.genes,
            train_conditions=split.train_conditions,
            validation_conditions=split.validation_conditions,
            gene2go=normalized_go,
            control_label=self.config.control_label,
            alpha_grid=self.config.alpha_grid,
            noise_grid=self.config.noise_grid,
            gate_confidence=self.config.gate_confidence,
        )
        return self

    def predict(self, conditions: Sequence[str]) -> PredictionBatch:
        """Predict condition means for singles or two-endpoint combinations."""
        conditions = tuple(map(str, conditions))
        if not conditions:
            raise ValidationError("conditions must not be empty")
        for condition in conditions:
            parts = parse_condition(condition)
            if condition != self.config.control_label and self.config.control_label in parts:
                raise ValidationError("control_label cannot be an endpoint of a combination")
        effects, factorized_effects = self.state.predict_effects(conditions)
        control = self.state.control_mean
        return PredictionBatch(
            conditions=conditions,
            genes=self.state.genes,
            means=effects + control,
            effects=effects,
            factorized_means=factorized_effects + control,
            factorized_effects=factorized_effects,
        )

    def diagnostics(self) -> dict[str, Any]:
        """Return JSON-compatible fitted gate and calibration diagnostics."""
        state = self.state
        endpoint = state.endpoint_head
        baseline = endpoint.baseline_head
        return {
            "package_contract": self.config.contract_id,
            "algorithm": "WitnessCell v14",
            "prediction_object": "condition_mean",
            "known_single_count": len(state.known_single_effects),
            "training_double_count": len(state.train_doubles),
            "validation_double_count": len(state.validation_doubles),
            "dense_active": baseline.dense_active,
            "sparse_active": baseline.sparse_active,
            "amplitude_active": endpoint.correction_active,
            "dense_all_gene_lcb": baseline.dense_all_gene_gain_lower,
            "joint_all_gene_lcb": baseline.final_all_gene_gain_lower,
            "joint_top100_pcc_lcb": baseline.final_top100_pcc_delta_lower,
            "amplitude_all_gene_lcb": endpoint.all_gene_upgrade_lower,
            "amplitude_top100_pcc_lcb": endpoint.top100_pcc_upgrade_lower,
            "alpha": state.alpha,
            "noise_ratio": state.noise_ratio,
            "gamma": state.gamma,
            "validation_mse_factorized": state.validation_mse_factorized,
            "validation_mse_witness": state.validation_mse_witness,
        }

    def save(self, path: str | Path) -> Path:
        """Save a safe, versioned model bundle."""
        from .serialization import save_model

        return save_model(self, path)

    @classmethod
    def load(cls, path: str | Path) -> "WitnessCell":
        """Load a model bundle created by :meth:`save`."""
        from .serialization import load_model

        return load_model(path)
