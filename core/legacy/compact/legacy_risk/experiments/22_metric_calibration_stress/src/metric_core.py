"""Source-locked metric primitives for Experiment 22.

The numerical definitions mirror the immutable upstream source recorded in
``SOURCE_LOCK.json``.  This module adds only fail-closed validation and explicit
gene alignment; it does not change the upstream arithmetic for evaluable units.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.metrics import r2_score


class MetricInputError(ValueError):
    """Raised when a formal metric unit is not safely evaluable."""


@dataclass(frozen=True)
class WeightVector:
    """A named, aligned vector of non-negative signal weights."""

    genes: np.ndarray
    values: np.ndarray

    def __post_init__(self) -> None:
        genes = np.asarray(self.genes, dtype=str)
        values = np.asarray(self.values, dtype=float)
        if genes.ndim != 1 or values.ndim != 1 or genes.size != values.size:
            raise MetricInputError("weight genes and values must be same-length vectors")
        if len(set(genes.tolist())) != genes.size:
            raise MetricInputError("aligned evaluation genes must be unique")
        if not np.all(np.isfinite(values)):
            raise MetricInputError("aligned weights contain non-finite values")
        if np.any(values < 0):
            raise MetricInputError("aligned weights contain negative values")
        if not float(values.sum()) > 0.0:
            raise MetricInputError("zero-sum weights: formal unit is non-evaluable")
        object.__setattr__(self, "genes", genes)
        object.__setattr__(self, "values", values)

    @property
    def normalized(self) -> np.ndarray:
        return self.values / np.sum(self.values)


def _as_finite_vector(name: str, value: Sequence[float]) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim != 1:
        raise MetricInputError(f"{name} must be a vector")
    if not np.all(np.isfinite(array)):
        raise MetricInputError(f"{name} contains non-finite values")
    return array


def assert_gene_order(expected: Sequence[str], observed: Sequence[str]) -> None:
    """Require exact named-gene order; silent set-based alignment is forbidden."""

    expected_array = np.asarray(expected, dtype=str)
    observed_array = np.asarray(observed, dtype=str)
    if expected_array.ndim != 1 or observed_array.ndim != 1:
        raise MetricInputError("gene identifiers must be vectors")
    if len(set(expected_array.tolist())) != expected_array.size:
        raise MetricInputError("expected evaluation genes are not unique")
    if len(set(observed_array.tolist())) != observed_array.size:
        raise MetricInputError("observed prediction genes are not unique")
    if not np.array_equal(expected_array, observed_array):
        raise MetricInputError("gene names/order mismatch; refusing positional scoring")


def source_weight_transform(
    scores: Sequence[float],
    score_genes: Sequence[str],
    evaluation_genes: Sequence[str],
) -> WeightVector:
    """Apply the locked abs→min-max→square→duplicate-max→align transform.

    This follows ``DataManager._precompute_deg_weights`` at the normative
    commit.  Duplicate score identifiers are collapsed by maximum transformed
    weight before alignment.  A constant/invalid score vector becomes a
    zero-sum vector and therefore fails closed through :class:`WeightVector`.
    """

    score_array = np.asarray(scores, dtype=float)
    score_gene_array = np.asarray(score_genes, dtype=str)
    eval_gene_array = np.asarray(evaluation_genes, dtype=str)
    if score_array.ndim != 1 or score_gene_array.ndim != 1:
        raise MetricInputError("scores and score_genes must be vectors")
    if score_array.size != score_gene_array.size:
        raise MetricInputError("scores and score_genes have different lengths")
    if eval_gene_array.ndim != 1 or len(set(eval_gene_array.tolist())) != eval_gene_array.size:
        raise MetricInputError("evaluation genes must be a unique vector")
    if score_array.size == 0:
        raise MetricInputError("empty score vector")

    absolute = np.abs(score_array)
    finite = absolute[np.isfinite(absolute)]
    if finite.size == 0:
        raise MetricInputError("all source scores are non-finite")
    min_value = np.min(absolute)
    max_value = np.max(absolute)
    with np.errstate(divide="ignore", invalid="ignore"):
        transformed = (absolute - min_value) / (max_value - min_value)
    transformed = np.nan_to_num(transformed, nan=0.0)
    transformed = np.square(transformed)

    collapsed: dict[str, float] = {}
    for gene, weight in zip(score_gene_array.tolist(), transformed.tolist(), strict=True):
        previous = collapsed.get(gene)
        if previous is None or weight > previous:
            collapsed[gene] = float(weight)
    aligned = np.asarray([collapsed.get(gene, 0.0) for gene in eval_gene_array], dtype=float)
    return WeightVector(eval_gene_array, aligned)


def wmse(prediction: Sequence[float], truth: Sequence[float], weights: WeightVector) -> float:
    """Locked weighted mean-squared error with explicit evaluability checks."""

    pred = _as_finite_vector("prediction", prediction)
    true = _as_finite_vector("truth", truth)
    if pred.shape != true.shape or pred.size != weights.values.size:
        raise MetricInputError("prediction, truth and weights have incompatible shapes")
    return float(np.sum(weights.normalized * np.square(pred - true)))


def weighted_delta_r2(
    prediction: Sequence[float],
    truth: Sequence[float],
    training_condition_mean: Sequence[float],
    weights: WeightVector,
) -> float:
    """Locked ``weighted_r2_deltapert`` relative to the training-condition mean."""

    pred = _as_finite_vector("prediction", prediction)
    true = _as_finite_vector("truth", truth)
    baseline = _as_finite_vector("training_condition_mean", training_condition_mean)
    if pred.shape != true.shape or pred.shape != baseline.shape or pred.size != weights.values.size:
        raise MetricInputError("weighted delta R2 inputs have incompatible shapes")
    value = r2_score(true - baseline, pred - baseline, sample_weight=weights.values)
    if not np.isfinite(value):
        raise MetricInputError("weighted delta R2 is non-finite")
    return float(value)


def delta_r2(
    prediction: Sequence[float], truth: Sequence[float], training_condition_mean: Sequence[float]
) -> float:
    """Unweighted context R² on perturbation-mean deltas."""

    pred = _as_finite_vector("prediction", prediction)
    true = _as_finite_vector("truth", truth)
    baseline = _as_finite_vector("training_condition_mean", training_condition_mean)
    if pred.shape != true.shape or pred.shape != baseline.shape:
        raise MetricInputError("delta R2 inputs have incompatible shapes")
    value = r2_score(true - baseline, pred - baseline)
    if not np.isfinite(value):
        raise MetricInputError("delta R2 is non-finite")
    return float(value)


def nir(
    predictions: np.ndarray,
    truths: np.ndarray,
    identities: Sequence[str],
    covariates: Sequence[str] | None = None,
) -> dict[str, float]:
    """Locked strict identity-retrieval score, evaluated within covariate.

    Ties are losses exactly as in the normative implementation.
    """

    pred = np.asarray(predictions, dtype=float)
    true = np.asarray(truths, dtype=float)
    ids = np.asarray(identities, dtype=str)
    if pred.ndim != 2 or true.ndim != 2 or pred.shape != true.shape:
        raise MetricInputError("NIR predictions and truths must be equal-shaped matrices")
    if pred.shape[0] != ids.size or len(set(ids.tolist())) != ids.size:
        raise MetricInputError("NIR identities must be unique and match matrix rows")
    if not np.all(np.isfinite(pred)) or not np.all(np.isfinite(true)):
        raise MetricInputError("NIR matrices contain non-finite values")
    if covariates is None:
        covs = np.asarray([identity.split("_")[0] for identity in ids], dtype=str)
    else:
        covs = np.asarray(covariates, dtype=str)
        if covs.shape != ids.shape:
            raise MetricInputError("NIR covariates must match identities")

    output: dict[str, float] = {}
    for covariate in dict.fromkeys(covs.tolist()):
        indices = np.flatnonzero(covs == covariate)
        if indices.size < 2:
            continue
        distances = cdist(pred[indices], true[indices], metric="euclidean")
        for local_i, global_i in enumerate(indices.tolist()):
            correct = distances[local_i, local_i]
            comparisons = [
                1.0 if correct < distances[local_i, local_j] else 0.0
                for local_j in range(indices.size)
                if local_j != local_i
            ]
            output[str(ids[global_i])] = float(np.mean(comparisons))
    return output


def drf(baseline_performance: float, duplicate_performance: float, *, higher_better: bool) -> float:
    """Locked direction-aware Dynamic Range Fraction arithmetic."""

    baseline = float(baseline_performance)
    duplicate = float(duplicate_performance)
    if not np.isfinite(baseline) or not np.isfinite(duplicate):
        return float("nan")
    perfect = 1.0 if higher_better else 0.0
    if higher_better:
        if baseline > perfect:
            raise MetricInputError("higher-is-better baseline exceeds perfect value")
        value = (duplicate - baseline) / ((perfect - baseline) + 1e-6)
    else:
        if baseline < perfect:
            raise MetricInputError("lower-is-better baseline beats perfect value")
        value = (baseline - duplicate) / (baseline + 1e-6)
    return float(np.clip(value, -1.0, 1.0))


def mse(prediction: Sequence[float], truth: Sequence[float]) -> float:
    pred = _as_finite_vector("prediction", prediction)
    true = _as_finite_vector("truth", truth)
    if pred.shape != true.shape:
        raise MetricInputError("MSE vectors have incompatible shapes")
    return float(np.mean(np.square(pred - true)))


def pearson_control_delta(
    prediction: Sequence[float], truth: Sequence[float], control: Sequence[float]
) -> float:
    pred = _as_finite_vector("prediction", prediction)
    true = _as_finite_vector("truth", truth)
    ctrl = _as_finite_vector("control", control)
    if pred.shape != true.shape or pred.shape != ctrl.shape:
        raise MetricInputError("Pearson vectors have incompatible shapes")
    value = np.corrcoef(pred - ctrl, true - ctrl)[0, 1]
    if not np.isfinite(value):
        raise MetricInputError("control-delta Pearson is non-finite")
    return float(value)

