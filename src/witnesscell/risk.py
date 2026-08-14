"""Exact and estimated target-conditioned Witness Risk."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .core import EPS, parse_condition
from .exceptions import NotFittedError, ValidationError


@dataclass(frozen=True)
class ExactRiskResult:
    """Minimum-risk linear-unbiased weights and risk decomposition."""

    weights: np.ndarray
    adequacy_term: float
    geometry_term: float
    risk: float
    direct_quadratic_risk: float

    def __post_init__(self) -> None:
        weights = np.asarray(self.weights, dtype=float)
        scalars = (
            self.adequacy_term,
            self.geometry_term,
            self.risk,
            self.direct_quadratic_risk,
        )
        if weights.ndim != 1 or not np.all(np.isfinite(weights)):
            raise ValidationError("exact-risk weights must be a finite vector")
        if not np.all(np.isfinite(scalars)):
            raise ValidationError("exact-risk decomposition must be finite")
        frozen = weights.copy()
        frozen.setflags(write=False)
        object.__setattr__(self, "weights", frozen)


def exact_witness_risk(
    design: np.ndarray,
    target: np.ndarray,
    covariance: np.ndarray,
    cross_covariance: np.ndarray,
    target_variance: float,
) -> ExactRiskResult:
    """Compute the exact BLUP/universal-kriging Witness Risk formula."""
    design = np.asarray(design, dtype=float)
    target = np.asarray(target, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    cross_covariance = np.asarray(cross_covariance, dtype=float)
    if design.ndim != 2 or design.shape[0] < 1 or design.shape[1] < 1:
        raise ValidationError("design must be a non-empty two-dimensional matrix")
    rows, columns = design.shape
    if target.shape != (columns,):
        raise ValidationError("target must contain one value per design column")
    if covariance.shape != (rows, rows):
        raise ValidationError("covariance must be square over design rows")
    if cross_covariance.shape != (rows,):
        raise ValidationError("cross_covariance must contain one value per design row")
    if not np.all(np.isfinite(design)) or not np.all(np.isfinite(target)):
        raise ValidationError("design and target must contain only finite values")
    if not np.all(np.isfinite(covariance)) or not np.all(np.isfinite(cross_covariance)):
        raise ValidationError("covariance inputs must contain only finite values")
    if not np.isfinite(target_variance) or target_variance < 0:
        raise ValidationError("target_variance must be finite and non-negative")
    if not np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-12):
        raise ValidationError("covariance must be symmetric")
    if np.linalg.matrix_rank(design) != columns:
        raise ValidationError("design must have full column rank")
    joint = np.block(
        [
            [covariance, cross_covariance[:, None]],
            [cross_covariance[None, :], np.asarray([[target_variance]])],
        ]
    )
    tolerance = 1e-10 * max(1.0, float(np.linalg.norm(joint, ord=2)))
    if float(np.linalg.eigvalsh(joint).min()) < -tolerance:
        raise ValidationError("covariance and target terms must be jointly positive semidefinite")
    try:
        np.linalg.cholesky(covariance)
        inverse_design = np.linalg.solve(covariance, design)
        inverse_cross = np.linalg.solve(covariance, cross_covariance)
        information = design.T @ inverse_design
        effective_target = target - design.T @ inverse_cross
        information_solution = np.linalg.solve(information, effective_target)
    except np.linalg.LinAlgError as exc:
        raise ValidationError("covariance and design information must be nonsingular") from exc
    weights = inverse_cross + inverse_design @ information_solution
    adequacy = float(target_variance - cross_covariance @ inverse_cross)
    geometry = float(effective_target @ information_solution)
    risk = adequacy + geometry
    direct = float(
        weights @ covariance @ weights
        - 2.0 * weights @ cross_covariance
        + target_variance
    )
    return ExactRiskResult(weights, adequacy, geometry, risk, direct)


def _incidence_edges(edges: np.ndarray, n_nodes: int) -> np.ndarray:
    matrix = np.zeros((len(edges), n_nodes), dtype=float)
    matrix[np.arange(len(edges)), edges[:, 0]] = 1.0
    matrix[np.arange(len(edges)), edges[:, 1]] = 1.0
    return matrix


def _safe_edge_split(
    indices: np.ndarray,
    edges: np.ndarray,
    n_nodes: int,
    fraction: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    local_edges = edges[indices]
    degree = np.bincount(local_edges.ravel(), minlength=n_nodes)
    held_local: list[int] = []
    requested = max(1, round(fraction * len(indices)))
    for local in rng.permutation(len(indices)):
        left, right = local_edges[local]
        if degree[left] > 1 and degree[right] > 1:
            held_local.append(int(local))
            degree[left] -= 1
            degree[right] -= 1
        if len(held_local) >= requested:
            break
    if not held_local:
        raise ValidationError(
            "training edge graph cannot create an endpoint-preserving inner split"
        )
    held = np.asarray(sorted(held_local), dtype=int)
    fit = np.setdiff1d(np.arange(len(indices)), held)
    return indices[fit], indices[held]


def _pair_features(latent: np.ndarray, edges: np.ndarray) -> np.ndarray:
    upper = np.triu_indices(latent.shape[1])
    rows = []
    for left, right in edges:
        outer = 0.5 * (
            np.outer(latent[left], latent[right])
            + np.outer(latent[right], latent[left])
        )
        rows.append(
            np.r_[
                latent[left] + latent[right],
                np.abs(latent[left] - latent[right]),
                outer[upper],
            ]
        )
    return np.asarray(rows, dtype=float)


@dataclass(frozen=True)
class _Kernels:
    geometry_train: np.ndarray
    geometry_cross: np.ndarray
    discrepancy_train: np.ndarray
    discrepancy_cross: np.ndarray


def _kernel_matrices(
    train_design: np.ndarray,
    target_design: np.ndarray,
    train_features: np.ndarray,
    target_features: np.ndarray,
    length_factor: float,
) -> _Kernels:
    mean = train_features.mean(axis=0)
    scale = train_features.std(axis=0)
    scale[scale <= 1e-10] = 1.0
    train = (train_features - mean) / scale
    target = (target_features - mean) / scale
    within = cdist(train, train, metric="sqeuclidean")
    upper = within[np.triu_indices(len(train), 1)]
    reference = float(np.median(upper[upper > 0])) if np.any(upper > 0) else 1.0
    reference = max(reference, EPS)
    cross = cdist(train, target, metric="sqeuclidean")
    return _Kernels(
        geometry_train=train_design @ train_design.T / 2.0,
        geometry_cross=train_design @ target_design.T / 2.0,
        discrepancy_train=np.exp(-within / (length_factor * reference)),
        discrepancy_cross=np.exp(-cross / (length_factor * reference)),
    )


@dataclass(frozen=True)
class RiskPrediction:
    """Estimated residual mean and target-conditioned self-risk."""

    conditions: tuple[str, ...]
    residual_means: np.ndarray
    risks: np.ndarray
    weights: np.ndarray

    def __post_init__(self) -> None:
        conditions = tuple(map(str, self.conditions))
        residual = np.asarray(self.residual_means, dtype=float)
        risks = np.asarray(self.risks, dtype=float)
        weights = np.asarray(self.weights, dtype=float)
        if not conditions:
            raise ValidationError("risk prediction conditions must not be empty")
        if residual.ndim != 2 or residual.shape[0] != len(conditions):
            raise ValidationError("residual_means must be conditions × output_features")
        if risks.shape != (len(conditions),) or np.any(risks < 0):
            raise ValidationError("risks must be one non-negative value per condition")
        if weights.ndim != 2 or weights.shape[1] != len(conditions):
            raise ValidationError("weights must be training_pairs × conditions")
        if not all(np.all(np.isfinite(array)) for array in (residual, risks, weights)):
            raise ValidationError("risk predictions contain NaN or infinite values")
        object.__setattr__(self, "conditions", conditions)
        for name, array in (
            ("residual_means", residual),
            ("risks", risks),
            ("weights", weights),
        ):
            frozen = array.copy()
            frozen.setflags(write=False)
            object.__setattr__(self, name, frozen)


def _fit_predict(
    train_response: np.ndarray,
    kernels: _Kernels,
    rho: float,
    noise_ratio: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    signal_train = (
        (1.0 - rho) * kernels.geometry_train
        + rho * kernels.discrepancy_train
    )
    signal_cross = (
        (1.0 - rho) * kernels.geometry_cross
        + rho * kernels.discrepancy_cross
    )
    covariance = signal_train + noise_ratio * np.eye(len(signal_train))
    solved_cross = np.linalg.solve(covariance, signal_cross)
    solved_response = np.linalg.solve(covariance, train_response)
    scale = float(np.sum(train_response * solved_response) / train_response.size)
    scale = max(scale, EPS)
    unit_risk = 1.0 + noise_ratio - np.sum(signal_cross * solved_cross, axis=0)
    return (
        solved_cross.T @ train_response,
        scale * np.maximum(unit_risk, EPS),
        solved_cross,
        scale,
    )


@dataclass(frozen=True)
class WitnessRiskConfig:
    """Frozen Gate 07 descriptor and nested-CV hyperparameter grid."""

    pca_components: int = 12
    inner_splits: int = 3
    inner_fraction: float = 0.18
    random_state: int = 0
    length_grid: tuple[float, ...] = (0.25, 0.5, 1.0, 2.0, 4.0)
    rho_grid: tuple[float, ...] = (0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0)
    noise_grid: tuple[float, ...] = (0.03, 0.10, 0.30, 1.0, 3.0, 10.0)

    def __post_init__(self) -> None:
        if type(self.pca_components) is not int:
            raise ValidationError("pca_components must be a positive integer")
        if type(self.inner_splits) is not int:
            raise ValidationError("inner_splits must be a positive integer")
        if int(self.pca_components) < 1:
            raise ValidationError("pca_components must be a positive integer")
        if int(self.inner_splits) < 1:
            raise ValidationError("inner_splits must be a positive integer")
        if type(self.random_state) is not int:
            raise ValidationError("random_state must be a non-negative integer")
        if int(self.random_state) < 0:
            raise ValidationError("random_state must be a non-negative integer")
        try:
            fraction = float(self.inner_fraction)
            length_grid = tuple(map(float, self.length_grid))
            rho_grid = tuple(map(float, self.rho_grid))
            noise_grid = tuple(map(float, self.noise_grid))
        except (TypeError, ValueError) as exc:
            raise ValidationError("risk configuration grids must be numeric") from exc
        if not np.isfinite(fraction) or not 0.0 < fraction < 1.0:
            raise ValidationError("inner_fraction must be strictly between 0 and 1")
        if not length_grid or not np.all(np.isfinite(length_grid)) or min(length_grid) <= 0:
            raise ValidationError("length_grid must contain finite positive values")
        if not rho_grid or not np.all(np.isfinite(rho_grid)) or min(rho_grid) < 0 or max(rho_grid) > 1:
            raise ValidationError("rho_grid must contain finite values in [0, 1]")
        if not noise_grid or not np.all(np.isfinite(noise_grid)) or min(noise_grid) <= 0:
            raise ValidationError("noise_grid must contain finite positive values")
        for name, grid in (
            ("length_grid", length_grid),
            ("rho_grid", rho_grid),
            ("noise_grid", noise_grid),
        ):
            if len(grid) != len(set(grid)):
                raise ValidationError(f"{name} must not contain duplicates")
            object.__setattr__(self, name, grid)
        object.__setattr__(self, "pca_components", int(self.pca_components))
        object.__setattr__(self, "inner_splits", int(self.inner_splits))
        object.__setattr__(self, "inner_fraction", fraction)
        object.__setattr__(self, "random_state", int(self.random_state))


class WitnessRiskEstimator:
    """Training-only estimated Witness Risk for pair targets.

    The estimator is intentionally separate from :class:`WitnessCell` because
    Gate 21 validates this score as WitnessCell self-risk, not as a universal
    uncertainty model or arbitrary-model router.
    """

    def __init__(self, config: WitnessRiskConfig | None = None) -> None:
        self.config = config or WitnessRiskConfig()
        self._node_names: tuple[str, ...] | None = None
        self._latent: np.ndarray | None = None
        self._train_pairs: tuple[str, ...] | None = None
        self._train_edges: np.ndarray | None = None
        self._train_design: np.ndarray | None = None
        self._train_features: np.ndarray | None = None
        self._train_response: np.ndarray | None = None
        self._selected: dict[str, float] | None = None
        self._tuning_rows: tuple[dict[str, float], ...] = ()

    def fit(
        self,
        *,
        node_names: Sequence[str],
        single_profiles: np.ndarray,
        train_pairs: Sequence[str],
        residual_response: np.ndarray,
        reliability: Sequence[float] | None = None,
    ) -> "WitnessRiskEstimator":
        nodes = tuple(map(str, node_names))
        if any(not node.strip() for node in nodes) or len(nodes) != len(set(nodes)) or len(nodes) < 3:
            raise ValidationError("node_names must contain at least three unique nodes")
        profiles = np.asarray(single_profiles, dtype=float)
        response = np.asarray(residual_response, dtype=float)
        pairs = tuple(map(str, train_pairs))
        if profiles.ndim != 2 or profiles.shape[0] != len(nodes) or profiles.shape[1] < 1:
            raise ValidationError("single_profiles must be nodes × features")
        if response.ndim != 2 or response.shape[0] != len(pairs) or response.shape[1] < 1:
            raise ValidationError("residual_response must be pairs × output_features")
        if len(pairs) < 6:
            raise ValidationError("at least six training pairs are required")
        if not np.all(np.isfinite(profiles)) or not np.all(np.isfinite(response)):
            raise ValidationError("risk inputs contain NaN or infinite values")
        node_id = {node: index for index, node in enumerate(nodes)}
        parsed_pairs: list[tuple[str, ...]] = []
        try:
            parsed_pairs = [parse_condition(pair) for pair in pairs]
            edges = np.asarray(
                [tuple(node_id[node] for node in parts) for parts in parsed_pairs],
                dtype=int,
            )
        except KeyError as exc:
            raise ValidationError(f"pair endpoint is absent from node_names: {exc}") from exc
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
        if edges.shape != (len(pairs), 2):
            raise ValidationError("train_pairs must all be two-endpoint combinations")
        canonical_pairs = [tuple(sorted(parts)) for parts in parsed_pairs]
        if len(canonical_pairs) != len(set(canonical_pairs)):
            raise ValidationError("train_pairs contains duplicate endpoint pairs")
        weights = (
            np.ones(len(pairs), dtype=float)
            if reliability is None
            else np.asarray(reliability, dtype=float)
        )
        if (
            weights.shape != (len(pairs),)
            or not np.all(np.isfinite(weights))
            or np.any(weights < 0)
            or weights.sum() <= 0
        ):
            raise ValidationError("reliability must be finite, non-negative, and have positive sum")
        standardized = StandardScaler().fit_transform(profiles)
        components = min(self.config.pca_components, len(nodes) - 1, profiles.shape[1])
        latent = PCA(n_components=components, random_state=0).fit_transform(standardized)
        features = _pair_features(latent, edges)
        design = _incidence_edges(edges, len(nodes))
        indices = np.arange(len(edges))
        splits = [
            _safe_edge_split(
                indices,
                edges,
                len(nodes),
                self.config.inner_fraction,
                np.random.default_rng(self.config.random_state + 1009 * repeat),
            )
            for repeat in range(self.config.inner_splits)
        ]
        rows: list[dict[str, float]] = []
        for rho in self.config.rho_grid:
            local_lengths = self.config.length_grid[:1] if rho == 0.0 else self.config.length_grid
            for length_factor in local_lengths:
                for noise_ratio in self.config.noise_grid:
                    scores = []
                    for fit, validation in splits:
                        kernels = _kernel_matrices(
                            design[fit],
                            design[validation],
                            features[fit],
                            features[validation],
                            length_factor,
                        )
                        prediction, _risk, _weights, _scale = _fit_predict(
                            response[fit], kernels, rho, noise_ratio
                        )
                        per_pair = np.mean((prediction - response[validation]) ** 2, axis=1)
                        local_weight = np.maximum(weights[validation], 0.0)
                        if local_weight.sum() <= EPS:
                            local_weight = np.ones_like(local_weight)
                        scores.append(float(np.sum(per_pair * local_weight) / local_weight.sum()))
                    rows.append(
                        {
                            "length_factor": float(length_factor),
                            "rho": float(rho),
                            "noise_ratio": float(noise_ratio),
                            "inner_weighted_mse": float(np.mean(scores)),
                        }
                    )
        rows.sort(
            key=lambda row: (
                row["inner_weighted_mse"],
                row["rho"],
                row["noise_ratio"],
                row["length_factor"],
            )
        )
        self._node_names = nodes
        self._latent = latent
        self._train_pairs = pairs
        self._train_edges = edges
        self._train_design = design
        self._train_features = features
        self._train_response = response.copy()
        self._selected = dict(rows[0])
        self._tuning_rows = tuple(rows)
        return self

    def predict(self, conditions: Sequence[str]) -> RiskPrediction:
        if self._selected is None:
            raise NotFittedError("call fit before predict")
        assert self._node_names is not None
        assert self._latent is not None
        assert self._train_design is not None
        assert self._train_features is not None
        assert self._train_response is not None
        conditions = tuple(map(str, conditions))
        if not conditions:
            raise ValidationError("risk target conditions must not be empty")
        node_id = {node: index for index, node in enumerate(self._node_names)}
        try:
            target_edges = np.asarray(
                [tuple(node_id[node] for node in parse_condition(c)) for c in conditions],
                dtype=int,
            )
        except KeyError as exc:
            raise ValidationError(f"target endpoint is absent from node_names: {exc}") from exc
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
        if target_edges.shape != (len(conditions), 2):
            raise ValidationError("risk targets must be two-endpoint combinations")
        target_design = _incidence_edges(target_edges, len(self._node_names))
        target_features = _pair_features(self._latent, target_edges)
        kernels = _kernel_matrices(
            self._train_design,
            target_design,
            self._train_features,
            target_features,
            self._selected["length_factor"],
        )
        prediction, risk, weights, _scale = _fit_predict(
            self._train_response,
            kernels,
            self._selected["rho"],
            self._selected["noise_ratio"],
        )
        return RiskPrediction(conditions, prediction, risk, weights)

    def diagnostics(self) -> dict[str, object]:
        if self._selected is None:
            raise NotFittedError("call fit before diagnostics")
        return {
            "contract": "estimated-witness-gate07",
            "selected": dict(self._selected),
            "training_pairs": len(self._train_pairs or ()),
            "candidate_count": len(self._tuning_rows),
        }
