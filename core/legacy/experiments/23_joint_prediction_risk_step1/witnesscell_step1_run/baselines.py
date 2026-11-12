"""Frozen Step 1 prediction controls; no fitting or test-label loading.

``build_baselines`` accepts whitelist dictionaries for the allowed training and
validation means.  Both dictionaries may include ``control``.  It deliberately
does not accept a whole-dataset cache.  Target identities have no labels here.

The independently selected saturation-only control uses the *final corrected
v14 endpoint estimates* and all-gene MSE on validation combinations.  The saved
factorized prediction retains the final model's jointly selected alpha and is
reported separately.  No model fitting occurs in this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping, Sequence

import numpy as np


ALPHA_GRID = (0.0, 0.02, 0.05, 0.1, 0.2, 0.5)


@dataclass(frozen=True)
class BaselinePredictions:
    """Matrices follow ``conditions`` and the supplied ``genes`` exactly."""

    conditions: tuple[str, ...]
    predictions: Mapping[str, np.ndarray]
    saturation_alpha: float
    validation_scores: tuple[dict[str, float], ...]
    metadata: Mapping[str, Any]

    def __getitem__(self, method: str) -> np.ndarray:
        return self.predictions[method]


def _conditions(values: Sequence[str], name: str) -> tuple[str, ...]:
    result = tuple(str(value) for value in values)
    if len(result) != len(set(result)):
        raise ValueError(f"{name} contains duplicate conditions")
    for condition in result:
        _endpoints(condition)
    return result


def _endpoints(condition: str) -> tuple[str, ...]:
    if condition == "control":
        return ()
    result = tuple(condition.split("+"))
    if not 1 <= len(result) <= 2 or any(
        not node or node in {"control", "ctrl"} for node in result
    ):
        raise ValueError(f"Expected a clean single/pair condition: {condition!r}")
    if len(set(result)) != len(result):
        raise ValueError(f"Repeated endpoint in condition: {condition!r}")
    return result


def _vector(value: Any, n_genes: int, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (n_genes,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite vector of length {n_genes}")
    return result


def _check_whitelist(
    means: Mapping[str, np.ndarray], allowed: Sequence[str], name: str
) -> None:
    permitted = set(allowed) | {"control"}
    extra = set(means) - permitted
    missing = set(allowed) - set(means)
    if extra or missing:
        raise ValueError(
            f"{name} is not the allowed-label whitelist; "
            f"unexpected={sorted(extra)}, missing={sorted(missing)}"
        )


def endpoint_effect_from_fit(
    core_fit: Any,
    train_means: Mapping[str, np.ndarray],
    genes: Sequence[str],
    train_conditions: Sequence[str],
    gene2go: Mapping,
) -> Callable[[str], np.ndarray]:
    """Expose the final v14 corrected endpoint rule without refitting it.

    The training-single insertion order is retained, matching v4 core exactly.
    ``core_fit`` must be the final ``WitnessFitV2`` from this same allowed train
    set and gene axis; it must not be reused across leave-out boundaries.
    """
    gene_names = tuple(map(str, genes))
    train = _conditions(train_conditions, "train_conditions")
    _check_whitelist(train_means, train, "train_means")
    control = _vector(train_means["control"], len(gene_names), "control")
    known = {
        condition: _vector(train_means[condition], len(gene_names), condition) - control
        for condition in train
        if condition != "control" and "+" not in condition
    }
    if set(core_fit.known_single_genes) != set(known):
        raise ValueError("core_fit training singles do not match this whitelist")
    expected_doubles = {condition for condition in train if "+" in condition}
    if set(core_fit.training_doubles) != expected_doubles:
        raise ValueError("core_fit training doubles do not match this whitelist")
    head = core_fit.identity_head
    if not hasattr(head, "correction_active") or not hasattr(head, "baseline_head"):
        raise ValueError("Expected the final multihead_v14 identity head")
    gene_index = {gene: index for index, gene in enumerate(gene_names)}
    if len(gene_index) != len(gene_names):
        raise ValueError("genes contains duplicates")
    memo: dict[str, np.ndarray] = {}

    def endpoint(node: str) -> np.ndarray:
        if node in known:
            return known[node]
        if node not in memo:
            value = (
                head.predict(node, known, gene_index, gene2go)
                if known else np.zeros_like(control)
            )
            memo[node] = _vector(value, len(gene_names), f"endpoint {node}")
        return memo[node]

    return endpoint


def build_baselines(
    train_means: Mapping[str, np.ndarray],
    validation_means: Mapping[str, np.ndarray],
    genes: Sequence[str],
    train_conditions: Sequence[str],
    validation_conditions: Sequence[str],
    target_conditions: Sequence[str],
    *,
    core_fit: Any | None = None,
    endpoint_effect: Callable[[str], np.ndarray] | None = None,
    factorized_prediction: np.ndarray | Mapping[str, np.ndarray] | None = None,
    gene2go: Mapping | None = None,
    alpha_grid: Sequence[float] = ALPHA_GRID,
) -> BaselinePredictions:
    """Return four label-safe prediction controls in target order.

    Method keys are ``trainMean``, ``shared_endpoint_additive``,
    ``saturation_only``, and ``factorized_prediction``.  ``endpoint_effect``
    returns an effect relative to control (not an absolute expression mean).
    If omitted, it is recovered from ``core_fit`` and ``gene2go``.  The saved
    official factorized matrix may be passed in target order; otherwise the
    same core fit's factorized prediction mapping is required.

    Independent alpha selection uses all validation combinations, all genes,
    equal condition weights, and the corrected endpoint rule.  Ties prefer
    smaller alpha, matching the original core ordering.  With no validation
    combinations alpha is zero and this limitation is recorded explicitly.
    """
    gene_names = tuple(map(str, genes))
    if not gene_names or len(set(gene_names)) != len(gene_names):
        raise ValueError("genes must be nonempty and unique")
    n_genes = len(gene_names)
    train = _conditions(train_conditions, "train_conditions")
    validation = _conditions(validation_conditions, "validation_conditions")
    targets = _conditions(target_conditions, "target_conditions")
    _check_whitelist(train_means, train, "train_means")
    _check_whitelist(validation_means, validation, "validation_means")
    train_labels = set(train) - {"control"}
    validation_labels = set(validation) - {"control"}
    if train_labels & validation_labels:
        raise ValueError("Training and validation labels overlap")
    if (set(targets) - {"control"}) & (train_labels | validation_labels):
        raise ValueError("Target identities overlap allowed training/validation labels")
    control = _vector(train_means["control"], n_genes, "control")
    if "control" in validation_means and not np.array_equal(
        control, _vector(validation_means["control"], n_genes, "validation control")
    ):
        raise ValueError("Training and validation control means differ")
    known = {
        condition: _vector(train_means[condition], n_genes, condition) - control
        for condition in train
        if condition != "control" and "+" not in condition
    }
    mean_effect = (
        np.mean(np.stack(list(known.values())), axis=0)
        if known else np.zeros_like(control)
    )
    if endpoint_effect is None:
        if core_fit is None or gene2go is None:
            raise ValueError("Provide endpoint_effect, or core_fit plus gene2go")
        endpoint_effect = endpoint_effect_from_fit(
            core_fit, train_means, gene_names, train, gene2go
        )
        endpoint_source = "core_fit.identity_head.predict; final corrected v14"
    else:
        endpoint_source = "caller-supplied frozen final corrected v14 endpoint_effect"

    endpoint_cache: dict[str, np.ndarray] = {}

    def final_endpoint(node: str) -> np.ndarray:
        if node not in endpoint_cache:
            # Known endpoints remain measured training singles in the v4 core.
            endpoint_cache[node] = known[node] if node in known else _vector(
                endpoint_effect(node), n_genes, f"endpoint {node}"
            )
        return endpoint_cache[node]

    def additive(condition: str, shared: bool) -> np.ndarray:
        nodes = _endpoints(condition)
        if not nodes:
            return np.zeros_like(control)
        return np.sum(
            [final_endpoint(node) if shared else known.get(node, mean_effect)
             for node in nodes], axis=0
        )

    alphas = tuple(float(alpha) for alpha in alpha_grid)
    if not alphas or any(not np.isfinite(a) or a < 0.0 for a in alphas):
        raise ValueError("alpha_grid must contain finite nonnegative values")
    val_doubles = tuple(c for c in validation if "+" in c)
    scores: list[dict[str, float]] = []
    if val_doubles:
        val_base = np.stack([additive(c, shared=True) for c in val_doubles])
        val_truth = np.stack([
            _vector(validation_means[c], n_genes, c) - control
            for c in val_doubles
        ])
        for alpha in alphas:
            prediction = val_base / (1.0 + alpha * np.abs(val_base))
            scores.append({
                "alpha": alpha,
                "validation_all_gene_mse": float(np.mean((prediction - val_truth) ** 2)),
            })
        selected = min(scores, key=lambda row: (row["validation_all_gene_mse"], row["alpha"]))
        selected_alpha = selected["alpha"]
        selection_status = "selected_on_allowed_validation_combinations"
    else:
        selected_alpha = 0.0
        selection_status = "no_validation_combinations; fixed_alpha_zero"

    def matrix(shared: bool, alpha: float | None = None) -> np.ndarray:
        result = []
        for condition in targets:
            effect = additive(condition, shared=shared)
            if alpha is not None and "+" in condition:
                effect = effect / (1.0 + alpha * np.abs(effect))
            result.append(control + effect)
        return np.stack(result) if result else np.empty((0, n_genes), dtype=np.float64)

    if factorized_prediction is None:
        if core_fit is None:
            raise ValueError("Provide existing factorized_prediction or core_fit")
        factorized_prediction = core_fit.factorized_predictions
        factorized_source = "core_fit.factorized_predictions"
    else:
        factorized_source = "caller-supplied official factorized_prediction"
    if isinstance(factorized_prediction, Mapping):
        factorized = np.stack([
            _vector(factorized_prediction[c], n_genes, f"factorized {c}")
            for c in targets
        ]) if targets else np.empty((0, n_genes), dtype=np.float64)
    else:
        factorized = np.asarray(factorized_prediction, dtype=np.float64)
    if factorized.shape != (len(targets), n_genes) or not np.all(np.isfinite(factorized)):
        raise ValueError("factorized_prediction has invalid target/gene shape or values")

    return BaselinePredictions(
        conditions=targets,
        predictions={
            "trainMean": matrix(shared=False),
            "shared_endpoint_additive": matrix(shared=True),
            "saturation_only": matrix(shared=True, alpha=selected_alpha),
            "factorized_prediction": factorized.copy(),
        },
        saturation_alpha=selected_alpha,
        validation_scores=tuple(scores),
        metadata={
            "trainMean_definition": "measured/equal-condition-mean-completed endpoints; additive",
            "training_single_count": len(known),
            "missing_endpoint_fallback": "equal_condition_mean_effect" if known else "zero_effect_no_training_singles",
            "endpoint_source": endpoint_source,
            "saturation_selection_status": selection_status,
            "saturation_validation_metric": "all_gene_mse; equal condition and gene weights",
            "saturation_uses_final_corrected_endpoints": True,
            "saturation_alpha_grid": list(alphas),
            "saturation_validation_conditions": list(val_doubles),
            "saturation_independently_selected": bool(val_doubles),
            "factorized_source": factorized_source,
            "factorized_alpha": float(core_fit.alpha) if core_fit is not None else None,
            "factorized_is_independently_tuned_saturation": False,
            "target_truth_accessed": False,
        },
    )
