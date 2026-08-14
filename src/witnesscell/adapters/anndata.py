"""AnnData-to-condition-moment adapter with an optional dependency."""

from __future__ import annotations

from typing import Any

import numpy as np

from ..data import ConditionMoments
from ..exceptions import ValidationError


def moments_from_anndata(
    adata: Any,
    *,
    condition_key: str,
    layer: str | None = None,
) -> ConditionMoments:
    """Aggregate cells into the condition-level moments WitnessCell consumes.

    Population variance (``ddof=0``) is used to match the frozen research
    contract.  The adapter never reads target labels beyond the supplied
    condition column and never serializes the input object.
    """
    try:
        conditions = np.asarray(adata.obs[condition_key]).astype(str)
        genes = np.asarray(adata.var_names).astype(str)
        matrix = adata.X if layer is None else adata.layers[layer]
    except (AttributeError, KeyError) as exc:
        raise ValidationError(f"invalid AnnData object or key: {exc}") from exc
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (len(conditions), len(genes)):
        raise ValidationError("AnnData expression matrix has an inconsistent shape")
    if not np.all(np.isfinite(matrix)):
        raise ValidationError("AnnData expression matrix contains NaN or infinite values")
    means: dict[str, np.ndarray] = {}
    variances: dict[str, np.ndarray] = {}
    counts: dict[str, int] = {}
    for condition in sorted(set(conditions)):
        rows = matrix[conditions == condition]
        means[condition] = rows.mean(axis=0)
        variances[condition] = rows.var(axis=0, ddof=0)
        counts[condition] = len(rows)
    return ConditionMoments.from_mappings(
        genes=tuple(genes.tolist()),
        means=means,
        variances=variances,
        counts=counts,
    )
