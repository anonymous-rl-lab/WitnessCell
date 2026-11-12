"""Small intervention-masked model of target-RNA responses.

This module estimates *effective interventional response relations*.  It does
not identify a direct gene regulatory graph: omitted genes, intervention
strength, nonlinear response and RNA/activity differences remain assumptions.

Conventions
-----------
``U[k, j] == 1`` means that condition k directly intervenes on target j.
``Z[k, j]`` is the RAW response of target RNA j relative to a common control.
``fit_relation`` internally divides each column of Z by an explicitly supplied
positive scale.  The returned matrix has orientation ``W[j, i]`` for i -> j,
so that a normalized response row obeys ``z @ W.T ~= z`` for equations that
have not been directly intervened on.  No diagonal/self-regression is fitted.

``predict_pair_features`` takes normalized single responses (rows are the
intervened targets, columns are the target-RNA ports), in the SAME port order
and units as W.  It returns predicted-minus-additive target-port features.
No query outcomes, all-gene response labels or held-out data enter this module.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def _matrix(value: Any, name: str) -> np.ndarray:
    a = np.asarray(value, dtype=np.float64)
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError(f"{name} must be a finite two-dimensional array")
    return a


def _nonnegative_scalar(value: Any, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a nonnegative finite scalar")
    v = float(value)
    if not np.isfinite(v) or v < 0:
        raise ValueError(f"{name} must be a nonnegative finite scalar")
    return v


def _square_relation(value: Any) -> np.ndarray:
    w = _matrix(value, "W")
    if w.shape[0] != w.shape[1] or w.shape[0] < 2:
        raise ValueError("W must be square with at least two target ports")
    if np.any(np.diag(w) != 0):
        raise ValueError("W must have an exactly zero diagonal")
    return w


def response_scales(Z: Any) -> np.ndarray:
    """Return a training-only RMS normalization with no self-effect division.

    For each target port, use max(column RMS, pooled RMS).  The pooled RMS
    floor avoids expanding almost-unresponsive ports to arbitrarily large
    values.  If ALL input responses are exactly zero, return unit scales.
    Callers must provide only rows allowed for the current fit (or a common
    allowed single-intervention reference); this function does not split data.
    """
    z = _matrix(Z, "Z")
    if not z.shape[0] or not z.shape[1]:
        raise ValueError("Z must contain at least one observation and one port")
    pooled = float(np.sqrt(np.mean(z * z)))
    if not np.isfinite(pooled):
        raise ValueError("Z magnitude overflows its RMS calculation")
    if pooled == 0:
        return np.ones(z.shape[1], dtype=np.float64)
    return np.maximum(np.sqrt(np.mean(z * z, axis=0)), pooled)


def fit_relation(
    U: Any,
    Z: Any,
    scales: Any,
    ridge: float,
    mode: str = "masked",
    observation_weights: Any | None = None,
) -> dict[str, Any]:
    """Fit normalized target-RNA equations, excluding directly forced equations.

    Parameters
    ----------
    U, Z : arrays of shape (conditions, target_ports)
        U is a binary intervention indicator; Z is raw baseline-centered RNA.
        Include only single and double conditions allowed in the current fit.
    scales : positive vector of shape (target_ports,)
        Shared normalization used for both equation fitting and later queries.
    ridge : nonnegative scalar
        Each equation minimizes weighted-MEAN squared error + ridge*||w||^2.
        ``ridge=0`` uses a minimum-norm least-squares solution.
    mode : {"masked", "unmasked"}
        The unmasked comparator includes directly intervened equations but is
        otherwise identical.  It must not be called intervention-aware.
    observation_weights : optional (conditions,) or (conditions, target_ports)
        Nonnegative fitting weights.  Per-equation weights are normalized to
        sum to one after masking; zero-weight rows do not count as support.

    Returns
    -------
    dict with ``W``, ``scales`` and JSON-compatible ``audit``.  W uses normalized
    units and W[j,i] means source i, response equation j.  Unsupported equations
    are explicitly zero and identified by the audit rather than fabricated.
    """
    u = np.asarray(U)
    z = _matrix(Z, "Z")
    if u.shape != z.shape or u.ndim != 2 or not np.all((u == 0) | (u == 1)):
        raise ValueError("U must be a binary array with the same shape as Z")
    n, p = z.shape
    if n == 0 or p < 2:
        raise ValueError("Fit requires observations and at least two target ports")
    s = np.asarray(scales, dtype=np.float64)
    if s.shape != (p,) or not np.isfinite(s).all() or np.any(s <= 0):
        raise ValueError("scales must be a finite positive vector over target ports")
    penalty = _nonnegative_scalar(ridge, "ridge")
    if mode not in ("masked", "unmasked"):
        raise ValueError("mode must be masked or unmasked")
    if observation_weights is None:
        weights = np.ones((n, p), dtype=np.float64)
    else:
        weights = np.asarray(observation_weights, dtype=np.float64)
        if weights.shape == (n,):
            weights = np.broadcast_to(weights[:, None], (n, p))
        if weights.shape != (n, p) or not np.isfinite(weights).all() or np.any(weights < 0):
            raise ValueError("observation_weights must be nonnegative with shape (n,) or (n,p)")
    normalized = z / s
    if not np.isfinite(normalized).all():
        raise ValueError("Normalized responses are nonfinite")
    w = np.zeros((p, p), dtype=np.float64)
    rows = []
    for j in range(p):
        allowed = (u[:, j] == 0) if mode == "masked" else np.ones(n, dtype=bool)
        allowed = allowed & (weights[:, j] > 0)
        columns = np.arange(p) != j
        count = int(allowed.sum())
        if count == 0:
            rows.append({"target": j, "observations": 0, "weight_sum": 0.0,
                         "design_rank": 0, "weighted_fit_mse": None,
                         "coefficient_l2": 0.0, "unsupported": True})
            continue
        sample_weight = weights[allowed, j]
        weight_sum = float(sample_weight.sum())
        if not np.isfinite(weight_sum) or weight_sum <= 0:
            raise ValueError("Equation weight sum is invalid")
        root_weight = np.sqrt(sample_weight / weight_sum)
        x = normalized[allowed][:, columns]
        y = normalized[allowed, j]
        weighted_x = x * root_weight[:, None]
        weighted_y = y * root_weight
        if penalty > 0:
            design = np.vstack((weighted_x, np.sqrt(penalty) * np.eye(p - 1)))
            response = np.concatenate((weighted_y, np.zeros(p - 1)))
        else:
            design, response = weighted_x, weighted_y
        coef = np.linalg.lstsq(design, response, rcond=None)[0]
        w[j, columns] = coef
        error = x @ coef - y
        rows.append({"target": j, "observations": count, "weight_sum": weight_sum,
                     "design_rank": int(np.linalg.matrix_rank(weighted_x)),
                     "weighted_fit_mse": float(np.sum(sample_weight * error * error) / weight_sum),
                     "coefficient_l2": float(np.linalg.norm(coef)), "unsupported": False})
    if not np.isfinite(w).all():
        raise ArithmeticError("Fitted relation coefficients are nonfinite")
    return {"W": w, "scales": s.copy(), "audit": {
        "mode": mode, "ridge": penalty, "conditions": n, "target_ports": p,
        "single_conditions": int(np.sum(u.sum(axis=1) == 1)),
        "double_conditions": int(np.sum(u.sum(axis=1) == 2)),
        "other_conditions": int(np.sum((u.sum(axis=1) != 1) & (u.sum(axis=1) != 2))),
        "orientation": "W[response_port,source_port]; normalized state",
        "loss": "weighted mean equation squared error + ridge times coefficient squared norm",
        "directly_intervened_equations_excluded": mode == "masked",
        "spectral_radius_diagnostic": float(np.max(np.abs(np.linalg.eigvals(w)))),
        "spectral_radius_used_as_filter": False,
        "equations": rows,
        "interpretation": "effective interventional RNA response; not identified direct gene regulation",
    }}


def predict_pair_features(
    W: Any,
    single_ports: Any,
    edges: Any,
    solve_ridge: float = 1e-3,
) -> dict[str, Any]:
    """Predict pair-port states and return their nonadditive feature vectors.

    ``single_ports[i,j]`` MUST already equal raw single response i at port j
    divided by the fitting scale of port j.  Both axes must match W's order.

    For pair {a,b}, fix z[a]=single_ports[a,a] and z[b]=single_ports[b,b].
    On remaining ports F, minimize

      ||(I-W_FF)z_F - W_FC*z_C||^2
        + solve_ridge*||z_F - (single_ports[a]+single_ports[b])_F||^2.

    The two clamped values express an EXPLICIT transport assumption about
    target RNA, not a measured fact about a previously unseen combination.
    Augmented least squares avoids a normal-equation condition-number square.
    At solve_ridge=0, the minimum-change least-squares solution is used if the
    equations are singular.  There is no hidden spectral/stability threshold.
    AB and BA produce identical values; directed information remains inside W.
    """
    w = _square_relation(W)
    p = w.shape[0]
    singles = _matrix(single_ports, "single_ports")
    if singles.shape != (p, p):
        raise ValueError("single_ports must have one intervention and RNA column per port")
    pair_ids = np.asarray(edges)
    if pair_ids.ndim != 2 or pair_ids.shape[1] != 2 or pair_ids.dtype.kind not in "iu":
        raise ValueError("edges must be an integer array of shape (pairs,2)")
    if np.any(pair_ids < 0) or np.any(pair_ids >= p) or np.any(pair_ids[:, 0] == pair_ids[:, 1]):
        raise ValueError("edges must contain distinct target indices within range")
    penalty = _nonnegative_scalar(solve_ridge, "solve_ridge")
    # Canonical order removes even floating-point ordering differences for AB/BA.
    pair_ids = np.sort(pair_ids.astype(np.int64, copy=False), axis=1)
    additive = singles[pair_ids[:, 0]] + singles[pair_ids[:, 1]]
    prediction = additive.copy()
    records = []
    for k, fixed in enumerate(pair_ids):
        free = np.setdiff1d(np.arange(p), fixed, assume_unique=True)
        clamped = singles[fixed, fixed]
        prediction[k, fixed] = clamped
        if len(free):
            equation = np.eye(len(free)) - w[np.ix_(free, free)]
            rhs = w[np.ix_(free, fixed)] @ clamped
            residual_rhs = rhs - equation @ additive[k, free]
            if penalty > 0:
                design = np.vstack((equation, np.sqrt(penalty) * np.eye(len(free))))
                response = np.concatenate((residual_rhs, np.zeros(len(free))))
            else:
                design, response = equation, residual_rhs
            change = np.linalg.lstsq(design, response, rcond=None)[0]
            prediction[k, free] = additive[k, free] + change
            condition = float(np.linalg.cond(equation))
            equation_rms = float(np.sqrt(np.mean((equation @ prediction[k, free] - rhs) ** 2)))
            rank = int(np.linalg.matrix_rank(equation))
        else:
            condition, equation_rms, rank = 1.0, 0.0, 0
        records.append({"targets": fixed.tolist(), "free_equation_rank": rank,
                        "free_equation_condition": condition if np.isfinite(condition) else None,
                        "free_equation_condition_infinite": not np.isfinite(condition),
                        "free_equation_residual_rms": equation_rms})
    features = prediction - additive
    if not np.isfinite(prediction).all() or not np.isfinite(features).all():
        raise ArithmeticError("Pair response solve produced nonfinite output")
    return {"features": features, "prediction": prediction, "additive": additive,
            "audit": {"pair_count": len(pair_ids), "target_ports": p,
                      "solve_ridge": penalty, "feature_dimension": p,
                      "feature_definition": "normalized predicted target-port response minus additive single response",
                      "clamp_assumption": "each directly intervened port retains its own measured single-intervention RNA displacement",
                      "query_response_labels_used": False, "pairs": records}}


def permute_relation(W: Any, permutation: Any) -> np.ndarray:
    """Mismatch relation identities while preserving its graph/weight multiset.

    The CALLER must keep intervention identities and single-response rows and
    columns fixed.  Reordering those too would be a mere relabeling, not a
    negative control.  This permutation preserves graph shape and the weight
    multiset, but generally not each named gene's incoming/outgoing degree.
    """
    w = _square_relation(W)
    ids = np.asarray(permutation)
    p = len(w)
    if ids.shape != (p,) or ids.dtype.kind not in "iu" or not np.array_equal(np.sort(ids), np.arange(p)):
        raise ValueError("permutation must be an integer bijection over target ports")
    return w[np.ix_(ids, ids)].copy()
