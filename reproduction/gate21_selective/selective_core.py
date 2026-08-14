#!/usr/bin/env python3
"""Shared calculations for the frozen selective-prediction gate."""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression


EPS = 1e-12


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_bucket(text: str, salt: str, modulus: int = 10) -> int:
    payload = f"{salt}::{text}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big") % modulus


def pair_balanced_weights(frame: pd.DataFrame) -> np.ndarray:
    counts = frame.groupby("pair")["pair"].transform("size").to_numpy(float)
    return 1.0 / counts


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    return float(np.sum(values * weights) / np.sum(weights))


def weighted_quantile(values: np.ndarray, quantile: float, weights: np.ndarray) -> float:
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights) / np.sum(weights)
    index = int(np.searchsorted(cdf, quantile, side="left"))
    return float(values[min(index, len(values) - 1)])


def fit_calibration(
    frame: pd.DataFrame,
    score: str,
    loss: str,
    coverage_grid: Iterable[float],
) -> tuple[dict, pd.DataFrame]:
    weights = pair_balanced_weights(frame)
    x = frame[score].to_numpy(float)
    y = frame[loss].to_numpy(float)
    model = IsotonicRegression(increasing=True, out_of_bounds="clip", y_min=0.0)
    model.fit(x, y, sample_weight=weights)
    rows = []
    for coverage in coverage_grid:
        threshold = weighted_quantile(x, float(coverage), weights)
        accepted = x <= threshold
        rows.append({
            "target_coverage": float(coverage),
            "risk_threshold": threshold,
            "calibrated_loss_threshold": float(model.predict([threshold])[0]),
            "calibration_pair_balanced_coverage": weighted_mean(
                accepted.astype(float), weights
            ),
            "calibration_selective_mse": weighted_mean(y[accepted], weights[accepted]),
        })
    mapping = {
        "score": score,
        "loss": loss,
        "x_thresholds": [float(v) for v in model.X_thresholds_],
        "y_thresholds": [float(v) for v in model.y_thresholds_],
        "calibration_rows": int(len(frame)),
        "calibration_pairs": int(frame.pair.nunique()),
    }
    return mapping, pd.DataFrame(rows)


def evaluate_threshold(
    frame: pd.DataFrame,
    score: str,
    loss: str,
    threshold: float,
) -> dict:
    weights = pair_balanced_weights(frame)
    risk = frame[score].to_numpy(float)
    error = frame[loss].to_numpy(float)
    accepted = risk <= threshold
    rejected = ~accepted
    if not accepted.any() or not rejected.any():
        raise ValueError("threshold produced an empty accepted or rejected set")
    all_mse = weighted_mean(error, weights)
    accepted_mse = weighted_mean(error[accepted], weights[accepted])
    rejected_mse = weighted_mean(error[rejected], weights[rejected])
    return {
        "rows": int(len(frame)),
        "pairs": int(frame.pair.nunique()),
        "accepted_rows": int(accepted.sum()),
        "rejected_rows": int(rejected.sum()),
        "row_coverage": float(accepted.mean()),
        "pair_balanced_coverage": weighted_mean(accepted.astype(float), weights),
        "all_mse": all_mse,
        "accepted_mse": accepted_mse,
        "rejected_mse": rejected_mse,
        "accepted_over_all_mse": accepted_mse / max(all_mse, EPS),
        "accepted_over_rejected_mse": accepted_mse / max(rejected_mse, EPS),
        "score_loss_spearman": float(spearmanr(risk, error).statistic),
    }


def cluster_bootstrap_ratio(
    frame: pd.DataFrame,
    score: str,
    loss: str,
    threshold: float,
    replicates: int,
    seed: int,
) -> np.ndarray:
    groups = [local for _, local in frame.groupby("pair", sort=True)]
    all_mean = []
    accepted_numerator = []
    accepted_denominator = []
    for local in groups:
        error = local[loss].to_numpy(float)
        accepted = local[score].to_numpy(float) <= threshold
        n_rows = len(local)
        all_mean.append(float(np.mean(error)))
        accepted_numerator.append(float(np.sum(error[accepted]) / n_rows))
        accepted_denominator.append(float(np.sum(accepted) / n_rows))
    all_mean = np.asarray(all_mean, float)
    accepted_numerator = np.asarray(accepted_numerator, float)
    accepted_denominator = np.asarray(accepted_denominator, float)
    rng = np.random.default_rng(seed)
    sampled = rng.integers(0, len(groups), size=(replicates, len(groups)))
    accepted_den = accepted_denominator[sampled].sum(axis=1)
    valid = accepted_den > 0
    accepted_mse = accepted_numerator[sampled].sum(axis=1)[valid] / accepted_den[valid]
    overall_mse = all_mean[sampled].mean(axis=1)[valid]
    return accepted_mse / np.maximum(overall_mse, EPS)


def within_seed_permutation(
    frame: pd.DataFrame,
    score: str,
    loss: str,
    threshold: float,
    replicates: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    weights = pair_balanced_weights(frame)
    error = frame[loss].to_numpy(float)
    risk = frame[score].to_numpy(float)
    numerator = np.zeros(replicates, float)
    denominator = np.zeros(replicates, float)
    for _, indices in frame.groupby("seed").groups.items():
        indices = np.asarray(indices, int)
        accepted_count = int(np.sum(risk[indices] <= threshold))
        if accepted_count == 0:
            continue
        random_order = rng.random((replicates, len(indices)))
        selected_local = np.argpartition(
            random_order, accepted_count - 1, axis=1
        )[:, :accepted_count]
        selected = indices[selected_local]
        numerator += np.sum(weights[selected] * error[selected], axis=1)
        denominator += np.sum(weights[selected], axis=1)
    valid = denominator > 0
    return numerator[valid] / denominator[valid]
