"""Source-locked condition-mean comparators for Experiment 22."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from metric_core import MetricInputError


def regenerate_mean_baselines(
    targets: Sequence[str],
    train_conditions: Sequence[str],
    condition_names: Sequence[str],
    condition_means: np.ndarray,
    condition_counts: np.ndarray,
    control: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return exact trainMean/baseControl condition means from released source.

    The released implementation uses an equal-condition mean of training
    singles when substituting a missing endpoint in combinations, but a
    cell-count-weighted mean across training-single cells for unseen singles.
    Preserving that asymmetry is required for source parity.
    """

    names = np.asarray(condition_names, dtype=str)
    means = np.asarray(condition_means, dtype=float)
    counts = np.asarray(condition_counts, dtype=np.int64)
    ctrl = np.asarray(control, dtype=float)
    if means.ndim != 2 or means.shape[0] != names.size or counts.shape != (names.size,):
        raise MetricInputError("comparator condition moments have incompatible shapes")
    if means.shape[1:] != ctrl.shape or not np.all(np.isfinite(means)):
        raise MetricInputError("comparator means/control are invalid")
    index = {condition: row for row, condition in enumerate(names.tolist())}
    training_singles = [
        condition
        for condition in train_conditions
        if condition != "control" and "+" not in condition
    ]
    if not training_singles or any(condition not in index for condition in training_singles):
        raise MetricInputError("training-single comparator support is empty or unresolved")
    rows = np.asarray([index[condition] for condition in training_singles], dtype=int)
    single_means = means[rows]
    single_counts = counts[rows]
    if np.any(single_counts <= 0):
        raise MetricInputError("training-single comparator support has no cells")
    equal_condition_mean = np.mean(single_means, axis=0)
    cell_weighted_mean = np.average(single_means, axis=0, weights=single_counts)
    known = {condition: means[index[condition]] for condition in training_singles}

    predictions = []
    for target in targets:
        if "+" not in target:
            predictions.append(cell_weighted_mean)
            continue
        left, right = target.split("+", maxsplit=1)
        left_known = left in known
        right_known = right in known
        if left_known and right_known:
            value = known[left] + known[right] - ctrl
        elif left_known:
            value = known[left] + equal_condition_mean - ctrl
        elif right_known:
            value = known[right] + equal_condition_mean - ctrl
        else:
            value = 2.0 * equal_condition_mean - ctrl
        predictions.append(value)
    train_mean = np.stack(predictions).astype(np.float64)
    base_control = np.repeat(ctrl[None, :], len(targets), axis=0).astype(np.float64)
    return train_mean, base_control
