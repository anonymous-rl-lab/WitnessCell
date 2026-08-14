from __future__ import annotations

import numpy as np

from witnesscell import ConditionMoments
from witnesscell.adapters import moments_from_anndata

from ._fixture import synthetic_problem


def test_moment_npz_round_trip(tmp_path) -> None:
    moments, _split, _gene2go, _pairs = synthetic_problem()
    path = tmp_path / "moments.npz"
    moments.to_npz(path)
    restored = ConditionMoments.from_npz(path)
    assert restored.genes == moments.genes
    assert restored.counts == moments.counts
    for condition in moments.means:
        np.testing.assert_array_equal(restored.means[condition], moments.means[condition])


def test_anndata_adapter_without_import_dependency() -> None:
    class Dummy:
        X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
        obs = {"condition": np.array(["control", "control", "A"])}
        var_names = np.array(["A", "B"])
        layers: dict[str, np.ndarray] = {}

    moments = moments_from_anndata(Dummy(), condition_key="condition")
    np.testing.assert_array_equal(moments.means["control"], [2.0, 3.0])
    np.testing.assert_array_equal(moments.variances["control"], [1.0, 1.0])
    assert moments.counts == {"A": 1, "control": 2}
