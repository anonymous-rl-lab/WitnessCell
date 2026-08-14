from __future__ import annotations

from zipfile import ZipFile

import numpy as np
import pytest

from witnesscell import SerializationError, WitnessCell

from ._fixture import synthetic_problem


def test_round_trip_is_exact(tmp_path) -> None:
    moments, split, gene2go, pairs = synthetic_problem()
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    path = model.save(tmp_path / "model.wcell")
    restored = WitnessCell.load(path)
    expected = model.predict(pairs[10:12])
    actual = restored.predict(pairs[10:12])
    np.testing.assert_array_equal(expected.means, actual.means)
    np.testing.assert_array_equal(expected.factorized_means, actual.factorized_means)
    assert expected.conditions == actual.conditions
    assert model.diagnostics() == restored.diagnostics()


def test_bundle_is_deterministic(tmp_path) -> None:
    moments, split, gene2go, _pairs = synthetic_problem()
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    first = model.save(tmp_path / "first.wcell")
    second = model.save(tmp_path / "second.wcell")
    assert first.read_bytes() == second.read_bytes()


def test_tampering_is_rejected(tmp_path) -> None:
    moments, split, gene2go, _pairs = synthetic_problem()
    original = WitnessCell().fit(moments, split, gene2go=gene2go).save(tmp_path / "model.wcell")
    tampered = tmp_path / "tampered.wcell"
    with ZipFile(original) as source, ZipFile(tampered, "w") as destination:
        for name in source.namelist():
            data = source.read(name)
            if name == "metadata.json":
                data += b" "
            destination.writestr(name, data)
    with pytest.raises(SerializationError, match="integrity"):
        WitnessCell.load(tampered)


def test_wrong_suffix_is_rejected(tmp_path) -> None:
    moments, split, gene2go, _pairs = synthetic_problem()
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    with pytest.raises(SerializationError, match=".wcell"):
        model.save(tmp_path / "model.zip")
