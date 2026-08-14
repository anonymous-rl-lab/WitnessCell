from __future__ import annotations

import json
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest

from witnesscell import SerializationError, WitnessCell, serialization

from ._fixture import synthetic_problem


def _fitted_bundle(tmp_path: Path) -> Path:
    moments, split, gene2go, _pairs = synthetic_problem()
    return WitnessCell().fit(moments, split, gene2go=gene2go).save(
        tmp_path / "model.wcell"
    )


def _rewrite_bundle(
    source: Path,
    destination: Path,
    *,
    metadata_mutator: object | None = None,
    array_mutator: object | None = None,
) -> Path:
    with ZipFile(source) as bundle:
        metadata = json.loads(bundle.read("metadata.json"))
        arrays_bytes = bundle.read("arrays.npz")
    if callable(metadata_mutator):
        metadata_mutator(metadata)
    if callable(array_mutator):
        with np.load(BytesIO(arrays_bytes), allow_pickle=False) as archive:
            arrays = {name: np.asarray(archive[name]).copy() for name in archive.files}
        array_mutator(arrays)
        buffer = BytesIO()
        np.savez_compressed(buffer, **arrays)
        arrays_bytes = buffer.getvalue()
    metadata_bytes = json.dumps(
        metadata,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    manifest_bytes = json.dumps(
        {
            "entries": {
                "arrays.npz": {
                    "sha256": sha256(arrays_bytes).hexdigest(),
                    "size": len(arrays_bytes),
                },
                "metadata.json": {
                    "sha256": sha256(metadata_bytes).hexdigest(),
                    "size": len(metadata_bytes),
                },
            },
            "format": "witnesscell-model",
            "format_version": 1,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    with ZipFile(destination, "w") as bundle:
        bundle.writestr("manifest.json", manifest_bytes)
        bundle.writestr("metadata.json", metadata_bytes)
        bundle.writestr("arrays.npz", arrays_bytes)
    return destination


def test_self_consistent_nan_array_is_rejected(tmp_path: Path) -> None:
    source = _fitted_bundle(tmp_path)

    def inject_nan(arrays: dict[str, np.ndarray]) -> None:
        arrays["control_mean"][0] = np.nan

    malicious = _rewrite_bundle(
        source,
        tmp_path / "nan.wcell",
        array_mutator=inject_nan,
    )
    with pytest.raises(SerializationError, match="control_mean"):
        WitnessCell.load(malicious)


def test_self_consistent_control_and_contract_mismatch_are_rejected(tmp_path: Path) -> None:
    source = _fitted_bundle(tmp_path)

    def change_control(metadata: dict[str, object]) -> None:
        state = metadata["state"]
        assert isinstance(state, dict)
        state["control_label"] = "other-control"

    mismatched = _rewrite_bundle(
        source,
        tmp_path / "mismatch.wcell",
        metadata_mutator=change_control,
    )
    with pytest.raises(SerializationError, match="control label"):
        WitnessCell.load(mismatched)

    def change_contract(metadata: dict[str, object]) -> None:
        metadata["algorithm_contract"] = "wrong-contract"

    wrong_contract = _rewrite_bundle(
        source,
        tmp_path / "contract.wcell",
        metadata_mutator=change_contract,
    )
    with pytest.raises(SerializationError, match="algorithm contract"):
        WitnessCell.load(wrong_contract)


def test_bundle_and_nested_archive_limits_are_enforced(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = _fitted_bundle(tmp_path)
    monkeypatch.setattr(serialization, "MAX_BUNDLE_BYTES", 1)
    with pytest.raises(SerializationError, match="safe size"):
        WitnessCell.load(source)

    monkeypatch.setattr(serialization, "MAX_BUNDLE_BYTES", 320 * 1024**2)
    monkeypatch.setattr(serialization, "MAX_TOTAL_ARRAY_ELEMENTS", 1)
    with pytest.raises(SerializationError, match="element limit"):
        WitnessCell.load(source)


def test_model_save_revalidates_mutated_state(tmp_path: Path) -> None:
    moments, split, gene2go, _pairs = synthetic_problem()
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    invalid = model.state.control_mean.copy()
    invalid[0] = np.nan
    object.__setattr__(model.state, "control_mean", invalid)
    with pytest.raises(SerializationError, match="control_mean"):
        model.save(tmp_path / "invalid.wcell")
