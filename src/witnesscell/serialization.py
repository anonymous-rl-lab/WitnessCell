"""Deterministic, pickle-free WitnessCell model bundles."""

from __future__ import annotations

import json
import os
from collections.abc import Mapping
from dataclasses import asdict
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from tempfile import NamedTemporaryFile
from types import MappingProxyType
from typing import Any, cast
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import numpy as np

from ._version import __version__
from .core import (
    DenseSparseEndpointHead,
    FittedWitnessCore,
    IncrementalAmplitudeHead,
    parse_condition,
)
from .exceptions import SerializationError
from .model import WitnessCell, WitnessCellConfig

MODEL_FORMAT = "witnesscell-model"
MODEL_FORMAT_VERSION = 1
MAX_BUNDLE_BYTES = 320 * 1024**2
MAX_MANIFEST_BYTES = 64 * 1024
MAX_METADATA_BYTES = 16 * 1024**2
MAX_ARRAY_ARCHIVE_BYTES = 256 * 1024**2
MAX_NPZ_MEMBER_BYTES = 256 * 1024**2
MAX_TOTAL_NPZ_BYTES = 512 * 1024**2
MAX_TOTAL_ARRAY_ELEMENTS = 64_000_000
_EXPECTED_ENTRIES = frozenset({"manifest.json", "metadata.json", "arrays.npz"})
_EXPECTED_ARRAYS = frozenset({"control_mean", "known_single_effects", "train_residual"})
_ENTRY_LIMITS = {
    "manifest.json": MAX_MANIFEST_BYTES,
    "metadata.json": MAX_METADATA_BYTES,
    "arrays.npz": MAX_ARRAY_ARCHIVE_BYTES,
}


def _exact_keys(payload: Any, expected: set[str] | frozenset[str], label: str) -> dict[str, Any]:
    if not isinstance(payload, dict) or set(payload) != set(expected):
        raise SerializationError(f"{label} has missing or unexpected fields")
    return payload


def _strict_bool(value: Any, label: str) -> bool:
    if type(value) is not bool:
        raise SerializationError(f"{label} must be a boolean")
    return value


def _strict_int(value: Any, label: str, *, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise SerializationError(f"{label} must be an integer >= {minimum}")
    return value


def _float_tuple(value: Any, width: int, label: str) -> tuple[float, ...]:
    if not isinstance(value, (list, tuple)) or len(value) != width:
        raise SerializationError(f"{label} must contain exactly {width} values")
    if any(isinstance(item, bool) for item in value):
        raise SerializationError(f"{label} must be numeric")
    try:
        result = tuple(map(float, value))
    except (TypeError, ValueError) as exc:
        raise SerializationError(f"{label} must be numeric") from exc
    if not np.all(np.isfinite(result)):
        raise SerializationError(f"{label} must contain only finite values")
    return result


def _metric(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise SerializationError(f"{label} must be numeric")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise SerializationError(f"{label} must be numeric") from exc
    if np.isnan(result) or np.isposinf(result):
        raise SerializationError(f"{label} must be finite or negative infinity")
    return result


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise SerializationError(f"{label} must be numeric")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise SerializationError(f"{label} must be numeric") from exc
    if not np.isfinite(result):
        raise SerializationError(f"{label} must be finite")
    return result


def _json_safe(value: Any) -> Any:
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if np.isnan(number):
            return {"__float__": "nan"}
        if np.isposinf(number):
            return {"__float__": "+inf"}
        if np.isneginf(number):
            return {"__float__": "-inf"}
        return number
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_safe(item) for item in value]
    return value


def _json_restore(value: Any) -> Any:
    if isinstance(value, dict) and set(value) == {"__float__"}:
        encoded = value["__float__"]
        return {"nan": float("nan"), "+inf": float("inf"), "-inf": float("-inf")}[encoded]
    if isinstance(value, dict):
        return {key: _json_restore(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_restore(item) for item in value]
    return value


def _json_bytes(value: Any) -> bytes:
    return json.dumps(
        _json_safe(value), ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _zip_info(name: str) -> ZipInfo:
    info = ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
    info.compress_type = ZIP_DEFLATED
    info.external_attr = 0o100644 << 16
    return info


def _head_metadata(head: IncrementalAmplitudeHead) -> dict[str, Any]:
    baseline = head.baseline_head
    return {
        "baseline": asdict(baseline),
        "correction_active": head.correction_active,
        "correction_weights": head.correction_weights,
        "all_gene_upgrade_mean": head.all_gene_upgrade_mean,
        "all_gene_upgrade_lower": head.all_gene_upgrade_lower,
        "top100_mse_upgrade_mean": head.top100_mse_upgrade_mean,
        "top100_mse_upgrade_lower": head.top100_mse_upgrade_lower,
        "top100_pcc_upgrade_mean": head.top100_pcc_upgrade_mean,
        "top100_pcc_upgrade_lower": head.top100_pcc_upgrade_lower,
        "known_single_count": head.known_single_count,
        "records": head.records,
        "ridge": head.ridge,
    }


def _restore_baseline(payload: dict[str, Any]) -> DenseSparseEndpointHead:
    _exact_keys(
        payload,
        {
            "dense_active",
            "sparse_active",
            "dense_weights",
            "joint_weights",
            "dense_all_gene_gain_mean",
            "dense_all_gene_gain_lower",
            "final_all_gene_gain_mean",
            "final_all_gene_gain_lower",
            "final_top100_pcc_delta_mean",
            "final_top100_pcc_delta_lower",
            "known_single_count",
            "records",
            "go_top_k",
        },
        "endpoint baseline",
    )
    if not isinstance(payload["records"], list) or not all(
        isinstance(record, dict) for record in payload["records"]
    ):
        raise SerializationError("endpoint baseline records must be JSON objects")
    return DenseSparseEndpointHead(
        dense_active=_strict_bool(payload["dense_active"], "dense_active"),
        sparse_active=_strict_bool(payload["sparse_active"], "sparse_active"),
        dense_weights=cast(
            tuple[float, float],
            _float_tuple(payload["dense_weights"], 2, "dense_weights"),
        ),
        joint_weights=cast(
            tuple[float, float, float],
            _float_tuple(payload["joint_weights"], 3, "joint_weights"),
        ),
        dense_all_gene_gain_mean=_metric(
            payload["dense_all_gene_gain_mean"], "dense_all_gene_gain_mean"
        ),
        dense_all_gene_gain_lower=_metric(
            payload["dense_all_gene_gain_lower"], "dense_all_gene_gain_lower"
        ),
        final_all_gene_gain_mean=_metric(
            payload["final_all_gene_gain_mean"], "final_all_gene_gain_mean"
        ),
        final_all_gene_gain_lower=_metric(
            payload["final_all_gene_gain_lower"], "final_all_gene_gain_lower"
        ),
        final_top100_pcc_delta_mean=_metric(
            payload["final_top100_pcc_delta_mean"], "final_top100_pcc_delta_mean"
        ),
        final_top100_pcc_delta_lower=_metric(
            payload["final_top100_pcc_delta_lower"], "final_top100_pcc_delta_lower"
        ),
        known_single_count=_strict_int(
            payload["known_single_count"], "baseline known_single_count"
        ),
        records=tuple(dict(record) for record in payload["records"]),
        go_top_k=_strict_int(payload["go_top_k"], "go_top_k", minimum=1),
    )


def _restore_head(payload: dict[str, Any]) -> IncrementalAmplitudeHead:
    _exact_keys(
        payload,
        {
            "baseline",
            "correction_active",
            "correction_weights",
            "all_gene_upgrade_mean",
            "all_gene_upgrade_lower",
            "top100_mse_upgrade_mean",
            "top100_mse_upgrade_lower",
            "top100_pcc_upgrade_mean",
            "top100_pcc_upgrade_lower",
            "known_single_count",
            "records",
            "ridge",
        },
        "endpoint head",
    )
    if not isinstance(payload["baseline"], dict):
        raise SerializationError("endpoint baseline must be a JSON object")
    if not isinstance(payload["records"], list) or not all(
        isinstance(record, dict) for record in payload["records"]
    ):
        raise SerializationError("endpoint head records must be JSON objects")
    ridge = float(payload["ridge"])
    if not np.isfinite(ridge) or ridge <= 0:
        raise SerializationError("endpoint ridge must be finite and positive")
    return IncrementalAmplitudeHead(
        baseline_head=_restore_baseline(payload["baseline"]),
        correction_active=_strict_bool(payload["correction_active"], "correction_active"),
        correction_weights=cast(
            tuple[float, float, float],
            _float_tuple(payload["correction_weights"], 3, "correction_weights"),
        ),
        all_gene_upgrade_mean=_metric(
            payload["all_gene_upgrade_mean"], "all_gene_upgrade_mean"
        ),
        all_gene_upgrade_lower=_metric(
            payload["all_gene_upgrade_lower"], "all_gene_upgrade_lower"
        ),
        top100_mse_upgrade_mean=_metric(
            payload["top100_mse_upgrade_mean"], "top100_mse_upgrade_mean"
        ),
        top100_mse_upgrade_lower=_metric(
            payload["top100_mse_upgrade_lower"], "top100_mse_upgrade_lower"
        ),
        top100_pcc_upgrade_mean=_metric(
            payload["top100_pcc_upgrade_mean"], "top100_pcc_upgrade_mean"
        ),
        top100_pcc_upgrade_lower=_metric(
            payload["top100_pcc_upgrade_lower"], "top100_pcc_upgrade_lower"
        ),
        known_single_count=_strict_int(
            payload["known_single_count"], "head known_single_count"
        ),
        records=tuple(dict(record) for record in payload["records"]),
        ridge=ridge,
    )


def _validate_double_conditions(
    conditions: tuple[str, ...], control_label: str, label: str
) -> None:
    if len(conditions) != len(set(conditions)):
        raise SerializationError(f"{label} contains duplicates")
    for condition in conditions:
        try:
            parts = parse_condition(condition)
        except ValueError as exc:
            raise SerializationError(f"invalid {label}: {exc}") from exc
        if len(parts) != 2:
            raise SerializationError(f"{label} must contain only two-endpoint combinations")
        if control_label in parts:
            raise SerializationError(f"{label} cannot contain the control label")


def _validate_model_state(
    config: WitnessCellConfig,
    state: FittedWitnessCore,
    algorithm_contract: str,
) -> None:
    if algorithm_contract != config.contract_id:
        raise SerializationError("algorithm contract does not match model configuration")
    if state.control_label != config.control_label:
        raise SerializationError("state control label does not match model configuration")
    if not state.genes or any(not gene.strip() for gene in state.genes):
        raise SerializationError("model genes must contain non-empty names")
    if len(state.genes) != len(set(state.genes)):
        raise SerializationError("model genes must be unique")
    width = len(state.genes)
    control = np.asarray(state.control_mean, dtype=float)
    if control.shape != (width,) or not np.all(np.isfinite(control)):
        raise SerializationError("control_mean must be a finite gene vector")

    known_nodes = tuple(state.known_single_effects)
    if len(known_nodes) != len(set(known_nodes)):
        raise SerializationError("known single nodes must be unique")
    for node, effect in state.known_single_effects.items():
        try:
            parts = parse_condition(node)
        except ValueError as exc:
            raise SerializationError(f"invalid known single node: {exc}") from exc
        if len(parts) != 1 or node == state.control_label:
            raise SerializationError("known single nodes must be non-control singles")
        array = np.asarray(effect, dtype=float)
        if array.shape != (width,) or not np.all(np.isfinite(array)):
            raise SerializationError("known single effects must be finite gene vectors")

    if not isinstance(state.gene2go, Mapping):
        raise SerializationError("gene2go state must be a mapping")
    for gene, terms in state.gene2go.items():
        if not isinstance(gene, str) or not gene:
            raise SerializationError("gene2go keys must be non-empty strings")
        if not isinstance(terms, tuple) or not all(isinstance(term, str) for term in terms):
            raise SerializationError("gene2go values must be tuples of strings")

    _validate_double_conditions(state.train_doubles, state.control_label, "train_doubles")
    _validate_double_conditions(
        state.validation_doubles, state.control_label, "validation_doubles"
    )
    overlap = set(state.train_doubles) & set(state.validation_doubles)
    if overlap:
        raise SerializationError("training and validation doubles overlap")
    residual = np.asarray(state.train_residual, dtype=float)
    if residual.shape != (len(state.train_doubles), width) or not np.all(
        np.isfinite(residual)
    ):
        raise SerializationError("train_residual must be a finite train_doubles × genes matrix")

    if not np.isfinite(state.alpha) or state.alpha < 0:
        raise SerializationError("alpha must be finite and non-negative")
    if not np.isfinite(state.noise_ratio) or state.noise_ratio < 0:
        raise SerializationError("noise_ratio must be finite and non-negative")
    if not np.isfinite(state.gamma) or not 0.0 <= state.gamma <= 1.0:
        raise SerializationError("gamma must be finite and in [0, 1]")
    if state.train_doubles and state.noise_ratio <= 0:
        raise SerializationError("models with training doubles require positive noise_ratio")
    if not state.train_doubles and (state.noise_ratio != 0 or state.gamma != 0):
        raise SerializationError("models without training doubles require zero noise and gamma")

    validation_metrics = (
        float(state.validation_mse_factorized),
        float(state.validation_mse_witness),
    )
    if state.validation_doubles:
        if not np.all(np.isfinite(validation_metrics)) or min(validation_metrics) < 0:
            raise SerializationError("validation MSE values must be finite and non-negative")
    elif not all(np.isnan(value) for value in validation_metrics):
        raise SerializationError("models without validation doubles require NaN validation MSE")

    head = state.endpoint_head
    baseline = head.baseline_head
    if baseline.sparse_active and not baseline.dense_active:
        raise SerializationError("sparse endpoint gate cannot be active without dense gate")
    if baseline.known_single_count != head.known_single_count:
        raise SerializationError("endpoint heads disagree on known_single_count")
    if not 0 <= head.known_single_count <= len(known_nodes):
        raise SerializationError("endpoint known_single_count is inconsistent")
    if baseline.go_top_k < 1 or not np.isfinite(head.ridge) or head.ridge <= 0:
        raise SerializationError("endpoint head hyperparameters are invalid")
    dense_weights = np.asarray(baseline.dense_weights, dtype=float)
    joint_weights = np.asarray(baseline.joint_weights, dtype=float)
    correction_weights = np.asarray(head.correction_weights, dtype=float)
    if dense_weights.shape != (2,) or joint_weights.shape != (3,):
        raise SerializationError("endpoint baseline weights have invalid dimensions")
    if correction_weights.shape != (3,):
        raise SerializationError("endpoint correction weights have invalid dimensions")
    if not all(
        np.all(np.isfinite(weights))
        for weights in (dense_weights, joint_weights, correction_weights)
    ):
        raise SerializationError("endpoint weights must be finite")
    if np.any(dense_weights < 0) or np.any(dense_weights > 2):
        raise SerializationError("dense endpoint weights are outside the fitted range")
    if np.any(joint_weights < 0) or np.any(joint_weights > 2):
        raise SerializationError("joint endpoint weights are outside the fitted range")
    if np.any(correction_weights < -1.5) or np.any(correction_weights > 1.5):
        raise SerializationError("amplitude weights are outside the fitted range")


def save_model(model: WitnessCell, path: str | Path) -> Path:
    """Atomically save a fitted estimator as a deterministic ``.wcell`` ZIP."""
    state = model.state
    _validate_model_state(model.config, state, model.config.contract_id)
    destination = Path(path)
    if destination.suffix != ".wcell":
        raise SerializationError("WitnessCell model paths must end in .wcell")
    destination.parent.mkdir(parents=True, exist_ok=True)

    nodes = tuple(sorted(state.known_single_effects))
    arrays_buffer = BytesIO()
    np.savez_compressed(
        arrays_buffer,
        control_mean=np.asarray(state.control_mean, dtype=np.float64),
        known_single_effects=(
            np.stack([state.known_single_effects[node] for node in nodes]).astype(np.float64)
            if nodes
            else np.zeros((0, len(state.genes)), dtype=np.float64)
        ),
        train_residual=np.asarray(state.train_residual, dtype=np.float64),
    )
    arrays_bytes = arrays_buffer.getvalue()
    metadata = {
        "format": MODEL_FORMAT,
        "format_version": MODEL_FORMAT_VERSION,
        "package_version": __version__,
        "algorithm_contract": model.config.contract_id,
        "config": asdict(model.config),
        "state": {
            "genes": state.genes,
            "control_label": state.control_label,
            "known_single_nodes": nodes,
            "gene2go": state.gene2go,
            "endpoint_head": _head_metadata(state.endpoint_head),
            "alpha": state.alpha,
            "noise_ratio": state.noise_ratio,
            "gamma": state.gamma,
            "train_doubles": state.train_doubles,
            "validation_doubles": state.validation_doubles,
            "validation_mse_factorized": state.validation_mse_factorized,
            "validation_mse_witness": state.validation_mse_witness,
        },
    }
    metadata_bytes = _json_bytes(metadata)
    manifest_bytes = _json_bytes(
        {
            "format": MODEL_FORMAT,
            "format_version": MODEL_FORMAT_VERSION,
            "entries": {
                "metadata.json": {"sha256": sha256(metadata_bytes).hexdigest(), "size": len(metadata_bytes)},
                "arrays.npz": {"sha256": sha256(arrays_bytes).hexdigest(), "size": len(arrays_bytes)},
            },
        }
    )
    try:
        with NamedTemporaryFile(dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
        with ZipFile(temporary, "w") as bundle:
            bundle.writestr(_zip_info("manifest.json"), manifest_bytes)
            bundle.writestr(_zip_info("metadata.json"), metadata_bytes)
            bundle.writestr(_zip_info("arrays.npz"), arrays_bytes)
        os.replace(temporary, destination)
    except Exception as exc:
        if "temporary" in locals():
            temporary.unlink(missing_ok=True)
        raise SerializationError(f"could not save model bundle: {exc}") from exc
    return destination


def _read_limited(bundle: ZipFile, name: str, limit: int) -> bytes:
    info = bundle.getinfo(name)
    if info.file_size > limit:
        raise SerializationError(f"model entry is too large: {name}")
    with bundle.open(info, "r") as handle:
        payload = handle.read(limit + 1)
    if len(payload) > limit or len(payload) != info.file_size:
        raise SerializationError(f"model entry size is invalid: {name}")
    return payload


def _json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise SerializationError(f"JSON object contains duplicate key: {key}")
        result[key] = value
    return result


def _load_json(payload: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(payload, object_pairs_hook=_json_object)
    except SerializationError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SerializationError(f"{label} is not valid UTF-8 JSON") from exc
    if not isinstance(value, dict):
        raise SerializationError(f"{label} must contain a JSON object")
    return cast(dict[str, Any], _json_restore(value))


def _load_array_archive(arrays_bytes: bytes) -> dict[str, np.ndarray]:
    try:
        with ZipFile(BytesIO(arrays_bytes), "r") as archive_zip:
            infos = archive_zip.infolist()
            names = [info.filename for info in infos]
            expected_names = {f"{name}.npy" for name in _EXPECTED_ARRAYS}
            if len(names) != len(set(names)) or set(names) != expected_names:
                raise SerializationError("model array archive has unexpected entries")
            total = 0
            for info in infos:
                if info.file_size > MAX_NPZ_MEMBER_BYTES:
                    raise SerializationError(
                        f"model array member is too large: {info.filename}"
                    )
                total += info.file_size
            if total > MAX_TOTAL_NPZ_BYTES:
                raise SerializationError("model array archive expands beyond the safe limit")

        result: dict[str, np.ndarray] = {}
        total_elements = 0
        with np.load(BytesIO(arrays_bytes), allow_pickle=False) as archive:
            if set(archive.files) != set(_EXPECTED_ARRAYS):
                raise SerializationError("model array archive has unexpected entries")
            for name in sorted(_EXPECTED_ARRAYS):
                raw = archive[name]
                if raw.dtype.kind != "f" or raw.dtype.itemsize != 8:
                    raise SerializationError(f"model array must use float64: {name}")
                total_elements += int(raw.size)
                if total_elements > MAX_TOTAL_ARRAY_ELEMENTS:
                    raise SerializationError("model arrays exceed the safe element limit")
                array = np.asarray(raw, dtype=np.float64).copy()
                array.setflags(write=False)
                result[name] = array
        return result
    except SerializationError:
        raise
    except Exception as exc:
        raise SerializationError(f"invalid model array archive: {exc}") from exc


def _read_bundle(path: Path) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    try:
        if path.stat().st_size > MAX_BUNDLE_BYTES:
            raise SerializationError("model bundle exceeds the safe size limit")
        with ZipFile(path, "r") as bundle:
            names = bundle.namelist()
            if len(names) != len(set(names)) or frozenset(names) != _EXPECTED_ENTRIES:
                raise SerializationError("model bundle has unexpected or duplicate entries")
            manifest_bytes = _read_limited(
                bundle, "manifest.json", _ENTRY_LIMITS["manifest.json"]
            )
            metadata_bytes = _read_limited(
                bundle, "metadata.json", _ENTRY_LIMITS["metadata.json"]
            )
            arrays_bytes = _read_limited(
                bundle, "arrays.npz", _ENTRY_LIMITS["arrays.npz"]
            )
    except SerializationError:
        raise
    except Exception as exc:
        raise SerializationError(f"could not read model bundle: {exc}") from exc

    try:
        manifest = _load_json(manifest_bytes, "manifest.json")
        _exact_keys(manifest, {"format", "format_version", "entries"}, "manifest")
        if (
            manifest["format"] != MODEL_FORMAT
            or type(manifest["format_version"]) is not int
            or manifest["format_version"] != MODEL_FORMAT_VERSION
        ):
            raise SerializationError("unsupported model format or format version")
        expected = manifest["entries"]
        _exact_keys(expected, {"metadata.json", "arrays.npz"}, "manifest entries")
        for name, payload in (("metadata.json", metadata_bytes), ("arrays.npz", arrays_bytes)):
            entry = expected[name]
            _exact_keys(entry, {"size", "sha256"}, f"manifest entry {name}")
            digest = entry["sha256"]
            if (
                type(entry["size"]) is not int
                or entry["size"] < 0
                or not isinstance(digest, str)
                or len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)
            ):
                raise SerializationError(f"manifest entry is invalid for {name}")
            if entry["size"] != len(payload) or digest != sha256(payload).hexdigest():
                raise SerializationError(f"model integrity check failed for {name}")
        metadata = _load_json(metadata_bytes, "metadata.json")
        if (
            metadata.get("format") != MODEL_FORMAT
            or type(metadata.get("format_version")) is not int
            or metadata.get("format_version") != MODEL_FORMAT_VERSION
        ):
            raise SerializationError("metadata declares an unsupported model format")
        arrays = _load_array_archive(arrays_bytes)
    except SerializationError:
        raise
    except Exception as exc:
        raise SerializationError(f"invalid model metadata or arrays: {exc}") from exc
    return metadata, arrays


def load_model(path: str | Path) -> WitnessCell:
    """Load and integrity-check a pickle-free ``.wcell`` model bundle."""
    source = Path(path)
    metadata, arrays = _read_bundle(source)
    try:
        _exact_keys(
            metadata,
            {
                "format",
                "format_version",
                "package_version",
                "algorithm_contract",
                "config",
                "state",
            },
            "metadata",
        )
        if not isinstance(metadata["package_version"], str) or not metadata[
            "package_version"
        ]:
            raise SerializationError("package_version must be a non-empty string")
        if not isinstance(metadata["algorithm_contract"], str):
            raise SerializationError("algorithm_contract must be a string")
        config_payload = metadata["config"]
        _exact_keys(
            config_payload,
            {"control_label", "alpha_grid", "noise_grid", "gate_confidence"},
            "model config",
        )
        if not isinstance(config_payload["control_label"], str):
            raise SerializationError("config control_label must be a string")
        for grid_name in ("alpha_grid", "noise_grid"):
            grid = config_payload[grid_name]
            if not isinstance(grid, list) or any(
                isinstance(item, bool) or not isinstance(item, (int, float)) for item in grid
            ):
                raise SerializationError(f"config {grid_name} must be a numeric array")
        config = WitnessCellConfig(
            control_label=config_payload["control_label"],
            alpha_grid=tuple(map(float, config_payload["alpha_grid"])),
            noise_grid=tuple(map(float, config_payload["noise_grid"])),
            gate_confidence=_number(
                config_payload["gate_confidence"], "config gate_confidence"
            ),
        )
        payload = metadata["state"]
        _exact_keys(
            payload,
            {
                "genes",
                "control_label",
                "known_single_nodes",
                "gene2go",
                "endpoint_head",
                "alpha",
                "noise_ratio",
                "gamma",
                "train_doubles",
                "validation_doubles",
                "validation_mse_factorized",
                "validation_mse_witness",
            },
            "model state",
        )
        for field in ("genes", "known_single_nodes", "train_doubles", "validation_doubles"):
            if not isinstance(payload[field], list) or not all(
                isinstance(item, str) for item in payload[field]
            ):
                raise SerializationError(f"state {field} must be an array of strings")
        if not isinstance(payload["control_label"], str):
            raise SerializationError("state control_label must be a string")
        if not isinstance(payload["gene2go"], dict) or not all(
            isinstance(gene, str)
            and isinstance(terms, list)
            and all(isinstance(term, str) for term in terms)
            for gene, terms in payload["gene2go"].items()
        ):
            raise SerializationError("state gene2go must map strings to string arrays")
        if not isinstance(payload["endpoint_head"], dict):
            raise SerializationError("state endpoint_head must be a JSON object")

        genes = tuple(payload["genes"])
        nodes = tuple(payload["known_single_nodes"])
        known_matrix = arrays["known_single_effects"]
        if arrays["control_mean"].shape != (len(genes),) or known_matrix.shape != (len(nodes), len(genes)):
            raise SerializationError("model arrays do not match declared genes and nodes")
        train_doubles = tuple(payload["train_doubles"])
        if arrays["train_residual"].shape != (len(train_doubles), len(genes)):
            raise SerializationError("training residual matrix has an invalid shape")
        known_effects = MappingProxyType(
            {node: known_matrix[index] for index, node in enumerate(nodes)}
        )
        gene2go = MappingProxyType(
            {
                gene: tuple(terms)
                for gene, terms in payload["gene2go"].items()
            }
        )
        state = FittedWitnessCore(
            genes=genes,
            control_label=payload["control_label"],
            control_mean=arrays["control_mean"],
            known_single_effects=known_effects,
            gene2go=gene2go,
            endpoint_head=_restore_head(payload["endpoint_head"]),
            alpha=_number(payload["alpha"], "state alpha"),
            noise_ratio=_number(payload["noise_ratio"], "state noise_ratio"),
            gamma=_number(payload["gamma"], "state gamma"),
            train_doubles=train_doubles,
            train_residual=arrays["train_residual"],
            validation_doubles=tuple(payload["validation_doubles"]),
            validation_mse_factorized=float(payload["validation_mse_factorized"]),
            validation_mse_witness=float(payload["validation_mse_witness"]),
        )
        _validate_model_state(
            config,
            state,
            metadata["algorithm_contract"],
        )
        model = WitnessCell(config)
        model._state = state
        return model
    except SerializationError:
        raise
    except Exception as exc:
        raise SerializationError(f"model state could not be reconstructed: {exc}") from exc
