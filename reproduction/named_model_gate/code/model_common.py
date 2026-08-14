"""Shared contracts for the GSE146194 named-model runners."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


ARCHIVE_KEYS = {
    "model", "model_version", "data_seed", "model_seed", "conditions",
    "genes", "pred_effect", "protocol_sha256", "training_budget_sha256",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_condition(condition: str) -> str:
    tokens = [
        token for token in str(condition).split("+")
        if token and token.lower() not in {"ctrl", "control", "ntc", "negctrl"}
    ]
    if not tokens:
        return "CTRL"
    return "+".join(sorted(set(tokens)))


def load_contract(protocol_path: Path, budget_path: Path) -> tuple[dict, dict, str, str]:
    protocol = json.loads(protocol_path.read_text())
    budget = json.loads(budget_path.read_text())
    if not protocol.get("frozen_before_named_model_predictions"):
        raise ValueError("scientific protocol is not frozen")
    if not budget.get("frozen_before_any_named_target_prediction"):
        raise ValueError("training budget is not frozen")
    return protocol, budget, sha256(protocol_path), sha256(budget_path)


def requested_conditions(protocol: dict) -> np.ndarray:
    return np.asarray(
        protocol["condition_split"]["witness_calibration_pairs"]
        + protocol["condition_split"]["final_test_pairs"], dtype=str,
    )


def save_prediction_archive(
    path: Path,
    *,
    model: str,
    model_version: str,
    data_seed: int,
    model_seed: int,
    conditions: np.ndarray,
    genes: np.ndarray,
    pred_effect: np.ndarray,
    protocol_sha256: str,
    training_budget_sha256: str,
) -> None:
    conditions = np.asarray(conditions, dtype=str)
    genes = np.asarray(genes, dtype=str)
    pred_effect = np.asarray(pred_effect, dtype=np.float32)
    if pred_effect.shape != (len(conditions), len(genes)):
        raise ValueError("prediction shape differs from conditions x genes")
    if len(set(conditions.tolist())) != len(conditions):
        raise ValueError("duplicate output conditions")
    if len(set(genes.tolist())) != len(genes):
        raise ValueError("duplicate output genes")
    if not np.all(np.isfinite(pred_effect)):
        raise ValueError("non-finite model prediction")
    payload = {
        "model": np.asarray(model),
        "model_version": np.asarray(model_version),
        "data_seed": np.asarray(data_seed),
        "model_seed": np.asarray(model_seed),
        "conditions": conditions,
        "genes": genes,
        "pred_effect": pred_effect,
        "protocol_sha256": np.asarray(protocol_sha256),
        "training_budget_sha256": np.asarray(training_budget_sha256),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **payload)
    with np.load(path, allow_pickle=False) as stored:
        if set(stored.files) != ARCHIVE_KEYS:
            raise RuntimeError("written archive schema changed")
        for key in stored.files:
            if stored[key].dtype.kind == "O":
                raise RuntimeError(f"object dtype in prediction archive: {key}")
