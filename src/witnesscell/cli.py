"""Command-line interface for fit, predict, inspect, and selective decisions."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from ._version import __version__
from .data import ConditionMoments, SplitSpec
from .exceptions import WitnessCellError
from .model import WitnessCell
from .selective import SelectivePolicy


def _read_json(path: str | Path) -> Any:
    try:
        with Path(path).open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read JSON from {path}: {exc}") from exc


def _write_json(value: Any) -> None:
    def compatible(item: Any) -> Any:
        if isinstance(item, (bool, np.bool_)):
            return bool(item)
        if isinstance(item, (float, np.floating)):
            number = float(item)
            return number if np.isfinite(number) else None
        if isinstance(item, (int, np.integer)):
            return int(item)
        if isinstance(item, dict):
            return {str(key): compatible(child) for key, child in item.items()}
        if isinstance(item, (tuple, list)):
            return [compatible(child) for child in item]
        return item

    print(json.dumps(compatible(value), ensure_ascii=False, allow_nan=False, indent=2, sort_keys=True))


def _conditions(arguments: argparse.Namespace) -> tuple[str, ...]:
    if arguments.conditions:
        return tuple(arguments.conditions)
    payload = _read_json(arguments.conditions_file)
    if not isinstance(payload, list) or not all(isinstance(item, str) for item in payload):
        raise ValueError("conditions file must be a JSON array of strings")
    return tuple(payload)


def _fit(arguments: argparse.Namespace) -> int:
    moments = ConditionMoments.from_npz(arguments.moments)
    split_payload = _read_json(arguments.split)
    if not isinstance(split_payload, dict):
        raise ValueError("split file must contain a JSON object")
    split = SplitSpec.create(
        split_payload.get("train_conditions", ()),
        split_payload.get("validation_conditions", ()),
    )
    gene2go = _read_json(arguments.gene2go)
    if not isinstance(gene2go, dict):
        raise ValueError("gene2go file must contain a JSON object")
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    output = model.save(arguments.output)
    result = model.diagnostics()
    result["model"] = str(output)
    _write_json(result)
    return 0


def _predict(arguments: argparse.Namespace) -> int:
    model = WitnessCell.load(arguments.model)
    conditions = _conditions(arguments)
    output = Path(arguments.output)
    model.predict(conditions).to_npz(output)
    _write_json({"conditions": len(conditions), "genes": len(model.state.genes), "output": str(output)})
    return 0


def _inspect(arguments: argparse.Namespace) -> int:
    _write_json(WitnessCell.load(arguments.model).diagnostics())
    return 0


def _validate(arguments: argparse.Namespace) -> int:
    moments = ConditionMoments.from_npz(arguments.moments)
    _write_json({"valid": True, "conditions": len(moments.means), "genes": len(moments.genes)})
    return 0


def _decide(arguments: argparse.Namespace) -> int:
    risks = np.asarray(np.load(arguments.risks, allow_pickle=False), dtype=float)
    if risks.ndim != 1:
        raise ValueError("risk file must be a one-dimensional .npy array")
    if arguments.frozen_norman:
        policy = SelectivePolicy.frozen_norman()
    else:
        policy = SelectivePolicy(float(arguments.threshold), "cli-user-supplied")
    decisions = policy.decide(risks)
    np.savez_compressed(
        arguments.output,
        risks=risks,
        decisions=np.asarray(decisions, dtype=str),
        threshold=np.asarray(policy.threshold, dtype=float),
        provenance=np.asarray(policy.provenance, dtype=str),
    )
    _write_json({"accepted": decisions.count("accept"), "abstained": decisions.count("abstain"), "output": str(arguments.output), "threshold": policy.threshold})
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="witnesscell", description="Evidence-gated condition-mean prediction")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    fit = subparsers.add_parser("fit", help="fit and save a .wcell model")
    fit.add_argument("--moments", required=True, help="condition-moment .npz file")
    fit.add_argument("--split", required=True, help="training/validation split JSON")
    fit.add_argument("--gene2go", required=True, help="gene-to-GO JSON mapping")
    fit.add_argument("--output", required=True, help="output path ending in .wcell")
    fit.set_defaults(handler=_fit)

    predict = subparsers.add_parser("predict", help="predict condition means")
    predict.add_argument("--model", required=True)
    source = predict.add_mutually_exclusive_group(required=True)
    source.add_argument("--conditions", nargs="+")
    source.add_argument("--conditions-file")
    predict.add_argument("--output", required=True, help="pickle-free output .npz")
    predict.set_defaults(handler=_predict)

    inspect = subparsers.add_parser("inspect", help="show fitted diagnostics")
    inspect.add_argument("--model", required=True)
    inspect.set_defaults(handler=_inspect)

    validate = subparsers.add_parser("validate", help="validate a moment archive")
    validate.add_argument("--moments", required=True)
    validate.set_defaults(handler=_validate)

    decide = subparsers.add_parser("decide", help="apply an accept/abstain threshold")
    decide.add_argument("--risks", required=True, help="one-dimensional .npy risk vector")
    threshold = decide.add_mutually_exclusive_group(required=True)
    threshold.add_argument("--threshold", type=float)
    threshold.add_argument("--frozen-norman", action="store_true", help="retrospective Gate 21 threshold; no transport guarantee")
    decide.add_argument("--output", required=True)
    decide.set_defaults(handler=_decide)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    try:
        arguments = parser.parse_args(argv)
        return int(arguments.handler(arguments))
    except (WitnessCellError, ValueError, OSError) as exc:
        parser.exit(2, f"witnesscell: error: {exc}\n")


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
