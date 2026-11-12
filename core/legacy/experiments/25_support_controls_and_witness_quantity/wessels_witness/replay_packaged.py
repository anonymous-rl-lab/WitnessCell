#!/usr/bin/env python3
"""Verify packaged predictions/losses from frozen numerical state; NumPy only.

No GO assets, endpoint training code, old repository or network are needed.
Inputs are read-only. An optional new JSON report can be written with --out.
Exit 0 means PASS; exit 1 means verification failed.
"""
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np

PREDICTION_ATOL = 1e-6
PREDICTION_RTOL = 1e-6
LOSS_ATOL = 1e-12
LOSS_RTOL = 1e-10


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def local(root, relative):
    path = (root / relative).resolve()
    path.relative_to(root)
    return path


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def npz(path):
    with np.load(path, allow_pickle=False) as z:
        return {key: z[key] for key in z.files}


def pair(condition):
    return "+".join(sorted(str(condition).split("+")))


def verify(root):
    t = time.monotonic()
    spec = read_json(root / "protocol.json")
    protocol_lock = read_json(root / "protocol_freeze.json")
    require(sha(root / "protocol.json") == protocol_lock["protocol_sha256"], "Protocol hash mismatch")
    require(sha(root / "subsets.json") == spec["subsets_sha256"], "Subset hash mismatch")
    freeze = read_json(root / "predictions_frozen.json")
    require(freeze["config_hash"] == protocol_lock["protocol_sha256"], "Prediction configuration mismatch")
    for relative, expected in freeze["inventory"].items():
        require(sha(local(root, relative)) == expected, f"Frozen inventory hash mismatch: {relative}")
    eval_manifest = read_json(root / "evaluation/evaluation_manifest.json")
    require(sha(root / "predictions_frozen.json") == eval_manifest["prediction_freeze_sha256"], "Evaluation prediction-freeze hash mismatch")
    for relative, expected in eval_manifest["output_hashes"].items():
        require(sha(local(root, "evaluation/" + relative)) == expected, f"Saved evaluation hash mismatch: {relative}")
    copied = read_json(root / "evaluation_inputs/source_manifest.json")
    require(copied["protocol_sha256"] == protocol_lock["protocol_sha256"], "Copied evaluation configuration mismatch")
    eval_inputs = {}
    for entry in copied["files"]:
        seed = int(entry["seed"])
        expected = next(v["original_evaluation_sha256"] for v in spec["inputs"] if int(v["seed"]) == seed)
        path = local(root, entry["relative_path"])
        require(entry["sha256"] == expected and sha(path) == expected, f"Evaluation input hash mismatch: seed{seed}")
        require(path.stat().st_size == entry["size_bytes"], f"Evaluation byte-size mismatch: seed{seed}")
        eval_inputs[seed] = path
    require(set(eval_inputs) == {1, 2, 3}, "Three evaluation inputs are required")
    subsets = {(int(r["seed"]), int(r["order"]), float(r["fraction"])): r for r in read_json(root / "subsets.json")}
    require(len(subsets) == 27, "Expected 27 unique logical subsets")
    logical = freeze["logical_predictions"]
    require(len(logical) == 27, "Expected 27 logical predictions")
    states = {seed: npz(root / f"seed{seed}/frozen_state.npz") for seed in (1, 2, 3)}
    predictions, prediction_checks = {}, []
    for setting in logical:
        seed, order, fraction = int(setting["seed"]), int(setting["order"]), float(setting["fraction"])
        subset = subsets[(seed, order, fraction)]
        state = states[seed]
        training = state["training_combinations"].astype(str).tolist()
        chosen = subset["chosen_conditions"]
        selected = set(chosen)
        indices = np.array([i for i, condition in enumerate(training) if condition in selected], dtype=int)
        require([training[i] for i in indices] == chosen, "Chosen witness order does not match fixed original row order")
        require(len(indices) == subset["size"], "Subset-size mismatch")
        relative = setting["prediction"]
        if relative in predictions:
            require(predictions[relative]["selected"] == chosen, "A reused prediction maps to different witnesses")
            continue
        design = state["train_design"][indices]
        train_kernel = design @ design.T / 2.0
        cross = design @ state["target_design"].T / 2.0
        coefficient = np.linalg.solve(train_kernel + float(state["noise"]) * np.eye(len(indices)), cross).T
        residual_correction = float(state["gamma"]) * (coefficient @ state["residual"][indices])
        replay = (state["backbone"] + residual_correction).astype(np.float32)
        saved = npz(local(root, relative))
        require(np.array_equal(saved["conditions"], state["conditions"]), "Prediction query axis differs from frozen state")
        require(np.array_equal(saved["genes"], state["genes"]), "Prediction gene axis differs from frozen state")
        require(saved["selected_conditions"].astype(str).tolist() == chosen, "Prediction's selected witness identities differ")
        require(replay.shape == saved["prediction"].shape, "Prediction shape mismatch")
        equal = np.isclose(replay, saved["prediction"], atol=PREDICTION_ATOL, rtol=PREDICTION_RTOL)
        require(bool(equal.all()), f"Prediction elementwise mismatch: {relative}")
        witnesses = (cross > 0).sum(axis=0).astype(int)
        require(np.array_equal(witnesses, saved["witness_count"]), f"Incident-witness count mismatch: {relative}")
        require(np.all(residual_correction[witnesses == 0] == 0), "Zero-incident query has a nonzero residual correction")
        predictions[relative] = {"prediction": replay, "selected": chosen, "witness_count": witnesses}
        prediction_checks.append(dict(path=relative, shape=list(replay.shape),
            float32_exact=bool(np.array_equal(replay, saved["prediction"])),
            exact_unequal_elements=int(np.count_nonzero(replay != saved["prediction"])),
            max_abs_error=float(np.max(np.abs(replay.astype(float) - saved["prediction"].astype(float)))),
            all_elements_within_tolerance=True))
    require(len(predictions) == 21, "Expected 21 unique stored predictions")
    rows = read_csv(root / "evaluation/query_losses.csv")
    comparator_rows = read_csv(root / "evaluation/fixed_comparator_losses.csv")
    row_map = {(int(r["seed"]), int(r["order"]), float(r["fraction"]), r["pair"]): r for r in rows}
    comparator_map = {(int(r["seed"]), r["method"], r["pair"]): r for r in comparator_rows}
    require(len(rows) == len(row_map) == 1629, "Expected 1629 unique query loss rows")
    require(len(comparator_rows) == len(comparator_map) == 543, "Expected 543 unique comparator loss rows")
    loss_max_difference = {"loss_mse": 0.0, "loss_wmse": 0.0}
    seen_rows, seen_comparators = set(), set()
    for seed in (1, 2, 3):
        state = states[seed]
        evaluation = npz(eval_inputs[seed])
        require(np.array_equal(evaluation["conditions"], state["conditions"]), "Evaluation query axis mismatch")
        require(np.array_equal(evaluation["genes"], state["genes"]), "Evaluation gene axis mismatch")
        truth, weights = evaluation["truth"].astype(float), evaluation["weights"].astype(float)
        conditions = state["conditions"].astype(str).tolist()
        gi = {str(g): i for i, g in enumerate(state["genes"])}
        mask = np.ones(truth.shape, dtype=bool)
        for i, condition in enumerate(conditions):
            for gene in condition.split("+"):
                mask[i, gi[gene]] = False
        require(np.all(mask.sum(axis=1) == truth.shape[1] - 2), "Both target coordinates must be excluded")
        retained_weight = weights * mask
        denom = retained_weight.sum(axis=1)
        require(np.all(denom > 0), "Nonpositive retained evaluation weight")

        def losses(prediction):
            error = (prediction.astype(float) - truth) ** 2
            return {"loss_mse": (error * mask).sum(axis=1) / mask.sum(axis=1),
                    "loss_wmse": (error * retained_weight).sum(axis=1) / denom}

        def compare_losses(actual, stored, index, label):
            for metric in loss_max_difference:
                value, expected = float(actual[metric][index]), float(stored[metric])
                loss_max_difference[metric] = max(loss_max_difference[metric], abs(value - expected))
                require(np.isclose(value, expected, atol=LOSS_ATOL, rtol=LOSS_RTOL), f"Loss mismatch: {label}/{metric}")

        for setting in logical:
            if int(setting["seed"]) != seed:
                continue
            replay = predictions[setting["prediction"]]
            actual = losses(replay["prediction"])
            for i, condition in enumerate(conditions):
                key = (seed, int(setting["order"]), float(setting["fraction"]), pair(condition))
                stored = row_map[key]
                require(stored["prediction_path"] == setting["prediction"] and int(stored["prediction_row"]) == i,
                        "CSV prediction-file/row mapping mismatch")
                require(int(stored["witness_count"]) == int(replay["witness_count"][i]), "CSV incident-witness count mismatch")
                compare_losses(actual, stored, i, str(key))
                seen_rows.add(key)
        for method, prediction in (("factorized_prediction", state["original_factorized"]),
                                   ("shared_endpoint_additive", state["shared_additive"]),
                                   ("original_final_v14", state["original_final"])):
            actual = losses(prediction)
            for i, condition in enumerate(conditions):
                key = (seed, method, pair(condition))
                compare_losses(actual, comparator_map[key], i, str(key))
                seen_comparators.add(key)
    require(len(seen_rows) == 1629 and len(seen_comparators) == 543, "Incomplete loss replay")
    return dict(status="PASS", scope="Prediction and loss reproduction from frozen numerical state only",
        endpoint_training_reproduced=False,
        endpoint_training_requirement="Reproducing the original endpoint training still requires the old training/GO/source assets; this verifier does not validate that training stage",
        old_repo_or_go_access_required=False, input_files_modified=False,
        unique_predictions_replayed=len(predictions), logical_settings=len(logical),
        prediction_elements_replayed=sum(np.prod(r["shape"]).item() for r in prediction_checks),
        all_predictions_float32_exact=all(r["float32_exact"] for r in prediction_checks),
        prediction_checks=prediction_checks, query_loss_rows_verified=len(seen_rows),
        fixed_comparator_loss_rows_verified=len(seen_comparators),
        loss_max_abs_difference=loss_max_difference,
        tolerances=dict(prediction_atol=PREDICTION_ATOL, prediction_rtol=PREDICTION_RTOL,
                        loss_atol=LOSS_ATOL, loss_rtol=LOSS_RTOL),
        evaluation_input_bytes=sum(path.stat().st_size for path in eval_inputs.values()),
        elapsed_seconds=time.monotonic()-t)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--out", type=Path, help="Optional new JSON report; existing files are never overwritten")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.out:
        args.out = args.out.resolve()
        args.out.relative_to(root)
        if args.out.exists():
            parser.error("Refusing to overwrite an existing report or package input")
    started = time.time()
    try:
        report = verify(root)
        code = 0
    except Exception as error:
        report = dict(status="FAIL", error_type=type(error).__name__, error=str(error), input_files_modified=False)
        code = 1
    report.update(started_unix=started, finished_unix=time.time(), source_sha256=sha(Path(__file__).resolve()),
                  numpy_version=np.__version__, python_version=sys.version)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "prediction_checks"}, ensure_ascii=False, allow_nan=False))
    raise SystemExit(code)


if __name__ == "__main__":
    main()
