#!/usr/bin/env python3
"""Score one paired GEARS/CPA run under the frozen six-way contract."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


EPS = 1e-12
STRATEGIES = (
    "Saturation", "GEARS", "CPA", "GEARS+Witness", "CPA+Witness", "WitnessCell"
)


def scalar_text(array: np.ndarray) -> str:
    return str(np.asarray(array).item())


def scalar_int(array: np.ndarray) -> int:
    return int(np.asarray(array).item())


def load_archive(path: Path, expected_model: str) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        output = {key: archive[key] for key in archive.files}
    required = {
        "model", "model_version", "data_seed", "model_seed",
        "conditions", "genes", "pred_effect", "protocol_sha256",
        "training_budget_sha256",
    }
    if set(output) != required:
        raise ValueError(f"{path}: archive keys differ from frozen schema")
    if scalar_text(output["model"]) != expected_model:
        raise ValueError(f"{path}: expected model {expected_model}")
    if output["pred_effect"].shape != (len(output["conditions"]), len(output["genes"])):
        raise ValueError(f"{path}: pred_effect shape mismatch")
    if output["genes"].dtype.kind == "O" or output["conditions"].dtype.kind == "O":
        raise ValueError(f"{path}: object strings are forbidden")
    if not np.isfinite(output["pred_effect"]).all():
        raise ValueError(f"{path}: non-finite prediction")
    return output


def align_prediction(archive: dict[str, np.ndarray], conditions: np.ndarray, genes: np.ndarray) -> np.ndarray:
    condition_id = {value: index for index, value in enumerate(archive["conditions"].astype(str))}
    gene_id = {value: index for index, value in enumerate(archive["genes"].astype(str))}
    missing_conditions = [value for value in conditions if value not in condition_id]
    missing_genes = [value for value in genes if value not in gene_id]
    if missing_conditions or missing_genes:
        raise ValueError(
            f"archive alignment failure: missing conditions={missing_conditions}, "
            f"missing genes={missing_genes[:10]}"
        )
    rows = np.asarray([condition_id[value] for value in conditions], dtype=int)
    columns = np.asarray([gene_id[value] for value in genes], dtype=int)
    return archive["pred_effect"][rows][:, columns].astype(float)


def closed_form_fusion_lambda(base: np.ndarray, witness: np.ndarray, truth: np.ndarray) -> float:
    direction = witness - base
    denominator = float(np.sum(direction * direction))
    if denominator <= EPS:
        return 0.0
    return float(np.clip(np.sum(direction * (truth - base)) / denominator, 0.0, 1.0))


def row_cosine(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    numerator = np.sum(left * right, axis=1)
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    return np.divide(
        numerator, denominator, out=np.zeros_like(numerator, dtype=float), where=denominator > EPS
    )


def summarize(residual: np.ndarray, truth: np.ndarray, baseline: np.ndarray) -> dict:
    full = baseline + residual
    full_truth = baseline + truth
    return {
        "residual_mse": float(np.mean((residual - truth) ** 2)),
        "residual_cosine": float(np.sum(residual * truth) / max(
            np.linalg.norm(residual) * np.linalg.norm(truth), EPS
        )),
        "full_effect_mse": float(np.mean((full - full_truth) ** 2)),
        "full_effect_cosine": float(np.sum(full * full_truth) / max(
            np.linalg.norm(full) * np.linalg.norm(full_truth), EPS
        )),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gears", type=Path, required=True)
    parser.add_argument("--cpa", type=Path, required=True)
    parser.add_argument("--anchor", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, default=Path(__file__).resolve().parent / "FROZEN_SCIENTIFIC_PROTOCOL.json")
    parser.add_argument("--budget", type=Path, default=Path(__file__).resolve().parent / "config/FROZEN_TRAINING_BUDGET.json")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    import hashlib
    protocol_hash = hashlib.sha256(args.protocol.read_bytes()).hexdigest()
    budget_hash = hashlib.sha256(args.budget.read_bytes()).hexdigest()
    gears = load_archive(args.gears, "GEARS")
    cpa = load_archive(args.cpa, "CPA")
    if scalar_int(gears["data_seed"]) != scalar_int(cpa["data_seed"]):
        raise ValueError("GEARS and CPA data seeds differ")
    if scalar_int(gears["model_seed"]) != scalar_int(cpa["model_seed"]):
        raise ValueError("GEARS and CPA model seeds differ")
    if scalar_text(gears["protocol_sha256"]) != protocol_hash:
        raise ValueError("GEARS protocol hash mismatch")
    if scalar_text(cpa["protocol_sha256"]) != protocol_hash:
        raise ValueError("CPA protocol hash mismatch")
    if scalar_text(gears["training_budget_sha256"]) != budget_hash:
        raise ValueError("GEARS training-budget hash mismatch")
    if scalar_text(cpa["training_budget_sha256"]) != budget_hash:
        raise ValueError("CPA training-budget hash mismatch")

    with np.load(args.anchor, allow_pickle=False) as stored:
        anchor = {key: stored[key] for key in stored.files}
    data_seed = scalar_int(gears["data_seed"])
    model_seed = scalar_int(gears["model_seed"])
    if scalar_int(anchor["data_seed"]) != data_seed:
        raise ValueError("anchor data seed differs from model archives")
    genes = anchor["genes"].astype(str)
    calibration_pairs = anchor["calibration_pairs"].astype(str)
    final_pairs = anchor["final_pairs"].astype(str)
    baseline_cal = anchor["baseline_calibration"].astype(float)
    truth_cal = anchor["truth_residual_calibration"].astype(float)
    witness_cal = float(anchor["gamma_witness"]) * anchor["witness_raw_calibration"].astype(float)
    baseline_final = anchor["baseline_final"].astype(float)
    truth_final = anchor["truth_residual_final"].astype(float)
    witness_final = anchor["witness_calibrated_final"].astype(float)

    gears_cal_effect = align_prediction(gears, calibration_pairs, genes)
    gears_final_effect = align_prediction(gears, final_pairs, genes)
    cpa_cal_effect = align_prediction(cpa, calibration_pairs, genes)
    cpa_final_effect = align_prediction(cpa, final_pairs, genes)
    gears_cal = gears_cal_effect - baseline_cal
    gears_final = gears_final_effect - baseline_final
    cpa_cal = cpa_cal_effect - baseline_cal
    cpa_final = cpa_final_effect - baseline_final

    lambda_gears = closed_form_fusion_lambda(gears_cal, witness_cal, truth_cal)
    lambda_cpa = closed_form_fusion_lambda(cpa_cal, witness_cal, truth_cal)
    gears_witness = (1 - lambda_gears) * gears_final + lambda_gears * witness_final
    cpa_witness = (1 - lambda_cpa) * cpa_final + lambda_cpa * witness_final
    residuals = {
        "Saturation": np.zeros_like(truth_final),
        "GEARS": gears_final,
        "CPA": cpa_final,
        "GEARS+Witness": gears_witness,
        "CPA+Witness": cpa_witness,
        "WitnessCell": witness_final,
    }

    state_metrics = {strategy: summarize(value, truth_final, baseline_final) for strategy, value in residuals.items()}
    target_rows: list[dict] = []
    for strategy, prediction in residuals.items():
        residual_mse = np.mean((prediction - truth_final) ** 2, axis=1)
        full_prediction = baseline_final + prediction
        full_truth = baseline_final + truth_final
        for index, pair in enumerate(final_pairs):
            target_rows.append({
                "data_seed": data_seed,
                "model_seed": model_seed,
                "pair": pair,
                "strategy": strategy,
                "residual_mse": float(residual_mse[index]),
                "residual_cosine": float(row_cosine(prediction[index:index+1], truth_final[index:index+1])[0]),
                "full_effect_mse": float(np.mean((full_prediction[index] - full_truth[index]) ** 2)),
                "full_effect_cosine": float(row_cosine(
                    full_prediction[index:index+1], full_truth[index:index+1]
                )[0]),
            })

    calibration_losses = {
        "GEARS": float(np.mean((gears_cal - truth_cal) ** 2)),
        "GEARS+Witness": float(np.mean(((1-lambda_gears)*gears_cal + lambda_gears*witness_cal - truth_cal) ** 2)),
        "CPA": float(np.mean((cpa_cal - truth_cal) ** 2)),
        "CPA+Witness": float(np.mean(((1-lambda_cpa)*cpa_cal + lambda_cpa*witness_cal - truth_cal) ** 2)),
        "WitnessCell": float(np.mean((witness_cal - truth_cal) ** 2)),
    }
    if calibration_losses["GEARS+Witness"] > calibration_losses["GEARS"] + 1e-12:
        raise RuntimeError("closed-form GEARS fusion is not optimal on calibration data")
    if calibration_losses["CPA+Witness"] > calibration_losses["CPA"] + 1e-12:
        raise RuntimeError("closed-form CPA fusion is not optimal on calibration data")

    verdict = {
        "status": "SIX_WAY_RUN_SCORED",
        "data_seed": data_seed,
        "model_seed": model_seed,
        "protocol_sha256": protocol_hash,
        "training_budget_sha256": budget_hash,
        "models": {
            "GEARS": scalar_text(gears["model_version"]),
            "CPA": scalar_text(cpa["model_version"]),
        },
        "lambda_gears": lambda_gears,
        "lambda_cpa": lambda_cpa,
        "gamma_witness": float(anchor["gamma_witness"]),
        "calibration_losses": calibration_losses,
        "strategies": state_metrics,
        "claim_boundary": "One paired model-seed/data-seed run. Formal inference requires all 3x3 runs and never treats cells or model seeds as biological replicates."
    }
    pd.DataFrame(target_rows).to_csv(args.out / "per_target.csv", index=False)
    np.savez_compressed(
        args.out / "six_way_predictions.npz",
        pairs=final_pairs,
        genes=genes,
        truth_residual=truth_final.astype(np.float32),
        saturation_baseline=baseline_final.astype(np.float32),
        **{strategy.replace("+", "_plus_").lower(): value.astype(np.float32) for strategy, value in residuals.items()},
    )
    (args.out / "metrics.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
