#!/usr/bin/env python3
"""Independent, target-level audit of the frozen GSE146194 named-AIVC gate.

This module deliberately does not import the experiment scorer or summarizer.
It reloads the raw GEARS/CPA archives and CPU anchors, reconstructs all six
strategies, verifies every packaged target metric, and performs inference at
the level of the 14 biological target pairs.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


EPS = 1e-12
DATA_SEEDS = (41, 47, 59)
MODEL_SEEDS = (0, 1, 2)
STRATEGIES = (
    "Saturation",
    "GEARS",
    "CPA",
    "GEARS+Witness",
    "CPA+Witness",
    "WitnessCell",
)
FROZEN_COMPARISONS = (
    ("GEARS+Witness", "GEARS"),
    ("CPA+Witness", "CPA"),
    ("WitnessCell", "GEARS"),
    ("WitnessCell", "CPA"),
)
CONSERVATIVE_COMPARISONS = (
    ("WitnessCell", "Saturation"),
    ("GEARS+Witness", "Saturation"),
    ("CPA+Witness", "Saturation"),
    ("GEARS+Witness", "WitnessCell"),
    ("CPA+Witness", "WitnessCell"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scalar(array: np.ndarray):
    return np.asarray(array).item()


def load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        result = {key: archive[key] for key in archive.files}
    for key, value in result.items():
        if value.dtype.kind == "O":
            raise AssertionError(f"pickle/object array is forbidden: {path}:{key}")
    return result


def align(archive: dict[str, np.ndarray], conditions: np.ndarray, genes: np.ndarray) -> np.ndarray:
    condition_index = {str(value): index for index, value in enumerate(archive["conditions"])}
    gene_index = {str(value): index for index, value in enumerate(archive["genes"])}
    missing_conditions = [str(value) for value in conditions if str(value) not in condition_index]
    missing_genes = [str(value) for value in genes if str(value) not in gene_index]
    if missing_conditions or missing_genes:
        raise AssertionError(
            f"alignment failure: conditions={missing_conditions}, genes={missing_genes[:10]}"
        )
    rows = np.asarray([condition_index[str(value)] for value in conditions], dtype=int)
    cols = np.asarray([gene_index[str(value)] for value in genes], dtype=int)
    return np.asarray(archive["pred_effect"], dtype=float)[rows][:, cols]


def fusion_weight(base: np.ndarray, witness: np.ndarray, truth: np.ndarray) -> float:
    direction = witness - base
    denominator = float(np.sum(direction * direction))
    if denominator <= EPS:
        return 0.0
    return float(np.clip(np.sum(direction * (truth - base)) / denominator, 0.0, 1.0))


def cosine_rows(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    numerator = np.sum(left * right, axis=1)
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > EPS)


def exact_signflip_p(differences: np.ndarray) -> float:
    observed = float(np.mean(differences))
    signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=len(differences))))
    randomized = np.mean(signs * differences[None, :], axis=1)
    return float(np.mean(randomized >= observed - 1e-15))


def reconstruct_run(
    raw_root: Path,
    data_seed: int,
    model_seed: int,
    protocol_hash: str,
    budget_hash: str,
) -> tuple[pd.DataFrame, dict]:
    gears_path = raw_root / f"predictions/formal/GEARS_d{data_seed}_m{model_seed}.npz"
    cpa_path = raw_root / f"predictions/formal/CPA_d{data_seed}_m{model_seed}.npz"
    anchor_path = raw_root / f"cpu_anchor/truth_anchor_d{data_seed}.npz"
    gears = load_npz(gears_path)
    cpa = load_npz(cpa_path)
    anchor = load_npz(anchor_path)

    expected_keys = {
        "model", "model_version", "data_seed", "model_seed", "conditions", "genes",
        "pred_effect", "protocol_sha256", "training_budget_sha256",
    }
    for name, archive in (("GEARS", gears), ("CPA", cpa)):
        if set(archive) != expected_keys:
            raise AssertionError(f"{name} archive schema mismatch at d{data_seed}/m{model_seed}")
        assert str(scalar(archive["model"])) == name
        assert int(scalar(archive["data_seed"])) == data_seed
        assert int(scalar(archive["model_seed"])) == model_seed
        assert str(scalar(archive["protocol_sha256"])) == protocol_hash
        assert str(scalar(archive["training_budget_sha256"])) == budget_hash
        assert np.isfinite(archive["pred_effect"]).all()

    genes = anchor["genes"].astype(str)
    cal_pairs = anchor["calibration_pairs"].astype(str)
    final_pairs = anchor["final_pairs"].astype(str)
    baseline_cal = anchor["baseline_calibration"].astype(float)
    truth_cal = anchor["truth_residual_calibration"].astype(float)
    witness_cal = float(scalar(anchor["gamma_witness"])) * anchor["witness_raw_calibration"].astype(float)
    baseline_final = anchor["baseline_final"].astype(float)
    truth_final = anchor["truth_residual_final"].astype(float)
    witness_final = anchor["witness_calibrated_final"].astype(float)

    gears_cal = align(gears, cal_pairs, genes) - baseline_cal
    cpa_cal = align(cpa, cal_pairs, genes) - baseline_cal
    gears_final = align(gears, final_pairs, genes) - baseline_final
    cpa_final = align(cpa, final_pairs, genes) - baseline_final
    lambda_gears = fusion_weight(gears_cal, witness_cal, truth_cal)
    lambda_cpa = fusion_weight(cpa_cal, witness_cal, truth_cal)
    residuals = {
        "Saturation": np.zeros_like(truth_final),
        "GEARS": gears_final,
        "CPA": cpa_final,
        "GEARS+Witness": (1.0 - lambda_gears) * gears_final + lambda_gears * witness_final,
        "CPA+Witness": (1.0 - lambda_cpa) * cpa_final + lambda_cpa * witness_final,
        "WitnessCell": witness_final,
    }

    rows: list[dict] = []
    full_truth = baseline_final + truth_final
    for strategy, prediction in residuals.items():
        full_prediction = baseline_final + prediction
        for index, pair in enumerate(final_pairs):
            rows.append({
                "data_seed": data_seed,
                "model_seed": model_seed,
                "pair": str(pair),
                "strategy": strategy,
                "residual_mse": float(np.mean((prediction[index] - truth_final[index]) ** 2)),
                "residual_cosine": float(cosine_rows(prediction[index:index + 1], truth_final[index:index + 1])[0]),
                "full_effect_mse": float(np.mean((full_prediction[index] - full_truth[index]) ** 2)),
                "full_effect_cosine": float(cosine_rows(full_prediction[index:index + 1], full_truth[index:index + 1])[0]),
            })
    metadata = {
        "data_seed": data_seed,
        "model_seed": model_seed,
        "lambda_gears": lambda_gears,
        "lambda_cpa": lambda_cpa,
        "gamma_witness": float(scalar(anchor["gamma_witness"])),
        "gears_archive_sha256": sha256(gears_path),
        "cpa_archive_sha256": sha256(cpa_path),
        "anchor_sha256": sha256(anchor_path),
    }
    return pd.DataFrame(rows), metadata


def compare_rows(independent: pd.DataFrame, packaged: pd.DataFrame) -> float:
    keys = ["data_seed", "model_seed", "pair", "strategy"]
    metrics = ["residual_mse", "residual_cosine", "full_effect_mse", "full_effect_cosine"]
    joined = independent.merge(packaged, on=keys, how="outer", suffixes=("_audit", "_package"), indicator=True)
    if len(joined) != 756 or not (joined["_merge"] == "both").all():
        raise AssertionError("packaged and independently reconstructed target rows do not align")
    maximum = 0.0
    for metric in metrics:
        delta = np.abs(joined[f"{metric}_audit"] - joined[f"{metric}_package"])
        maximum = max(maximum, float(delta.max()))
    if maximum > 5e-10:
        raise AssertionError(f"target metric mismatch: max abs delta={maximum}")
    return maximum


def paired_comparison(
    targets: pd.DataFrame,
    candidate: str,
    comparator: str,
    bootstrap_seed: int,
    draws: int,
) -> dict:
    pair_mean = targets.groupby(["pair", "strategy"], as_index=False).agg(
        residual_mse=("residual_mse", "mean"),
        full_effect_cosine=("full_effect_cosine", "mean"),
    )
    mse = pair_mean.pivot(index="pair", columns="strategy", values="residual_mse")
    cosine = pair_mean.pivot(index="pair", columns="strategy", values="full_effect_cosine")
    base = mse[comparator].to_numpy(float)
    trial = mse[candidate].to_numpy(float)
    difference = base - trial
    improvement = float(difference.mean() / max(base.mean(), EPS))
    rng = np.random.default_rng(bootstrap_seed)
    selected = rng.integers(0, len(base), size=(draws, len(base)))
    sampled_base = base[selected].mean(axis=1)
    sampled_trial = trial[selected].mean(axis=1)
    sampled = (sampled_base - sampled_trial) / np.maximum(sampled_base, EPS)
    run_mean = targets.groupby(["data_seed", "model_seed", "strategy"], as_index=False).agg(
        residual_mse=("residual_mse", "mean")
    ).pivot(index=["data_seed", "model_seed"], columns="strategy", values="residual_mse")
    data_mean = targets.groupby(["data_seed", "strategy"], as_index=False).agg(
        residual_mse=("residual_mse", "mean")
    ).pivot(index="data_seed", columns="strategy", values="residual_mse")
    by_data = (data_mean[comparator] - data_mean[candidate]) / data_mean[comparator]
    return {
        "candidate": candidate,
        "comparator": comparator,
        "pairs": int(len(base)),
        "relative_residual_mse_improvement": improvement,
        "pair_bootstrap_ci95_low": float(np.quantile(sampled, 0.025)),
        "pair_bootstrap_ci95_high": float(np.quantile(sampled, 0.975)),
        "exact_pair_signflip_p_one_sided": exact_signflip_p(difference),
        "pair_win_rate": float(np.mean(trial < base)),
        "run_win_rate": float(np.mean(run_mean[candidate] < run_mean[comparator])),
        "data_seed_min_improvement": float(by_data.min()),
        "data_seed_max_improvement": float(by_data.max()),
        "full_effect_cosine_delta": float(np.mean(cosine[candidate] - cosine[comparator])),
    }


def strategy_table(targets: pd.DataFrame, draws: int = 200_000) -> pd.DataFrame:
    pair_mean = targets.groupby(["pair", "strategy"], as_index=False).agg(
        residual_mse=("residual_mse", "mean"),
        residual_cosine=("residual_cosine", "mean"),
        full_effect_cosine=("full_effect_cosine", "mean"),
    )
    rows: list[dict] = []
    rng = np.random.default_rng(2026080902)
    for strategy in STRATEGIES:
        local = pair_mean[pair_mean.strategy == strategy].set_index("pair").sort_index()
        selected = rng.integers(0, len(local), size=(draws, len(local)))
        row = {"strategy": strategy, "pairs": int(len(local))}
        for metric in ("residual_mse", "residual_cosine", "full_effect_cosine"):
            values = local[metric].to_numpy(float)
            sampled = values[selected].mean(axis=1)
            row[metric] = float(values.mean())
            row[f"{metric}_ci95_low"] = float(np.quantile(sampled, 0.025))
            row[f"{metric}_ci95_high"] = float(np.quantile(sampled, 0.975))
        rows.append(row)
    return pd.DataFrame(rows)


def check_engineering_audits(raw_root: Path) -> dict:
    gears_audits = sorted((raw_root / "work/formal").glob("gears_*/training_audit.json"))
    cpa_audits = sorted((raw_root / "work/formal").glob("cpa_*/training_audit.json"))
    query_audits = sorted((raw_root / "predictions/formal").glob("CPA_*_query_audit.json"))
    logs = sorted((raw_root / "logs/formal").glob("*.log"))
    assert len(gears_audits) == 9 and len(cpa_audits) == 9
    assert len(query_audits) == 9 and len(logs) == 18
    for path in gears_audits:
        audit = json.loads(path.read_text())
        assert audit["status"] == "PASS_GEARS_FORMAL_FIT"
        assert audit["epochs"] == 40
        assert audit["target_expression_used"] is False
        assert audit["test_rows_are_copied_source_controls"] is True
    for path in cpa_audits:
        audit = json.loads(path.read_text())
        assert audit["requested_max_epochs"] == 80
        assert audit["actual_epochs_completed"] == 80
        assert audit["full_schedule_complete"] is True
        assert audit["validation_schedule_complete"] is True
        assert audit["postwarm_validation_checks"] == 4
    for path in query_audits:
        audit = json.loads(path.read_text())
        assert audit["status"] == "PASS_CPA_FORMAL_FIT"
        assert audit["test_expression_used"] is False
        assert audit["test_condition_cell_count_used"] is False
        assert audit["counterfactual_registry_reencoded"] is True
        assert audit["prediction_geometry"]["collapsed"] is False
    return {
        "gears_training_audits": len(gears_audits),
        "cpa_training_audits": len(cpa_audits),
        "cpa_query_audits": len(query_audits),
        "training_logs": len(logs),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    experiment = args.experiment_root.resolve()
    raw = experiment / "raw_results"
    protocol_path = experiment / "FROZEN_SCIENTIFIC_PROTOCOL.json"
    budget_path = experiment / "code/config/FROZEN_TRAINING_BUDGET.json"
    protocol_hash = sha256(protocol_path)
    budget_hash = sha256(budget_path)
    assert protocol_hash == "459bb2b1663a78fefb4231d101941d2fe54a41c6ad0ceca1198a6d1dbe9e66b0"
    assert budget_hash == "6d64b639a60405df0cd2fbd0cc2f6869df256d4b22d95552abaeb60e8ebda7c7"

    engineering = check_engineering_audits(raw)
    frames: list[pd.DataFrame] = []
    run_metadata: list[dict] = []
    for data_seed in DATA_SEEDS:
        for model_seed in MODEL_SEEDS:
            frame, metadata = reconstruct_run(raw, data_seed, model_seed, protocol_hash, budget_hash)
            frames.append(frame)
            run_metadata.append(metadata)
    independent = pd.concat(frames, ignore_index=True)
    packaged = pd.read_csv(raw / "results/formal_summary/all_target_rows.csv")
    max_delta = compare_rows(independent, packaged)

    frozen = [
        paired_comparison(independent, candidate, comparator, 20260820 + index, 20_000)
        for index, (candidate, comparator) in enumerate(FROZEN_COMPARISONS)
    ]
    conservative = [
        paired_comparison(independent, candidate, comparator, 2026080903 + index, 200_000)
        for index, (candidate, comparator) in enumerate(CONSERVATIVE_COMPARISONS)
    ]
    strategy = strategy_table(independent)
    metadata_frame = pd.DataFrame(run_metadata)
    summary = {
        "status": "PASS_INDEPENDENT_AUDIT",
        "protocol_sha256": protocol_hash,
        "training_budget_sha256": budget_hash,
        "named_gpu_fits": 18,
        "paired_scoring_runs": 9,
        "biological_target_pairs": 14,
        "independently_reconstructed_target_rows": int(len(independent)),
        "max_abs_metric_delta_vs_packaged": max_delta,
        "engineering_audits": engineering,
        "fusion_weights": {
            "lambda_gears_mean": float(metadata_frame.lambda_gears.mean()),
            "lambda_gears_min": float(metadata_frame.lambda_gears.min()),
            "lambda_gears_max": float(metadata_frame.lambda_gears.max()),
            "lambda_cpa_mean": float(metadata_frame.lambda_cpa.mean()),
            "lambda_cpa_min": float(metadata_frame.lambda_cpa.min()),
            "lambda_cpa_max": float(metadata_frame.lambda_cpa.max()),
            "gamma_witness_mean": float(metadata_frame.gamma_witness.mean()),
        },
        "frozen_comparisons": frozen,
        "conservative_comparisons": conservative,
        "claim_boundary": (
            "Inference is over 14 strict-seen1 biological target pairs sharing FDPS. "
            "Cells, data seeds, and model seeds are not independent biological replicates."
        ),
    }

    if args.check_only:
        assert all(row["relative_residual_mse_improvement"] > 0 for row in frozen)
        assert conservative[0]["pair_win_rate"] == 1.0
        assert conservative[3]["pair_bootstrap_ci95_low"] > 0.0
        assert conservative[4]["pair_bootstrap_ci95_low"] < 0.0 < conservative[4]["pair_bootstrap_ci95_high"]
        print(json.dumps(summary, indent=2))
        return

    args.out.mkdir(parents=True, exist_ok=True)
    independent.to_csv(args.out / "independent_target_rows.csv", index=False)
    strategy.to_csv(args.out / "strategy_metrics.csv", index=False)
    pd.DataFrame(frozen).to_csv(args.out / "frozen_comparisons.csv", index=False)
    pd.DataFrame(conservative).to_csv(args.out / "conservative_comparisons.csv", index=False)
    metadata_frame.to_csv(args.out / "run_metadata.csv", index=False)
    (args.out / "audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
