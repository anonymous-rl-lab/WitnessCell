#!/usr/bin/env python3
"""Frozen-endpoint Wessels witness-library intervention; explicit fit/evaluate phases."""
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import csv
import hashlib
import json
import math
import pickle
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np

PROTOCOL_ID = "WitnessCell-step1-followup-v1-20261017"
SEEDS = (1, 2, 3)
ORDERS = (1, 2, 3)
FRACTIONS = (0.25, 0.5, 1.0)
EXPECTED = {1: (28, 63), 2: (42, 54), 3: (36, 64)}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def jd(obj):
    return json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n"


def write_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(jd(obj), encoding="utf-8")


def write_csv(path, rows):
    rows = list(rows)
    assert rows
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def canonical(condition):
    return "+".join(sorted(str(condition).split("+")))


def source_dir(args):
    return args.repo / "experiments/19_v14_incremental_amplitude_gate/src"


def input_dir(args, seed):
    return args.step1_root / "results" / f"Wessels_seed{seed}"


def load_allowed(args, seed):
    with np.load(input_dir(args, seed) / "allowed_input.npz", allow_pickle=False) as z:
        z = {k: z[k] for k in z.files}
    names = z["conditions"].astype(str).tolist()
    train = z["train"].astype(str).tolist()
    val = z["validation"].astype(str).tolist()
    targets = z["targets"].astype(str).tolist()
    doubles = [c for c in train if "+" in c]
    assert (len(doubles), len(targets)) == EXPECTED[seed]
    assert set(names) == set(["control"] + train + val)
    assert not set(names) & set(targets)
    assert len({canonical(c) for c in doubles}) == len(doubles)
    return z, names, train, val, targets, doubles


def sequence(seed, order, doubles):
    def key(c):
        text = f"{PROTOCOL_ID}|wessels_witness_quantity|{seed}|{order}|{canonical(c)}"
        return hashlib.sha256(text.encode()).hexdigest(), canonical(c)
    return sorted(doubles, key=key)


def freeze(args):
    destination = args.out / "protocol.json"
    if destination.exists():
        raise RuntimeError("Protocol already frozen; refuse replacement")
    spec = {
        "protocol_id": PROTOCOL_ID,
        "experiment": "Wessels fixed-endpoint witness-library quantity intervention",
        "design": "retrospective_development_and_diagnostic",
        "seeds": list(SEEDS), "orders": list(ORDERS), "fractions": list(FRACTIONS),
        "budget_total_fit_seconds": 3600,
        "pilot": {"seed": 1, "order": 1, "fractions": [1.0, 0.25], "test_truth_decoded": False},
        "ordering": "SHA256(protocol_id|wessels_witness_quantity|seed|order|canonical_pair), ascending",
        "nested_subset_size": "floor(fraction * number_of_original_training_combinations)",
        "subset_fitting_order": "Retain original training-condition row order inside each chosen set",
        "fixed": ["original training singles and all their moments", "control", "GO", "gene axis", "original final v14 endpoint rules", "original alpha/noise/gamma", "original test queries", "final factorized backbone"],
        "endpoint_restoration": "Fit original final endpoint head once per seed using original allowed training singles only; never refit on subset. Restore full-precision fixed backbone and verify its float32 identity/tolerance against original stored factorized predictions.",
        "training_residual": "Original core convention: observed training-combination effect minus saturate(sum of frozen baseline_head(corrected=False) endpoint effects, fixed alpha).",
        "subset_predictor": "Fixed corrected final factorized backbone + fixed gamma * solve(K_SS + fixed_noise*I, K_SQ).T @ fixed_training_residual[S]; K=incidence@incidence.T/2.",
        "100pct_gate": {"all_coordinates": True, "atol": 1e-6, "rtol": 1e-6, "checks": ["original frozen final prediction", "original frozen factorized", "original shared_endpoint_additive", "fixed selected parameters"]},
        "zero_direct_witness": "No selected training combination shares either query endpoint; count each selected combination once even if sharing both endpoints.",
        "primary_loss": "MSE excluding both query target-gene coordinates; E18 full-cell condition means",
        "secondary_loss": "WMSE excluding both query target coordinates using original Step1 evaluation weights, normalized by retained weight sum",
        "fixed_comparators": ["factorized_prediction", "shared_endpoint_additive"],
        "statistics": {
            "primary": "For each seed and fraction average per-query losses across the 3 orders, then across queries; divide by same-seed fixed factorized mean. rho_f is geometric mean of the 3 seed ratios; G_f=1-rho_f; Delta=G_100-G_25=rho_25-rho_100.",
            "bootstrap": "5000 common canonical-pair cluster resamples retaining all seeds/orders, recomputing rho and Delta; never treat 9 seed-order runs as independent n",
            "primary_metric": "loss_mse", "secondary_metric": "loss_wmse",
        },
        "phase_separation": "freeze protocol/source/subset manifests; run label-free pilot; generate ALL predictions; hash-freeze complete prediction inventory; only then decode original evaluation_only.npz and attach losses",
        "inputs": [], "source_hashes": {},
    }
    paths = list(sorted(source_dir(args).glob("*.py"))) + [
        args.step1_root / "fast_amplitude_bank.py", args.step1_root / "baselines.py",
        Path(__file__).resolve(),
    ]
    for p in paths:
        spec["source_hashes"][str(p.resolve())] = digest(p)
    go = args.repo / "experiments/15_scperturbench_sota/module/data/gears_assets/gene2go_all.pkl"
    spec["go_sha256"] = digest(go)
    subsets = []
    for seed in SEEDS:
        _, _, train, val, targets, doubles = load_allowed(args, seed)
        d = input_dir(args, seed)
        spec["inputs"].append({
            "seed": seed, "allowed_input_sha256": digest(d / "allowed_input.npz"),
            "full_replay_metadata_sha256": digest(d / "fits/full_replay.json"),
            "baseline_predictions_sha256": digest(d / "baseline_predictions.npz"),
            "original_evaluation_sha256": json.loads((d / "split_manifest.json").read_text())["evaluation_only_sha256"],
            "train_combinations": len(doubles), "targets": len(targets),
        })
        for order in ORDERS:
            ordered = sequence(seed, order, doubles)
            for fraction in FRACTIONS:
                take = math.floor(fraction * len(doubles))
                chosen = set(ordered[:take])
                subsets.append({"seed": seed, "order": order, "fraction": fraction,
                                "size": take, "chosen_conditions": [c for c in doubles if c in chosen],
                                "ordered_prefix": ordered[:take]})
    write_json(args.out / "subsets.json", subsets)
    spec["subsets_sha256"] = digest(args.out / "subsets.json")
    write_json(destination, spec)
    write_json(args.out / "protocol_freeze.json", {
        "protocol_sha256": digest(destination), "frozen_unix": time.time(),
        "test_truth_decoded": False, "test_losses_calculated": False,
    })
    print(json.dumps({"event": "protocol_frozen", "sha256": digest(destination), "subsets": len(subsets)}), flush=True)


def verify_lock(args):
    spec = json.loads((args.out / "protocol.json").read_text())
    lock = json.loads((args.out / "protocol_freeze.json").read_text())
    assert digest(args.out / "protocol.json") == lock["protocol_sha256"]
    assert digest(args.out / "subsets.json") == spec["subsets_sha256"]
    for path, expected in spec["source_hashes"].items():
        assert digest(path) == expected, path
    for row in spec["inputs"]:
        d = input_dir(args, row["seed"])
        assert digest(d / "allowed_input.npz") == row["allowed_input_sha256"]
        assert digest(d / "fits/full_replay.json") == row["full_replay_metadata_sha256"]
        assert digest(d / "baseline_predictions.npz") == row["baseline_predictions_sha256"]
    return spec, lock["protocol_sha256"]


def restore_state(args, seed, config_hash):
    out = args.out / f"seed{seed}"
    state_path = out / "frozen_state.npz"
    if state_path.exists():
        meta = json.loads((out / "restoration.json").read_text())
        assert meta["config_hash"] == config_hash and meta["state_sha256"] == digest(state_path)
        with np.load(state_path, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}, meta
    started = time.monotonic()
    sys.path.insert(0, str(source_dir(args)))
    sys.path.insert(0, str(args.step1_root))
    import witnesscell_v4_core as core
    from fast_amplitude_bank import install
    install()
    z, names, train, val, targets, doubles = load_allowed(args, seed)
    genes = z["genes"].astype(str).tolist()
    means = {c: z["means"][i].astype(float) for i, c in enumerate(names)}
    variances = {c: z["variances"][i].astype(float) for i, c in enumerate(names)}
    counts = {c: int(z["counts"][i]) for i, c in enumerate(names)}
    ctrl = means["control"]
    known = {c: means[c] - ctrl for c in train if c != "control" and "+" not in c}
    go_path = args.repo / "experiments/15_scperturbench_sota/module/data/gears_assets/gene2go_all.pkl"
    spec = json.loads((args.out / "protocol.json").read_text())
    assert digest(go_path) == spec["go_sha256"]
    with go_path.open("rb") as f:
        go = pickle.load(f)
    endpoint_start = time.monotonic()
    head = core.fit_incremental_amplitude_head(
        single_effects=known, single_means={c: means[c] for c in known},
        single_variances={c: variances[c] for c in known}, single_counts={c: counts[c] for c in known},
        control_mean=ctrl, control_variance=variances["control"], control_count=counts["control"],
        genes=genes, gene2go=go,
    )
    endpoint_seconds = time.monotonic() - endpoint_start
    par = json.loads((input_dir(args, seed) / "fits/full_replay.json").read_text())
    alpha, noise, gamma = [float(par[k]) for k in ("alpha", "noise", "gamma")]
    original_manifest = args.repo / f"experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo/Wessels/seed{seed}/manifest.json"
    original_par = json.loads(original_manifest.read_text())["selected_v2"]
    assert (alpha, noise, gamma) == (original_par["alpha"], original_par["noise_ratio"], original_par["gamma"])
    gi = {g: i for i, g in enumerate(genes)}
    nodes = sorted({g for c in train + val + targets if c != "control" for g in c.split("+")})
    baseline_endpoints, final_endpoints = {}, {}
    for node in nodes:
        baseline_endpoints[node] = known[node] if node in known else head.baseline_head.predict(node, known, gi, go)
        final_endpoints[node] = known[node] if node in known else head.predict(node, known, gi, go)
    training_base = np.stack([core.saturate(np.sum([baseline_endpoints[g] for g in c.split("+")], axis=0), alpha) for c in doubles])
    residual = np.stack([means[c] - ctrl for c in doubles]) - training_base
    additive_effect = np.stack([np.sum([final_endpoints[g] for g in c.split("+")], axis=0) for c in targets])
    additive = ctrl + additive_effect
    backbone = ctrl + core.saturate(additive_effect, alpha)
    train_design = core.incidence(doubles, nodes)
    target_design = core.incidence(targets, nodes)
    full = backbone + gamma * core.kernel_predict(train_design, target_design, residual, noise)
    with np.load(input_dir(args, seed) / "baseline_predictions.npz", allow_pickle=False) as b:
        assert np.array_equal(b["conditions"].astype(str), np.array(targets))
        original_additive = b["shared_endpoint_additive"].astype(float)
    checks = {}
    for name, restored, original in [
        ("full_prediction", full, z["frozen_prediction"]),
        ("factorized_prediction", backbone, z["frozen_factorized"]),
        ("shared_endpoint_additive", additive, original_additive),
    ]:
        restored32 = restored.astype(np.float32)
        ok = bool(np.allclose(restored32, original, atol=1e-6, rtol=1e-6))
        checks[name] = {"allclose": ok, "max_abs_error": float(np.max(np.abs(restored32.astype(float) - original.astype(float)))),
                        "float32_exact": bool(np.array_equal(restored32, original.astype(np.float32)))}
        if not ok:
            raise RuntimeError(f"100% gate failed: seed={seed}, comparison={name}, checks={checks[name]}")
    out.mkdir(parents=True, exist_ok=True)
    state = dict(genes=np.array(genes), conditions=np.array(targets), training_combinations=np.array(doubles),
                 nodes=np.array(nodes), train_design=train_design, target_design=target_design,
                 residual=residual, backbone=backbone, shared_additive=additive,
                 original_final=z["frozen_prediction"], original_factorized=z["frozen_factorized"],
                 control=ctrl, alpha=np.array(alpha), noise=np.array(noise), gamma=np.array(gamma))
    np.savez_compressed(state_path, **state)
    metadata = dict(seed=seed, config_hash=config_hash, alpha=alpha, noise=noise, gamma=gamma,
                    training_singles=len(known), training_combinations=len(doubles), queries=len(targets),
                    identity_restoration_seconds=endpoint_seconds, restoration_total_seconds=time.monotonic() - started,
                    checks=checks, state_sha256=digest(state_path), test_truth_decoded=False)
    write_json(out / "restoration.json", metadata)
    print(json.dumps({"event": "100pct_gate_passed", **metadata}), flush=True)
    return state, metadata


def prediction_name(order, fraction):
    return "prediction_full.npz" if fraction == 1.0 else f"prediction_order{order}_pct{round(100 * fraction):03d}.npz"


def predict_one(args, seed, order, fraction, state, config_hash):
    out = args.out / f"seed{seed}"
    path = out / prediction_name(order, fraction)
    meta_path = path.with_suffix(".json")
    if path.exists():
        md = json.loads(meta_path.read_text())
        assert md["config_hash"] == config_hash and md["prediction_sha256"] == digest(path)
        return md
    t = time.monotonic()
    doubles = state["training_combinations"].astype(str).tolist()
    n = math.floor(fraction * len(doubles))
    selected = set(sequence(seed, order, doubles)[:n])
    indices = np.array([i for i, c in enumerate(doubles) if c in selected], dtype=int)
    assert len(indices) == n
    train_design = state["train_design"][indices]
    target_design = state["target_design"]
    train_kernel = train_design @ train_design.T / 2.0
    cross_kernel = train_design @ target_design.T / 2.0
    coeff = np.linalg.solve(train_kernel + float(state["noise"]) * np.eye(n), cross_kernel).T
    correction = float(state["gamma"]) * (coeff @ state["residual"][indices])
    prediction = (state["backbone"] + correction).astype(np.float32)
    assert np.isfinite(prediction).all()
    witness_count = (cross_kernel > 0).sum(axis=0).astype(int)
    zero = witness_count == 0
    assert np.all(correction[zero] == 0)
    if fraction == 1.0:
        assert np.allclose(prediction, state["original_final"], atol=1e-6, rtol=1e-6)
    solve_seconds = time.monotonic() - t
    np.savez_compressed(path, prediction=prediction, conditions=state["conditions"], genes=state["genes"],
                        witness_count=witness_count, correction_rms=np.sqrt(np.mean(correction ** 2, axis=1)),
                        selected_conditions=np.array([doubles[i] for i in indices]))
    md = dict(seed=seed, order="all" if fraction == 1.0 else order, fraction=fraction, witness_library_size=n,
              queries=len(prediction), zero_direct_witness_queries=int(zero.sum()), zero_direct_witness_fraction=float(zero.mean()),
              fit_seconds=solve_seconds, total_seconds=time.monotonic() - t, config_hash=config_hash,
              prediction_sha256=digest(path), test_truth_decoded=False)
    write_json(meta_path, md)
    print(json.dumps({"event": "subset_prediction_saved", **md}), flush=True)
    return md


def fit(args, pilot):
    spec, config_hash = verify_lock(args)
    assert not (args.out / "evaluation").exists(), "Prediction generation cannot follow evaluation in this run"
    if (args.out / "predictions_frozen.json").exists():
        raise RuntimeError("Predictions already frozen; refuse changes")
    t = time.monotonic()
    records = []
    for seed in (1,) if pilot else SEEDS:
        state, restore_meta = restore_state(args, seed, config_hash)
        predict_one(args, seed, 1, 1.0, state, config_hash)
        for order in (1,) if pilot else ORDERS:
            for fraction in (0.25,) if pilot else (0.25, 0.5):
                if time.monotonic() - t > spec["budget_total_fit_seconds"]:
                    raise RuntimeError("3600 second fitting budget exhausted")
                records.append(predict_one(args, seed, order, fraction, state, config_hash))
    report = {"mode": "pilot" if pilot else "formal", "elapsed_seconds": time.monotonic() - t,
              "test_truth_decoded": False, "config_hash": config_hash, "new_or_reused_predictions": len(records)}
    write_json(args.out / ("pilot_report.json" if pilot else "formal_fit_report.json"), report)
    print(json.dumps({"event": "fit_phase_finished", **report}), flush=True)


def freeze_predictions(args):
    _, config_hash = verify_lock(args)
    assert not (args.out / "evaluation").exists()
    assert not (args.out / "predictions_frozen.json").exists()
    inventory, logical = {}, []
    for seed in SEEDS:
        for name in ("frozen_state.npz", "restoration.json"):
            p = args.out / f"seed{seed}" / name
            inventory[str(p.relative_to(args.out))] = digest(p)
        for order in ORDERS:
            for fraction in FRACTIONS:
                p = args.out / f"seed{seed}" / prediction_name(order, fraction)
                md = json.loads(p.with_suffix(".json").read_text())
                assert md["prediction_sha256"] == digest(p)
                inventory[str(p.relative_to(args.out))] = digest(p)
                inventory[str(p.with_suffix(".json").relative_to(args.out))] = digest(p.with_suffix(".json"))
                logical.append(dict(seed=seed, order=order, fraction=fraction, prediction=str(p.relative_to(args.out))))
    write_json(args.out / "predictions_frozen.json", {"frozen_unix": time.time(), "config_hash": config_hash,
               "inventory": inventory, "logical_predictions": logical, "test_truth_decoded": False,
               "logical_query_predictions": sum(EXPECTED[s][1] * 9 for s in SEEDS)})
    print(json.dumps({"event": "all_predictions_frozen", "logical_predictions": len(logical),
                      "logical_query_predictions": sum(EXPECTED[s][1] * 9 for s in SEEDS)}), flush=True)


def evaluate(args):
    spec, config_hash = verify_lock(args)
    freeze = json.loads((args.out / "predictions_frozen.json").read_text())
    assert freeze["config_hash"] == config_hash
    for path, expected in freeze["inventory"].items():
        assert digest(args.out / path) == expected
    out = args.out / "evaluation"
    assert not out.exists(), "Refuse overwriting an existing evaluation"
    rows, comparators, fit_counts = [], [], []
    evaluation_started = time.time()
    for seed in SEEDS:
        with np.load(args.out / f"seed{seed}/frozen_state.npz", allow_pickle=False) as z:
            state = {k: z[k] for k in z.files}
        expected = next(r["original_evaluation_sha256"] for r in spec["inputs"] if r["seed"] == seed)
        eval_path = input_dir(args, seed) / "evaluation_only.npz"
        assert digest(eval_path) == expected
        with np.load(eval_path, allow_pickle=False) as z:
            truth, weights = z["truth"].astype(float), z["weights"].astype(float)
            assert np.array_equal(z["conditions"], state["conditions"])
            assert np.array_equal(z["genes"], state["genes"])
        genes = state["genes"].astype(str).tolist()
        gi = {g: i for i, g in enumerate(genes)}
        conditions = state["conditions"].astype(str).tolist()
        masks = np.ones(truth.shape, dtype=bool)
        for i, condition in enumerate(conditions):
            for gene in condition.split("+"):
                masks[i, gi[gene]] = False
        retained_weights = weights * masks
        denom = retained_weights.sum(axis=1)
        assert np.all(denom > 0)
        def losses(pred):
            error = (pred.astype(float) - truth) ** 2
            return (error * masks).sum(axis=1) / masks.sum(axis=1), (error * retained_weights).sum(axis=1) / denom
        for method, pred in (("factorized_prediction", state["original_factorized"]),
                             ("shared_endpoint_additive", state["shared_additive"]),
                             ("original_final_v14", state["original_final"])):
            mse, wmse = losses(pred)
            for i, condition in enumerate(conditions):
                comparators.append(dict(dataset="Wessels", seed=seed, pair=canonical(condition), condition=condition,
                                        method=method, loss_mse=float(mse[i]), loss_wmse=float(wmse[i]), config_hash=config_hash))
        for order in ORDERS:
            for fraction in FRACTIONS:
                p = args.out / f"seed{seed}" / prediction_name(order, fraction)
                md = json.loads(p.with_suffix(".json").read_text())
                with np.load(p, allow_pickle=False) as z:
                    pred, count = z["prediction"], z["witness_count"]
                    assert np.array_equal(z["conditions"], state["conditions"])
                mse, wmse = losses(pred)
                fit_counts.append(dict(seed=seed, order=order, fraction=fraction, witness_library_size=md["witness_library_size"],
                                       query_count=len(conditions), zero_direct_witness_queries=md["zero_direct_witness_queries"],
                                       zero_direct_witness_fraction=md["zero_direct_witness_fraction"], fit_seconds=md["fit_seconds"]))
                for i, condition in enumerate(conditions):
                    rows.append(dict(dataset="Wessels", seed=seed, order=order, fraction=fraction,
                                     pair=canonical(condition), condition=condition, witness_library_size=md["witness_library_size"],
                                     witness_count=int(count[i]), zero_direct_witness=int(count[i] == 0),
                                     loss_mse=float(mse[i]), loss_wmse=float(wmse[i]), prediction_path=str(p.relative_to(args.out)),
                                     prediction_row=i, config_hash=config_hash))
    out.mkdir()
    write_csv(out / "query_losses.csv", rows)
    write_csv(out / "fixed_comparator_losses.csv", comparators)
    write_csv(out / "witness_counts.csv", fit_counts)
    point = {}
    for metric in ("loss_mse", "loss_wmse"):
        ratios, by_seed = {}, []
        for fraction in FRACTIONS:
            vals = []
            for seed in SEEDS:
                numerator = float(np.mean([r[metric] for r in rows if r["seed"] == seed and r["fraction"] == fraction]))
                denominator = float(np.mean([r[metric] for r in comparators if r["seed"] == seed and r["method"] == "factorized_prediction"]))
                vals.append(numerator / denominator)
                by_seed.append(dict(seed=seed, fraction=fraction, witness_mean_loss=numerator,
                                    factorized_mean_loss=denominator, ratio=numerator / denominator))
            ratios[str(fraction)] = float(np.exp(np.mean(np.log(vals))))
        point[metric] = dict(rho=ratios, G={k: 1-v for k, v in ratios.items()},
                             Delta_G100_minus_G25=ratios["0.25"] - ratios["1.0"], seed_results=by_seed)
    write_json(out / "point_summary.json", point)
    write_json(out / "evaluation_manifest.json", dict(protocol_id=PROTOCOL_ID, config_hash=config_hash,
               prediction_freeze_sha256=digest(args.out / "predictions_frozen.json"),
               evaluation_started_unix=evaluation_started, evaluation_finished_unix=time.time(),
               rows=len(rows), fixed_comparator_rows=len(comparators), loss_contract=spec["primary_loss"],
               inference="Bootstrap statistics delegated to root using frozen query_losses and fixed_comparator_losses; point summaries only here",
               output_hashes={p.name: digest(p) for p in sorted(out.glob("*")) if p.is_file()}))
    print(json.dumps({"event": "evaluation_complete", "rows": len(rows), "point_summary": point}), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "pilot", "predict", "freeze_predictions", "evaluate"))
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--step1-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.repo, args.step1_root, args.out = [p.resolve() for p in (args.repo, args.step1_root, args.out)]
    if args.out == args.step1_root or args.repo in args.out.parents:
        raise RuntimeError("Output must not modify original Step1 or original repository")
    args.out.mkdir(parents=True, exist_ok=True)
    {"freeze": freeze, "pilot": lambda a: fit(a, True), "predict": lambda a: fit(a, False),
     "freeze_predictions": freeze_predictions, "evaluate": evaluate}[args.mode](args)


if __name__ == "__main__":
    main()
