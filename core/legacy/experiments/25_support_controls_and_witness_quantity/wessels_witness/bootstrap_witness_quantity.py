#!/usr/bin/env python3
"""Common-pair bootstrap of already frozen Wessels witness intervention outputs."""
import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np

SEEDS = (1, 2, 3)
ORDERS = (1, 2, 3)
FRACTIONS = (0.25, 0.5, 1.0)
METRICS = ("loss_mse", "loss_wmse")
REFERENCES = ("factorized_prediction", "shared_endpoint_additive")
BOOTSTRAP_SEED = 2026101702
REPLICATES = 5000
PAIR_COUNT = 106
COUNTS = {1: 63, 2: 54, 3: 64}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_rows(path):
    with Path(path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows):
    rows = list(rows)
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def inputs(root):
    return [root / name for name in ("protocol.json", "protocol_freeze.json", "subsets.json",
            "predictions_frozen.json", "evaluation/query_losses.csv", "evaluation/fixed_comparator_losses.csv",
            "evaluation/witness_counts.csv", "evaluation/point_summary.json", "evaluation/evaluation_manifest.json")]


def lock(args):
    out = args.root / "statistics"
    if out.exists():
        raise RuntimeError("Statistics directory already exists; refuse replacing the analysis lock")
    plan = json.loads(args.analysis_plan.read_text())
    section = plan["wessels_witness_quantity"]
    assert section["bootstrap"]["seed"] == BOOTSTRAP_SEED
    assert section["bootstrap"]["replicates"] == REPLICATES
    assert tuple(section["seeds"]) == SEEDS and tuple(section["orders"]) == ORDERS
    assert tuple(section["fractions"]) == FRACTIONS
    out.mkdir()
    value = {
        "protocol_id": plan["protocol_id"], "analysis": "fixed-endpoint Wessels witness-library pair bootstrap",
        "locked_unix": time.time(), "source_path": str(Path(__file__).resolve()),
        "source_sha256": sha(Path(__file__).resolve()), "central_analysis_plan_path": str(args.analysis_plan),
        "central_analysis_plan_sha256": sha(args.analysis_plan),
        "inputs_sha256": {str(p.relative_to(args.root)): sha(p) for p in inputs(args.root)},
        "new_witness_experiment_point_results_previously_inspected": True,
        "historical_status": "retrospective follow-up; script/hash locked after inspecting point estimates and before calculating these intervals; not new independent analysis or prospective validation",
        "replicates": REPLICATES, "rng_seed": BOOTSTRAP_SEED, "rng": "numpy.random.default_rng, PCG64",
        "draw": "One shared rng.multinomial(106, repeat(1/106,106), size=5000) over lexically sorted canonical-pair union; retain all seed/order/fraction records for every sampled pair",
        "query_aggregation": "Average three order losses for each seed/fraction/query before applying shared pair multiplicities; then compute each seed weighted query mean",
        "seed_aggregation": "Geometric mean of the 3 seed mean-loss ratios against the same seed's fixed comparator",
        "gain": "G_fraction=1-rho_fraction",
        "primary_contrast": "Delta=G100-G25=rho25-rho100; calculate within the same bootstrap draw",
        "references": list(REFERENCES), "metrics": list(METRICS),
        "primary": {"metric": "loss_mse", "reference": "factorized_prediction", "contrast": "Delta_G100_minus_G25"},
        "interval": "Pointwise percentile 95%, numpy.percentile(method='linear'); Delta interval from draw-wise Delta, never subtracting endpoint interval limits",
        "invalid_draws": "A replicate is invalid for the relevant estimate if any seed has zero included queries, nonpositive comparator or numerator loss, or nonfinite ratios; report valid/invalid counts",
        "training_resampling": False, "orders_are_independent_replicates": False,
        "subset_namespace_execution_note": "Actual pre-evaluation experiment protocol locked SHA256(protocol_id|wessels_witness_quantity|seed|order|canonical_pair). Central plan abbreviated the namespace as dataset. Preserve the locked subsets without redrawing or changing them after results; the actual protocol/subsets hashes are recorded above.",
        "numpy_version": np.__version__, "python_version": sys.version,
    }
    write_json(out / "analysis_freeze.json", value)
    print(json.dumps({"event": "bootstrap_script_locked", "source_sha256": value["source_sha256"],
                      "analysis_freeze_sha256": sha(out / "analysis_freeze.json"), "point_results_previously_inspected": True}), flush=True)


def load_matrices(root):
    rows = read_rows(root / "evaluation/query_losses.csv")
    comparators = read_rows(root / "evaluation/fixed_comparator_losses.csv")
    assert len(rows) == 1629 and len(comparators) == 543
    pairs = sorted({r["pair"] for r in rows})
    assert len(pairs) == PAIR_COUNT
    for pair in pairs:
        assert pair == "+".join(sorted(pair.split("+"))) and len(pair.split("+")) == 2
    pi = {pair: i for i, pair in enumerate(pairs)}
    values = np.full((len(METRICS), len(SEEDS), len(ORDERS), len(FRACTIONS), len(pairs)), np.nan)
    baseline = np.full((len(METRICS), len(REFERENCES), len(SEEDS), len(pairs)), np.nan)
    presence = np.zeros((len(SEEDS), len(pairs)), dtype=bool)
    witnesses = np.full((len(SEEDS), len(ORDERS), len(FRACTIONS), len(pairs)), -1, dtype=int)
    libraries = {}
    for row in rows:
        assert row["dataset"] == "Wessels"
        s = SEEDS.index(int(row["seed"]))
        o = ORDERS.index(int(row["order"]))
        f = FRACTIONS.index(float(row["fraction"]))
        q = pi[row["pair"]]
        assert np.isnan(values[:, s, o, f, q]).all(), "Duplicate seed/order/fraction/pair"
        for m, metric in enumerate(METRICS):
            val = float(row[metric])
            assert np.isfinite(val) and val >= 0
            values[m, s, o, f, q] = val
        presence[s, q] = True
        witnesses[s, o, f, q] = int(row["witness_count"])
        assert int(row["zero_direct_witness"]) == int(int(row["witness_count"]) == 0)
        key = (s, f)
        if key in libraries:
            assert libraries[key] == int(row["witness_library_size"])
        libraries[key] = int(row["witness_library_size"])
    for row in comparators:
        assert row["dataset"] == "Wessels"
        if row["method"] not in REFERENCES:
            continue
        c = REFERENCES.index(row["method"])
        s = SEEDS.index(int(row["seed"]))
        q = pi[row["pair"]]
        assert np.isnan(baseline[:, c, s, q]).all()
        for m, metric in enumerate(METRICS):
            baseline[m, c, s, q] = float(row[metric])
    for s, seed in enumerate(SEEDS):
        assert int(presence[s].sum()) == COUNTS[seed]
        assert np.isfinite(values[:, s, :, :, :][..., presence[s]]).all()
        assert np.isfinite(baseline[:, :, s, :][..., presence[s]]).all()
        assert np.isnan(values[:, s, :, :, :][..., ~presence[s]]).all()
        assert np.isnan(baseline[:, :, s, :][..., ~presence[s]]).all()
        for o in range(len(ORDERS)):
            assert np.array_equal(values[:, s, o, 2, presence[s]], values[:, s, 0, 2, presence[s]])
    values = np.nan_to_num(values, nan=0.0)
    baseline = np.nan_to_num(baseline, nan=0.0)
    # Absent pairs are zero-filled only after asserting an identical complete query universe for every order/fraction.
    averaged = values.mean(axis=2)
    return pairs, presence, values, averaged, baseline, witnesses, libraries


def interval(values):
    values = np.asarray(values)
    good = np.isfinite(values)
    count = int(good.sum())
    ci = np.percentile(values[good], [2.5, 97.5], method="linear").tolist() if count else None
    return ci, count, int(len(values) - count)


def run(args):
    out = args.root / "statistics"
    lock_path = out / "analysis_freeze.json"
    locked = json.loads(lock_path.read_text())
    assert sha(Path(__file__).resolve()) == locked["source_sha256"]
    assert sha(args.analysis_plan) == locked["central_analysis_plan_sha256"]
    assert not (out / "summary.json").exists(), "Refuse overwriting an existing bootstrap result"
    for relative, expected in locked["inputs_sha256"].items():
        assert sha(args.root / relative) == expected, relative
    started = time.time()
    t = time.monotonic()
    pairs, presence, by_order, averaged, baseline, witnesses, libraries = load_matrices(args.root)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    assert type(rng.bit_generator).__name__ == "PCG64"
    multiplicity = rng.multinomial(PAIR_COUNT, np.full(PAIR_COUNT, 1.0 / PAIR_COUNT), size=REPLICATES)
    assert np.all(multiplicity.sum(axis=1) == PAIR_COUNT)
    sample_n = multiplicity @ presence.T
    original_n = presence.sum(axis=1)
    # rho has metric, comparator, bootstrap-draw, fraction axes.
    rho = np.full((len(METRICS), len(REFERENCES), REPLICATES, len(FRACTIONS)), np.nan)
    delta = np.full((len(METRICS), len(REFERENCES), REPLICATES), np.nan)
    fraction_rows, contrast_rows, seed_rows, order_rows = [], [], [], []
    fractions_json, contrasts_json = [], []
    point_original = json.loads((args.root / "evaluation/point_summary.json").read_text())
    point_max_difference = 0.0
    for m, metric in enumerate(METRICS):
        for c, reference in enumerate(REFERENCES):
            point_rhos = []
            boot_base = np.divide(multiplicity @ baseline[m, c].T, sample_n,
                                  out=np.full(sample_n.shape, np.nan), where=sample_n > 0)
            point_base = baseline[m, c].sum(axis=1) / original_n
            for f, fraction in enumerate(FRACTIONS):
                point_num = averaged[m, :, f].sum(axis=1) / original_n
                point_seed_ratios = point_num / point_base
                point_rho = float(np.exp(np.mean(np.log(point_seed_ratios))))
                point_rhos.append(point_rho)
                if reference == "factorized_prediction":
                    expected = point_original[metric]["rho"][str(fraction)]
                    point_max_difference = max(point_max_difference, abs(point_rho - expected))
                boot_num = np.divide(multiplicity @ averaged[m, :, f].T, sample_n,
                                     out=np.full(sample_n.shape, np.nan), where=sample_n > 0)
                valid = (np.isfinite(boot_base).all(axis=1) & np.isfinite(boot_num).all(axis=1)
                         & (boot_base > 0).all(axis=1) & (boot_num > 0).all(axis=1) & (sample_n > 0).all(axis=1))
                rho[m, c, valid, f] = np.exp(np.mean(np.log(boot_num[valid] / boot_base[valid]), axis=1))
                gain_draws = 1.0 - rho[m, c, :, f]
                ci_rho, valid_count, invalid_count = interval(rho[m, c, :, f])
                ci_gain, _, _ = interval(gain_draws)
                result = dict(metric=metric, reference=reference, fraction=fraction, rho=point_rho,
                              gain=1.0 - point_rho, rho_ci95=ci_rho, gain_ci95=ci_gain,
                              valid_bootstrap_replicates=valid_count, invalid_bootstrap_replicates=invalid_count)
                fractions_json.append(result)
                fraction_rows.append({k: v for k, v in result.items() if k not in ("rho_ci95", "gain_ci95")} |
                                     dict(rho_ci95_low=ci_rho[0], rho_ci95_high=ci_rho[1],
                                          gain_ci95_low=ci_gain[0], gain_ci95_high=ci_gain[1]))
                for s, seed in enumerate(SEEDS):
                    seed_rows.append(dict(metric=metric, reference=reference, seed=seed, fraction=fraction,
                                          mean_loss=float(point_num[s]), comparator_mean_loss=float(point_base[s]),
                                          ratio=float(point_seed_ratios[s]), gain=float(1.0 - point_seed_ratios[s]),
                                          original_query_count=int(original_n[s]), witness_library_size=libraries[(s, f)]))
            delta[m, c] = rho[m, c, :, 0] - rho[m, c, :, 2]
            ci, valid_count, invalid_count = interval(delta[m, c])
            point_delta = point_rhos[0] - point_rhos[2]
            contrast = dict(metric=metric, reference=reference, contrast="Delta_G100_minus_G25",
                            estimate=point_delta, ci95=ci, valid_bootstrap_replicates=valid_count,
                            invalid_bootstrap_replicates=invalid_count,
                            interval_construction="percentile of same-draw rho25-rho100")
            contrasts_json.append(contrast)
            contrast_rows.append({k: v for k, v in contrast.items() if k != "ci95"} |
                                 dict(ci95_low=ci[0], ci95_high=ci[1]))
            for s, seed in enumerate(SEEDS):
                for o, order in enumerate(ORDERS):
                    losses = by_order[m, s, o].sum(axis=1) / original_n[s]
                    ratios = losses / point_base[s]
                    order_rows.append(dict(metric=metric, reference=reference, seed=seed, order=order,
                                           loss25=float(losses[0]), loss50=float(losses[1]), loss100=float(losses[2]),
                                           G25=float(1-ratios[0]), G50=float(1-ratios[1]), G100=float(1-ratios[2]),
                                           Delta_G100_minus_G25=float(ratios[0]-ratios[2]),
                                           full_better_than_quarter=bool(losses[2] < losses[0]),
                                           half_better_than_quarter=bool(losses[1] < losses[0]),
                                           full_better_than_half=bool(losses[2] < losses[1])))
    assert point_max_difference < 1e-12, point_max_difference
    zero_rows, zero_summary = [], []
    for f, fraction in enumerate(FRACTIONS):
        total_zero, total_rows, seed_rates = 0, 0, []
        for s, seed in enumerate(SEEDS):
            order_rates = []
            for o, order in enumerate(ORDERS):
                counts = witnesses[s, o, f, presence[s]]
                assert (counts >= 0).all()
                nzero = int((counts == 0).sum())
                rate = float(nzero / len(counts))
                order_rates.append(rate)
                total_zero += nzero
                total_rows += len(counts)
                zero_rows.append(dict(seed=seed, order=order, fraction=fraction,
                                      witness_library_size=libraries[(s, f)], queries=len(counts),
                                      zero_incident_witness_queries=nzero, zero_incident_witness_fraction=rate,
                                      mean_incident_witness_count=float(counts.mean()), min_incident_witness_count=int(counts.min()),
                                      max_incident_witness_count=int(counts.max())))
            seed_rates.append(float(np.mean(order_rates)))
        zero_summary.append(dict(fraction=fraction, seed_equal_mean_zero_fraction=float(np.mean(seed_rates)),
                                 query_order_pooled_zero_fraction=total_zero/total_rows,
                                 zero_query_order_rows=total_zero, total_query_order_rows=total_rows,
                                 seed1_mean_zero_fraction=seed_rates[0], seed2_mean_zero_fraction=seed_rates[1],
                                 seed3_mean_zero_fraction=seed_rates[2],
                                 caveat="Order/query records describe repeated predictions; not independent sample sizes"))
    stability = []
    for metric in METRICS:
        rows = [r for r in order_rows if r["metric"] == metric and r["reference"] == "factorized_prediction"]
        seeds = []
        for seed in SEEDS:
            sub = [r for r in seed_rows if r["metric"] == metric and r["reference"] == "factorized_prediction" and r["seed"] == seed]
            losses = {r["fraction"]: r["mean_loss"] for r in sub}
            seeds.append(dict(seed=seed, full_better_than_quarter=losses[1.0] < losses[0.25],
                              monotone_mean_loss=losses[1.0] < losses[0.5] < losses[0.25]))
        stability.append(dict(metric=metric, seed_order_paths=len(rows),
                              full_better_than_quarter=sum(r["full_better_than_quarter"] for r in rows),
                              half_better_than_quarter=sum(r["half_better_than_quarter"] for r in rows),
                              full_better_than_half=sum(r["full_better_than_half"] for r in rows),
                              seeds=seeds, interpretation="Descriptive direction consistency; nine paths are not nine independent replicates"))
    write_csv(out / "fraction_results.csv", fraction_rows)
    write_csv(out / "contrasts.csv", contrast_rows)
    write_csv(out / "seed_results.csv", seed_rows)
    write_csv(out / "seed_order_stability.csv", order_rows)
    write_csv(out / "zero_witness_by_setting.csv", zero_rows)
    write_csv(out / "zero_witness_summary.csv", zero_summary)
    np.savez_compressed(out / "bootstrap_draws.npz", pairs=np.array(pairs), multiplicity=multiplicity.astype(np.uint16),
                        metrics=np.array(METRICS), references=np.array(REFERENCES), fractions=np.array(FRACTIONS),
                        rho=rho, delta=delta, seed_query_multiplicity=sample_n.astype(np.uint16))
    summary = dict(protocol_id=locked["protocol_id"], study_status=locked["historical_status"],
                   rng_seed=BOOTSTRAP_SEED, replicates=REPLICATES, pair_union_count=len(pairs),
                   original_queries_by_seed={seed: COUNTS[seed] for seed in SEEDS},
                   fraction_results=fractions_json, contrasts=contrasts_json, stability=stability,
                   zero_incident_witness_summary=zero_summary,
                   max_point_difference_from_already_seen_summary=point_max_difference,
                   interval_scope="Pointwise 95% percentile intervals; secondary metrics/references are not independent confirmations",
                   training_resampled=False, subset_orders_resampled=False,
                   primary_key="loss_mse / factorized_prediction / Delta_G100_minus_G25")
    write_json(out / "summary.json", summary)
    manifest = dict(started_unix=started, finished_unix=time.time(), elapsed_seconds=time.monotonic()-t,
                    analysis_freeze_sha256=sha(lock_path), point_results_previously_inspected=True,
                    source_sha256=sha(Path(__file__).resolve()),
                    outputs_sha256={p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "run_manifest.json"})
    write_json(out / "run_manifest.json", manifest)
    print(json.dumps({"event": "bootstrap_complete", "contrasts": contrasts_json, "stability": stability,
                      "zero_incident_witness_summary": zero_summary, "elapsed_seconds": manifest["elapsed_seconds"]}), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("lock", "run"))
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--analysis-plan", type=Path, required=True)
    args = parser.parse_args()
    args.root, args.analysis_plan = args.root.resolve(), args.analysis_plan.resolve()
    {"lock": lock, "run": run}[args.mode](args)


if __name__ == "__main__":
    main()
