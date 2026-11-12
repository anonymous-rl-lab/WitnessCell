#!/usr/bin/env python3
"""Summarize the frozen Step 1 prediction controls without reading expression data.

Inputs are only ROOT/*/baseline_losses.csv and sibling frozen_query_outputs.csv.
The upstream evaluator must already have removed both target-gene coordinates.
This module cannot establish that fact or the temporal freeze from CSVs alone.

Usage: python summarize_baselines.py --root results --out evaluation
Synthetic validation only: python summarize_baselines.py --self-test
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


PROTOCOL_ID = "WitnessCell-step1-v1-20261017"
EXPECTED = {"Norman": {1: 73, 2: 66, 3: 69},
            "Wessels": {1: 63, 2: 54, 3: 64}}
COMPARATORS = ("trainMean", "shared_endpoint_additive", "saturation_only",
               "factorized_prediction")
METHODS = ("WitnessCell",) + COMPARATORS
METRICS = ("loss_mse", "loss_wmse")
POOLS = ("all_queries", "S_accept_50pct", "seen2_all", "seen2_in_S")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tie_hash(dataset, pair):
    # Identical to evaluate.py: seed does not change a pair's tie-breaking key.
    return hashlib.sha256((PROTOCOL_ID + "|" + dataset + "|" + pair).encode()).hexdigest()


def _read_csv(path, required):
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Missing fields {sorted(missing)} in {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty input: {path}")
    return rows


def _canonical(pair):
    genes = pair.split("+")
    if len(genes) != 2 or len(set(genes)) != 2 or not all(genes):
        raise ValueError(f"Not a two-endpoint pair: {pair!r}")
    if pair != "+".join(sorted(genes)):
        raise ValueError(f"Pair is not canonical: {pair!r}")


def load_batches(root, allow_partial=False):
    """Join only the two stated CSV inputs; reject mismatched candidate pools."""
    files = sorted(Path(root).glob("*/baseline_losses.csv"))
    if not files:
        raise ValueError(f"No */baseline_losses.csv below {root}")
    batches, inventory = {}, []
    for path in files:
        score_path = path.with_name("frozen_query_outputs.csv")
        score_rows = _read_csv(score_path, {"dataset", "seed", "pair", "seen", "risk_S"})
        loss_rows = _read_csv(path, {"dataset", "seed", "pair", "seen", "method", *METRICS})
        identities = {(r["dataset"], int(r["seed"])) for r in score_rows + loss_rows}
        if len(identities) != 1:
            raise ValueError(f"Expected one matching dataset/seed across {path} and {score_path}")
        dataset, seed = next(iter(identities))
        if dataset not in EXPECTED or seed not in EXPECTED[dataset]:
            raise ValueError(f"Not a prespecified dataset/seed: {dataset}, {seed}")
        if (dataset, seed) in batches:
            raise ValueError(f"Duplicate batch: {dataset}, {seed}")
        by_pair = {}
        for row in score_rows:
            pair = row["pair"]
            _canonical(pair)
            if pair in by_pair:
                raise ValueError(f"Duplicate frozen pair in {score_path}: {pair}")
            seen = int(row["seen"])
            score = float(row["risk_S"])
            if seen not in (0, 1, 2) or not math.isfinite(score):
                raise ValueError(f"Invalid seen/risk_S for {pair} in {score_path}")
            by_pair[pair] = {"seen": seen, "score": score, "losses": {}}
        for row in loss_rows:
            pair, method = row["pair"], row["method"]
            if pair not in by_pair or method not in METHODS:
                raise ValueError(f"Unexpected pair/method {pair}/{method} in {path}")
            record = by_pair[pair]
            if int(row["seen"]) != record["seen"] or method in record["losses"]:
                raise ValueError(f"Inconsistent seen or duplicate method for {pair}/{method}")
            values = [float(row[field]) for field in METRICS]
            if any(math.isinf(value) or (math.isfinite(value) and value < 0) for value in values):
                raise ValueError(f"Infinite or negative prediction loss for {pair}/{method}")
            # NaN is retained and reported as undefined; never alter selection
            # or silently remove a method-specific row to improve a comparison.
            record["losses"][method] = values
        for pair, record in by_pair.items():
            if set(record["losses"]) != set(METHODS):
                raise ValueError(f"Incomplete five-method panel for {pair} in {path}")
        pairs = sorted(by_pair)
        if not allow_partial and len(pairs) != EXPECTED[dataset][seed]:
            raise ValueError(f"Unexpected query count for {dataset} seed {seed}: {len(pairs)}")
        scores = np.asarray([by_pair[p]["score"] for p in pairs], float)
        order = np.asarray(sorted(range(len(pairs)), key=lambda i: (
            scores[i], tie_hash(dataset, pairs[i]), pairs[i])), int)
        batches[(dataset, seed)] = {
            "dataset": dataset, "seed": seed, "pairs": pairs, "order": order,
            "seen": np.asarray([by_pair[p]["seen"] for p in pairs], int),
            "losses": np.asarray([[by_pair[p]["losses"][m] for m in METHODS] for p in pairs], float),
        }
        inventory.append({"dataset": dataset, "seed": seed, "query_count": len(pairs),
                          "baseline_losses": str(path.resolve()), "baseline_sha256": sha256(path),
                          "frozen_query_outputs": str(score_path.resolve()),
                          "frozen_query_sha256": sha256(score_path)})
    if not allow_partial:
        for dataset in {d for d, _ in batches}:
            if {s for d, s in batches if d == dataset} != {1, 2, 3}:
                raise ValueError(f"Formal {dataset} summary requires seeds 1, 2, 3")
    return batches, inventory


def selection_weights(multiplicities, order):
    """Reselect floor(n/2) copies by frozen S, allowing partial boundary weight."""
    m = np.atleast_2d(np.asarray(multiplicities, dtype=float))
    ordered = m[:, order]
    quota = np.floor(m.sum(axis=1) / 2.0)
    before = np.cumsum(ordered, axis=1) - ordered
    selected = np.minimum(ordered, np.maximum(quota[:, None] - before, 0.0))
    result = np.zeros_like(m)
    result[:, order] = selected
    if not np.array_equal(result.sum(axis=1), quota):
        raise AssertionError("S acceptance quota mismatch")
    return result


def pool_weights(batch, multiplicities, pool):
    m = np.atleast_2d(np.asarray(multiplicities, dtype=float))
    if pool in ("S_accept_50pct", "seen2_in_S"):
        m = selection_weights(m, batch["order"])
    if pool in ("seen2_all", "seen2_in_S"):
        m = m * (batch["seen"] == 2)[None, :]
    return m


def mean_losses(batch, weights):
    """Same denominator for all methods; any weighted NaN makes that risk NaN."""
    values = batch["losses"]
    finite = np.isfinite(values)
    total = weights.sum(axis=1)
    numerators = np.einsum("bi,ijk->bjk", weights, np.where(finite, values, 0.0), optimize=False)
    bad = np.einsum("bi,ijk->bjk", weights, (~finite).astype(float), optimize=False) > 0
    result = np.divide(numerators, total[:, None, None], out=np.full_like(numerators, np.nan),
                       where=total[:, None, None] > 0)
    result[bad] = np.nan
    return result


def geometric_ratio(numerator, denominator):
    """Equal seed log ratios; zero reference risk is undefined, without epsilon."""
    a, b = np.asarray(numerator, float), np.asarray(denominator, float)
    valid = np.all(np.isfinite(a) & np.isfinite(b) & (a >= 0) & (b > 0), axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.exp(np.mean(np.log(a / b), axis=0))
    return np.where(valid, ratio, np.nan)


def _interval(values):
    finite = np.asarray(values, float)
    finite = finite[np.isfinite(finite)]
    return np.quantile(finite, [.025, .975]).tolist() if len(finite) else None


def analyze_dataset(dataset, batches, replicates=5000, master_seed=1729):
    group = sorted([b for (d, _), b in batches.items() if d == dataset], key=lambda b: b["seed"])
    pairs = sorted({p for batch in group for p in batch["pairs"]})
    indices = {pair: i for i, pair in enumerate(pairs)}
    # Same deterministic cluster draws as evaluate.py. A drawn pair carries
    # all its seeds, methods and losses; seeds are not independent resamples.
    seed_text = PROTOCOL_ID + "|" + dataset + "|bootstrap|" + str(master_seed)
    bootstrap_seed = int(hashlib.sha256(seed_text.encode()).hexdigest()[:8], 16)
    rng = np.random.RandomState(bootstrap_seed)
    counts = rng.multinomial(len(pairs), np.full(len(pairs), 1 / len(pairs)), size=replicates)
    per_batch = [counts[:, [indices[p] for p in b["pairs"]]] for b in group]
    summaries, seed_results, selection_rows = [], [], []
    for batch in group:
        accepted = selection_weights(np.ones(len(batch["pairs"])), batch["order"])[0]
        for i, pair in enumerate(batch["pairs"]):
            selection_rows.append({"dataset": dataset, "seed": batch["seed"], "pair": pair,
                                   "seen": int(batch["seen"][i]), "accepted_S": int(accepted[i]),
                                   "tie_hash": tie_hash(dataset, pair)})
    for pool in POOLS:
        point_values, boot_values, pool_sizes = [], [], []
        for batch, multiplicities in zip(group, per_batch):
            point_weights = pool_weights(batch, np.ones(len(batch["pairs"])), pool)
            point = mean_losses(batch, point_weights)[0]
            boot = mean_losses(batch, pool_weights(batch, multiplicities, pool))
            point_values.append(point)
            boot_values.append(boot)
            pool_size = int(point_weights.sum())
            pool_sizes.append(pool_size)
            for j, comparator in enumerate(COMPARATORS, start=1):
                for k, metric in enumerate(METRICS):
                    wc, ref = float(point[0, k]), float(point[j, k])
                    ratio = wc / ref if math.isfinite(wc) and math.isfinite(ref) and ref > 0 else float("nan")
                    seed_results.append({"dataset": dataset, "seed": batch["seed"], "pool": pool,
                                         "metric": metric, "comparator": comparator,
                                         "query_count": len(batch["pairs"]), "pool_count": pool_size,
                                         "pool_fraction": pool_size / len(batch["pairs"]),
                                         "witnesscell_loss": wc, "comparator_loss": ref, "ratio": ratio})
        point_values, boot_values = np.asarray(point_values), np.asarray(boot_values)
        for j, comparator in enumerate(COMPARATORS, start=1):
            for k, metric in enumerate(METRICS):
                ratio = float(geometric_ratio(point_values[:, 0, k], point_values[:, j, k]))
                boot_ratio = geometric_ratio(boot_values[:, :, 0, k], boot_values[:, :, j, k])
                valid = int(np.isfinite(boot_ratio).sum())
                ci = _interval(boot_ratio)
                gain_ci = _interval(1.0 - boot_ratio)
                summaries.append({
                    "dataset": dataset, "pool": pool, "metric": metric, "comparator": comparator,
                    "baseline_definition": (
                        "measured_single_endpoint_additive" if comparator == "trainMean" and pool.startswith("seen2")
                        else "measured_or_equal_mean_completed_endpoint_additive" if comparator == "trainMean"
                        else comparator),
                    "seeds": [b["seed"] for b in group], "seed_pool_counts": pool_sizes,
                    "ratio": ratio, "gain": 1.0 - ratio, "ratio_ci95": ci, "gain_ci95": gain_ci,
                    "witnesscell_seed_equal_mean_loss": float(np.mean(point_values[:, 0, k])),
                    "comparator_seed_equal_mean_loss": float(np.mean(point_values[:, j, k])),
                    "valid_bootstrap_replicates": valid, "invalid_bootstrap_replicates": replicates - valid,
                    "status": "estimable" if math.isfinite(ratio) and valid == replicates else "undefined_values_require_review",
                })
    missing = {metric: {method: int(sum(np.count_nonzero(~np.isfinite(b["losses"][:, j, k])) for b in group))
                       for j, method in enumerate(METHODS)} for k, metric in enumerate(METRICS)}
    saturation = next(r for r in summaries if r["pool"] == "all_queries" and
                      r["metric"] == "loss_wmse" and r["comparator"] == "saturation_only")
    complete = ([b["seed"] for b in group] == [1, 2, 3]
                and replicates == 5000 and master_seed == 1729)
    simple_gate = "NOT_ESTIMABLE_OR_NONFORMAL"
    if complete and saturation["status"] == "estimable" and saturation["ratio_ci95"]:
        simple_gate = "PASS_1.05_UPPER_BOUND" if saturation["ratio_ci95"][1] <= 1.05 else "NOT_PASS_1.05_UPPER_BOUND"
    return ({"dataset": dataset, "seeds": [b["seed"] for b in group],
             "query_seed_units": sum(len(b["pairs"]) for b in group), "unique_pair_clusters": len(pairs),
             "bootstrap_seed": bootstrap_seed, "bootstrap_replicates": replicates,
             "nonfinite_loss_counts": missing, "results": summaries,
             "saturation_only_noninferiority_component": simple_gate,
             "scouter": {"status": "PENDING", "reason": "No locked condition-by-gene scouter predictions supplied"},
             "prediction_competitiveness_gate": "PENDING_SCOUTER"}, seed_results, selection_rows)


def json_safe(value):
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def _write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(json_safe(rows))


def write_outputs(out, report, seed_rows, selections):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "baseline_summary.json").write_text(
        json.dumps(json_safe(report), ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    csv_rows = []
    for dataset in report["datasets"].values():
        for row in dataset["results"]:
            flat = {k: v for k, v in row.items() if k not in {"ratio_ci95", "gain_ci95", "seeds", "seed_pool_counts"}}
            for name in ("ratio", "gain"):
                ci = row[name + "_ci95"]
                flat[name + "_ci_low"], flat[name + "_ci_high"] = ci if ci else (None, None)
            flat["seeds"] = ";".join(map(str, row["seeds"]))
            flat["seed_pool_counts"] = ";".join(map(str, row["seed_pool_counts"]))
            csv_rows.append(flat)
    _write_csv(out / "baseline_summary.csv", csv_rows)
    _write_csv(out / "baseline_seed_results.csv", seed_rows)
    _write_csv(out / "baseline_selection.csv", selections)


def self_test():
    """In-memory synthetic checks only; do not open actual result tables."""
    rng = np.random.RandomState(37)
    order = np.array([2, 0, 3, 1, 4])
    multiplicities = rng.multinomial(5, [.2] * 5, size=100)
    selected = selection_weights(multiplicities, order)
    for counts, weights in zip(multiplicities, selected):
        expanded = np.repeat(order, counts[order])[:int(counts.sum() // 2)]
        np.testing.assert_array_equal(weights, np.bincount(expanded, minlength=5))
    batches = {}
    for seed, pairs in [(1, ["A+B", "A+C", "A+D", "B+C", "B+D"]),
                        (2, ["A+B", "A+C", "A+D", "B+C"]),
                        (3, ["A+B", "A+C", "A+D", "B+D"])]:
        n = len(pairs)
        values = np.zeros((n, len(METHODS), len(METRICS)))
        for i in range(n):
            values[i, 0] = [i + seed + 1.0, (i + seed + 1.0) / 10]
            values[i, 1:] = values[i, 0] * 2
        batches[("Norman", seed)] = {"dataset": "Norman", "seed": seed, "pairs": pairs,
            "order": np.arange(n), "seen": np.array([2, 2] + [1] * (n - 2)), "losses": values}
    report, _, _ = analyze_dataset("Norman", batches, 5000)
    for row in report["results"]:
        np.testing.assert_allclose(row["ratio"], .5, rtol=0, atol=1e-15)
        if row["valid_bootstrap_replicates"]:
            np.testing.assert_allclose(row["ratio_ci95"], [.5, .5], rtol=0, atol=1e-15)
    assert report["unique_pair_clusters"] == 5
    assert report["scouter"]["status"] == "PENDING"
    assert float(geometric_ratio(np.array([0., 0.]), np.array([1., 2.]))) == 0.0
    assert math.isnan(float(geometric_ratio(np.array([0., 1.]), np.array([0., 2.]))))
    print("SYNTHETIC PASS: 5000 paired-cluster draws, common S quotas, equal-seed ratios, zero-loss handling")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--allow-partial", action="store_true", help="Diagnostic only; no formal competitiveness gate")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.root is None or args.out is None:
        parser.error("--root and --out are required")
    if args.bootstrap < 1:
        parser.error("--bootstrap must be positive")
    batches, inventory = load_batches(args.root, args.allow_partial)
    report = {
        "protocol_id": PROTOCOL_ID, "design": "retrospective_development_and_diagnostic",
        "aggregation": "exp(mean_over_three_seeds(log(mean_loss_WitnessCell / mean_loss_comparator)))",
        "ratio_direction": "less than one favors WitnessCell",
        "bootstrap": "canonical pair clusters shared across seeds and methods; S reranked and exact half quota reselected",
        "selection": "frozen risk_S ascending; ties by protocol/dataset/pair SHA256; floor(n/2) per seed",
        "seen2_in_S": "intersection of the full-batch S accepted set with seen2; no independent stratum reselection",
        "loss_contract": "upstream non-target MSE and E22-weighted WMSE using E18 full-cell condition means",
        "verification_limit": "Only two CSVs read; non-target construction and temporal score freeze are upstream responsibilities",
        "missing_value_rule": "no candidate removal or method-specific reselection; undefined weighted losses reported",
        "diagnostic_only": bool(args.allow_partial or args.bootstrap != 5000 or args.seed != 1729),
        "bootstrap_master_seed": args.seed,
        "input_inventory": inventory, "summarizer_sha256": sha256(__file__), "datasets": {},
        "prediction_competitiveness_gate": "PENDING_SCOUTER",
    }
    seed_rows, selections = [], []
    for dataset in sorted({d for d, _ in batches}):
        result, per_seed, per_query = analyze_dataset(dataset, batches, args.bootstrap, args.seed)
        if args.allow_partial:
            result["saturation_only_noninferiority_component"] = "DESCRIPTIVE_ONLY_PARTIAL_INPUT"
        report["datasets"][dataset] = result
        seed_rows.extend(per_seed)
        selections.extend(per_query)
    write_outputs(args.out, report, seed_rows, selections)
    print(f"BASELINE SUMMARY: {(args.out / 'baseline_summary.json').resolve()}")
    print("PREDICTION COMPETITIVENESS GATE: PENDING_SCOUTER")


if __name__ == "__main__":
    main()
