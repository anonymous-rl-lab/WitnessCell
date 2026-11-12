#!/usr/bin/env python3
"""Locked, descriptive follow-up diagnostics for WitnessCell step 1.

This program never fits or changes a predictor or risk score. It reads only
labelled scored_queries.csv files supplied AFTER the score freeze.

1. Leave one gene endpoint out of the test queries, independently reselect
   floor(50% of remaining queries) for S and B within each dataset/seed, and
   report exp(mean_seed(log(MSE_S/MSE_B))). No new significance tests/gates.
2. Compare including-target versus excluding-target losses on exactly the SAME
   original frozen S/B acceptance identities. A metric's evaluable pool is the
   intersection of finite values under both coordinate scopes; there is no
   reselection after inspecting a loss or after removing unavailable values.

Usage: python extra_diagnostics.py --root results --out evaluation
Synthetic checks: python extra_diagnostics.py --self-test
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


PROTOCOL_ID = "WitnessCell-step1-v1-20261017"
METRICS = ("loss_mse", "loss_wmse", "loss_top100", "weighted_delta_r2")


def file_hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def finite(value):
    try:
        x = float(value)
    except (ValueError, TypeError):
        return None
    return x if math.isfinite(x) else None


def tie_hash(dataset, pair):
    return hashlib.sha256((PROTOCOL_ID+"|"+dataset+"|"+pair).encode()).hexdigest()


def selected_indices(rows, scorer, allowed=None):
    ids = list(range(len(rows))) if allowed is None else list(allowed)
    order = sorted(ids, key=lambda i: (rows[i]["risk_"+scorer],
        tie_hash(rows[i]["dataset"], rows[i]["pair"]), rows[i]["pair"]))
    return order[:len(ids)//2]


def average(values):
    values = list(values)
    return sum(values)/len(values) if values else None


def ratio(numerator, denominator):
    if numerator is None or denominator is None or denominator <= 0 or numerator < 0:
        return None
    return numerator/denominator


def geometric_mean(values):
    values = list(values)
    if not values or any(x is None or not math.isfinite(x) or x < 0 for x in values):
        return None
    if any(x == 0 for x in values):
        return 0.
    return math.exp(sum(math.log(x) for x in values)/len(values))


def load(root):
    batches, inventory = {}, []
    for path in sorted(Path(root).glob("*/scored_queries.csv")):
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            required = {"dataset", "seed", "pair", "risk_S", "risk_B"}
            required.update(METRICS)
            required.update(m+"_including_targets" for m in METRICS)
            if not required.issubset(reader.fieldnames or []):
                raise ValueError("Missing diagnostic columns {} in {}".format(
                    required-set(reader.fieldnames or []), path))
            rows = list(reader)
        if not rows:
            raise ValueError("Empty batch {}".format(path))
        keys = {(r["dataset"], int(r["seed"])) for r in rows}
        if len(keys) != 1:
            raise ValueError("More than one dataset_seed in {}".format(path))
        key = next(iter(keys))
        if key[0] not in ("Norman", "Wessels"):
            raise ValueError("Unplanned dataset {}".format(key[0]))
        if key in batches:
            raise ValueError("Duplicate dataset_seed {}".format(key))
        pairs = set()
        for row in rows:
            row["seed"] = int(row["seed"])
            pair = row["pair"]
            genes = pair.split("+")
            if len(genes) != 2 or not all(genes) or pair != "+".join(sorted(genes)) or pair in pairs:
                raise ValueError("Noncanonical/duplicate pair {} in {}".format(pair, path))
            pairs.add(pair)
            for field in ("risk_S", "risk_B"):
                row[field] = finite(row[field])
                if row[field] is None:
                    raise ValueError("Nonfinite frozen score {} in {}".format(field, path))
            for metric in METRICS:
                for field in (metric, metric+"_including_targets"):
                    row[field] = finite(row[field])
                    if metric != "weighted_delta_r2" and row[field] is not None and row[field] < 0:
                        raise ValueError("Negative loss {} in {}".format(field, path))
            if row["loss_mse"] is None:
                raise ValueError("Primary non-target MSE must be finite for all queries: {}".format(path))
        meta_path = path.with_name("evaluation_metadata.json")
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("target_genes_removed") is not True:
            raise ValueError("Missing primary non-target assertion in {}".format(meta_path))
        if meta.get("loss_truth") != "E18_full_cell_condition_means":
            raise ValueError("Unexpected truth source in {}".format(meta_path))
        batches[key] = rows
        inventory.append({"path": str(path.resolve()), "sha256": file_hash(path),
                          "dataset": key[0], "seed": key[1], "queries": len(rows)})
    if not batches:
        raise ValueError("No */scored_queries.csv below {}".format(root))
    return batches, inventory


def leave_endpoint_out(batches):
    aggregate, seed_table, summaries = [], [], {}
    for dataset in sorted({d for d, _ in batches}):
        group = sorted((s, rows) for (d, s), rows in batches.items() if d == dataset)
        all_pairs = {r["pair"] for _, rows in group for r in rows}
        genes = sorted({g for p in all_pairs for g in p.split("+")})
        dataset_rows = []
        for excluded in [None]+genes:
            seed_ratios, removed_units = [], 0
            remaining_pairs = {p for p in all_pairs if excluded is None or excluded not in p.split("+")}
            for seed, rows in group:
                allowed = [i for i, row in enumerate(rows)
                           if excluded is None or excluded not in row["pair"].split("+")]
                n = len(allowed)
                removed_units += len(rows)-n
                risks = {scorer: average(rows[i]["loss_mse"] for i in selected_indices(rows, scorer, allowed))
                         for scorer in ("S", "B")}
                q = ratio(risks["S"], risks["B"])
                seed_ratios.append(q)
                seed_table.append({"dataset": dataset, "seed": seed,
                    "excluded_endpoint": excluded if excluded is not None else "__NONE__",
                    "original_queries": len(rows), "remaining_queries": n,
                    "removed_queries": len(rows)-n, "quota_each": n//2,
                    "coverage_actual": (n//2)/n if n else None,
                    "risk_S": risks["S"], "risk_B": risks["B"], "ratio_S_B": q})
            q = geometric_mean(seed_ratios)
            row = {"dataset": dataset, "excluded_endpoint": excluded if excluded is not None else "__NONE__",
                "seeds_required": len(group), "seeds_estimable": sum(x is not None for x in seed_ratios),
                "pair_clusters_original": len(all_pairs), "pair_clusters_remaining": len(remaining_pairs),
                "query_seed_units_removed": removed_units, "ratio_S_B": q,
                "gain": 1-q if q is not None else None}
            dataset_rows.append(row)
        baseline = dataset_rows[0]["ratio_S_B"]
        for row in dataset_rows:
            q = row["ratio_S_B"]
            row["ratio_change_from_no_exclusion"] = q-baseline if q is not None and baseline is not None else None
        evaluable = [r for r in dataset_rows[1:] if r["ratio_S_B"] is not None]
        summaries[dataset] = {"baseline_ratio": baseline, "endpoints": len(genes),
            "evaluable_exclusions": len(evaluable), "undefined_exclusions": len(genes)-len(evaluable),
            "minimum_leave_endpoint_out_ratio": min((r["ratio_S_B"] for r in evaluable), default=None),
            "maximum_leave_endpoint_out_ratio": max((r["ratio_S_B"] for r in evaluable), default=None),
            "exclusions_with_S_ratio_at_least_one": sum(r["ratio_S_B"] >= 1 for r in evaluable),
            "interpretation": "descriptive sensitivity only; no refitting, confidence interval, significance gate, or best-endpoint selection"}
        aggregate.extend(dataset_rows)
    return aggregate, seed_table, summaries


def target_coordinate_comparison(batches):
    seed_table, aggregate = [], []
    for (dataset, seed), rows in sorted(batches.items()):
        policies = {s: selected_indices(rows, s) for s in ("S", "B")}
        policies["all"] = list(range(len(rows)))
        for metric in METRICS:
            full = metric+"_including_targets"
            # The same evaluable query identities under both coordinate scopes.
            valid = {i for i, r in enumerate(rows) if r[metric] is not None and r[full] is not None}
            for policy, fixed_ids in policies.items():
                ids = [i for i in fixed_ids if i in valid]
                without = average(rows[i][metric] for i in ids)
                with_targets = average(rows[i][full] for i in ids)
                seed_table.append({"dataset": dataset, "seed": seed, "metric": metric, "policy": policy,
                    "total_queries": len(rows), "fixed_selected_n": len(fixed_ids),
                    "evaluable_selected_n": len(ids), "excluded_selected_n": len(fixed_ids)-len(ids),
                    "common_evaluable_pool_n": len(valid), "without_targets_mean": without,
                    "including_targets_mean": with_targets,
                    "including_minus_without": with_targets-without if without is not None and with_targets is not None else None})
    for dataset in sorted({d for d, _ in batches}):
        seeds = sorted(s for d, s in batches if d == dataset)
        for metric in METRICS:
            records = [r for r in seed_table if r["dataset"] == dataset and r["metric"] == metric]
            for scope in ("without_targets", "including_targets"):
                col = scope+"_mean"
                by_policy = {policy: [next(r[col] for r in records if r["seed"] == seed and r["policy"] == policy)
                    for seed in seeds] for policy in ("S", "B", "all")}
                means = {p: average(vals) if all(x is not None for x in vals) else None
                         for p, vals in by_policy.items()}
                differences = [a-b if a is not None and b is not None else None
                               for a, b in zip(by_policy["S"], by_policy["B"])]
                delta = average(differences) if all(x is not None for x in differences) else None
                ratios = [ratio(a,b) for a,b in zip(by_policy["S"], by_policy["B"])]
                # R2 can be negative: a geometric ratio is not a meaningful summary.
                q = geometric_mean(ratios) if metric != "weighted_delta_r2" else None
                aggregate.append({"dataset": dataset, "metric": metric, "coordinate_scope": scope,
                    "seeds_required": len(seeds), "mean_seed_S": means["S"], "mean_seed_B": means["B"],
                    "mean_seed_all": means["all"], "mean_seed_S_minus_B": delta,
                    "geometric_seed_ratio_S_B": q,
                    "seeds_with_unequal_evaluable_S_B_counts": sum(
                        next(r["evaluable_selected_n"] for r in records if r["seed"] == s and r["policy"] == "S") !=
                        next(r["evaluable_selected_n"] for r in records if r["seed"] == s and r["policy"] == "B") for s in seeds),
                    "selection_rule": "original full-batch S/B identities, unchanged across losses and coordinate scopes"})
    return seed_table, aggregate


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def self_test():
    rows = []
    for i, pair in enumerate(["A+B", "A+C", "B+C", "D+E", "D+F", "E+F"]):
        row = {"dataset": "Norman", "seed": 1, "pair": pair, "risk_S": float(i), "risk_B": float(i)}
        for metric in METRICS:
            row[metric] = float(i+1)
            row[metric+"_including_targets"] = 2*float(i+1)
        rows.append(row)
    batches = {("Norman", 1): rows}
    aggregate, details, _ = leave_endpoint_out(batches)
    assert all(r["ratio_S_B"] == 1 for r in aggregate)
    assert all(r["quota_each"] == r["remaining_queries"]//2 for r in details)
    by_seed, summary = target_coordinate_comparison(batches)
    assert all(r["geometric_seed_ratio_S_B"] == 1 for r in summary if r["metric"] != "weighted_delta_r2")
    assert all(r["including_targets_mean"] == 2*r["without_targets_mean"] for r in by_seed)
    # Missing optional losses do not change selection identities or primary MSE.
    rows[0]["loss_wmse_including_targets"] = None
    by_seed, _ = target_coordinate_comparison(batches)
    row = next(r for r in by_seed if r["policy"] == "S" and r["metric"] == "loss_wmse")
    assert row["fixed_selected_n"] == 3 and row["evaluable_selected_n"] == 2
    print("SYNTHETIC SELF-TEST PASS: endpoint removal, equal quota, fixed identities, paired coordinate masks")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(); return
    if args.root is None or args.out is None:
        parser.error("--root and --out are required")
    batches, inventory = load(args.root)
    endpoint_aggregate, endpoint_seeds, endpoint_summary = leave_endpoint_out(batches)
    coordinate_seeds, coordinate_summary = target_coordinate_comparison(batches)
    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out/"leave_endpoint_out_summary.csv", endpoint_aggregate)
    write_csv(args.out/"leave_endpoint_out_seed_details.csv", endpoint_seeds)
    write_csv(args.out/"target_coordinate_seed_details.csv", coordinate_seeds)
    write_csv(args.out/"target_coordinate_summary.csv", coordinate_summary)
    report = {"protocol_id": PROTOCOL_ID, "analysis_role": "prespecified_descriptive_diagnostics_no_new_significance_gate",
        "source_sha256": file_hash(__file__), "inputs": inventory,
        "endpoint_sensitivity": endpoint_summary, "target_coordinate_summary": coordinate_summary,
        "limitations": ["No predictor or score refitting", "Original acceptance identities are fixed for coordinate comparisons",
            "Only endpoint exclusion reselects half of remaining queries",
            "Missing optional metrics reduce their paired evaluable sets and may change evaluable S/B counts; counts are reported",
            "R2 is summarized by arithmetic seed means and S-minus-B, not a ratio",
            "These diagnostics do not introduce new significance tests or selection gates"]}
    (args.out/"extra_diagnostics.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print("EXTRA DIAGNOSTICS: {}".format((args.out/"extra_diagnostics.json").resolve()))


if __name__ == "__main__":
    main()
