#!/usr/bin/env python3
"""Independent aggregation for the locked WitnessCell step-1 experiment.

Does not fit a predictor, fit a risk score, or read expression matrices. Inputs:
  ROOT/{dataset}_seed{seed}/scored_queries.csv
  ROOT/{dataset}_seed{seed}/evaluation_metadata.json

Required CSV fields: dataset, seed, pair, risk_B, risk_S, effect_score,
variance_score, seen, loss_mse, true_effect. Optional loss fields are loss_wmse,
loss_relative, loss_cosine, loss_top100. `pair` must be sorted "GENE1+GENE2".
All loss columns describe the frozen FINAL predictor; loss_mse excludes both
target coordinates. Scores must be frozen before attaching test outcome labels.
The caller verifies that separation; this program records input hashes but cannot
infer temporal separation from an already labelled CSV.

Required metadata: target_genes_removed=true,
loss_truth="E18_full_cell_condition_means", oof_effect_quartile_edges=[...],
oof_effect_quartile_axis="effect_score". Quartiles are TRAIN-OOF cutpoints,
not cutpoints estimated from test outcomes. `condition` may retain the old name.

Usage: python evaluate.py --root results --out evaluation
Synthetic checks only: python evaluate.py --self-test
Requires Python >=3.8 and NumPy. No SciPy, pandas, pickle or model imports.
"""

import argparse
import csv
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


PROTOCOL_ID = "WitnessCell-step1-v1-20261017"
EXPECTED = {"Norman": {1: 73, 2: 66, 3: 69},
            "Wessels": {1: 63, 2: 54, 3: 64}}
SCORERS = {"S": "risk_S", "B": "risk_B",
           "effect": "effect_score", "variance": "variance_score"}
LOSS_FIELDS = ("loss_mse", "loss_wmse", "loss_relative", "loss_cosine", "loss_top100")
COVERAGES = (0.5, 0.3, 0.7, 0.9)


def quota_count(n, coverage):
    """Exact decimal coverage: avoid floor(0.7*90) incorrectly returning 62."""
    fraction = Fraction(str(coverage))
    return int(n)*fraction.numerator//fraction.denominator


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tie_hash(dataset, pair):
    return hashlib.sha256((PROTOCOL_ID + "|" + dataset + "|" + pair).encode()).hexdigest()


def finite_float(value, field, path):
    try:
        x = float(value)
    except (ValueError, TypeError):
        raise ValueError("Invalid {} in {}: {!r}".format(field, path, value))
    if not math.isfinite(x):
        raise ValueError("Nonfinite {} in {}".format(field, path))
    return x


def load_batches(root, allow_partial=False):
    files = sorted(Path(root).glob("*/scored_queries.csv"))
    if not files:
        raise ValueError("No */scored_queries.csv below {}".format(root))
    batches, inventory = {}, []
    for path in files:
        meta_path = path.with_name("evaluation_metadata.json")
        if not meta_path.is_file():
            raise ValueError("Missing metadata: {}".format(meta_path))
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("target_genes_removed") is not True:
            raise ValueError("Primary non-target MSE not asserted in {}".format(meta_path))
        if meta.get("loss_truth") != "E18_full_cell_condition_means":
            raise ValueError("Wrong loss truth in {}".format(meta_path))
        if meta.get("oof_effect_quartile_axis") != "effect_score":
            raise ValueError("OOF quartile axis must be effect_score: {}".format(meta_path))
        cuts = np.asarray(meta.get("oof_effect_quartile_edges", []), dtype=float)
        if cuts.shape != (3,) or not np.isfinite(cuts).all() or np.any(np.diff(cuts) < 0):
            raise ValueError("Expected three nondecreasing TRAIN-OOF cutpoints: {}".format(meta_path))
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            required = {"dataset", "seed", "pair", "seen", "loss_mse", "true_effect"}
            required.update(SCORERS.values())
            if not required.issubset(reader.fieldnames or []):
                raise ValueError("Missing fields {}: {}".format(required - set(reader.fieldnames or []), path))
            rows = list(reader)
        if not rows:
            raise ValueError("Empty batch: {}".format(path))
        identities = {(r["dataset"], int(r["seed"])) for r in rows}
        if len(identities) != 1:
            raise ValueError("CSV must contain one dataset_seed: {}".format(path))
        dataset, seed = next(iter(identities))
        if dataset not in EXPECTED:
            raise ValueError("Only prespecified Norman and Wessels supported; got {}".format(dataset))
        if (dataset, seed) in batches:
            raise ValueError("Duplicate dataset_seed {} {}".format(dataset, seed))
        pairs = [r["pair"] for r in rows]
        for pair in pairs:
            genes = pair.split("+")
            if len(genes) != 2 or not all(genes) or pair != "+".join(sorted(genes)):
                raise ValueError("Noncanonical double-perturbation pair {!r}".format(pair))
        if len(set(pairs)) != len(pairs):
            raise ValueError("Duplicate pair within a batch: {}".format(path))
        if not allow_partial and len(rows) != EXPECTED[dataset].get(seed):
            raise ValueError("Unexpected query count for {} seed {}: {}".format(dataset, seed, len(rows)))
        values = {}
        for field in tuple(SCORERS.values()) + ("true_effect",) + LOSS_FIELDS:
            present = field in rows[0]
            if not present and field in LOSS_FIELDS and field != "loss_mse":
                continue
            if field in LOSS_FIELDS and field != "loss_mse":
                optional = []
                for r in rows:
                    try:
                        x = float(r[field])
                    except (ValueError, TypeError):
                        x = float("nan")
                    optional.append(x if math.isfinite(x) else float("nan"))
                values[field] = np.asarray(optional)
            else:
                values[field] = np.array([finite_float(r[field], field, path) for r in rows])
            if (field in LOSS_FIELDS or field == "true_effect") and np.any(values[field] < 0):
                raise ValueError("Negative loss/effect is not silently clipped: {} {}".format(field, path))
        orders = {name: np.array(sorted(range(len(rows)), key=lambda i: (
            values[field][i], tie_hash(dataset, pairs[i]), pairs[i])), dtype=int)
            for name, field in SCORERS.items()}
        batches[(dataset, seed)] = dict(dataset=dataset, seed=seed, pairs=pairs,
            seen=[r["seen"] for r in rows], rows=rows, values=values, orders=orders,
            metadata=meta, cuts=cuts, path=str(path.resolve()))
        inventory.append(dict(path=str(path.resolve()), sha256=sha256(path),
            metadata_path=str(meta_path.resolve()), metadata_sha256=sha256(meta_path),
            dataset=dataset, seed=seed, rows=len(rows)))
    if not allow_partial:
        required_batches = {(d, s) for d in EXPECTED for s in EXPECTED[d]}
        if set(batches) != required_batches:
            raise ValueError("Formal aggregation requires all six batches; use --allow-partial only for diagnostics")
    return batches, inventory


def quota_weights(sorted_multiplicities, coverage):
    """Multiplicity-aware exact quotas, allowing partial boundary cluster weight.

    A bootstrap cluster copied m times gets weight m. If the boundary needs w<m
    copies, its inclusion fraction is w/m. No arbitrary copy is promoted, and
    every method has the same floor(coverage * resampled batch size) quota.
    """
    m = np.asarray(sorted_multiplicities, dtype=float)
    if m.ndim == 1:
        m = m[None, :]
    if np.any(m < 0):
        raise ValueError("Negative multiplicity")
    totals = m.sum(axis=1)
    fraction = Fraction(str(coverage))
    quotas = np.floor(fraction.numerator*totals/fraction.denominator + 1e-12)
    before = np.cumsum(m, axis=1) - m
    selected = np.minimum(m, np.maximum(quotas[:, None] - before, 0.0))
    if not np.allclose(selected.sum(axis=1), quotas, rtol=0, atol=1e-9):
        raise AssertionError("Selection quota mismatch")
    return selected, quotas


def selected_risk(batch, scorer, field, coverage, multiplicities=None):
    n = len(batch["pairs"])
    m = np.ones((1, n)) if multiplicities is None else multiplicities
    order = batch["orders"][scorer]
    weights, quotas = quota_weights(np.asarray(m)[:, order], coverage)
    numerators = np.einsum("bi,i->b", weights, batch["values"][field][order])
    return np.divide(numerators, quotas, out=np.full_like(numerators, np.nan), where=quotas > 0)


def geometric_ratio(risk_s, risk_b):
    """Seed-balanced geometric mean of ratios; no invented epsilon for zero loss."""
    a, b = np.asarray(risk_s, float), np.asarray(risk_b, float)
    valid = np.all(np.isfinite(a) & np.isfinite(b) & (a >= 0) & (b > 0), axis=0)
    result = np.full(a.shape[1:] or (), np.nan, dtype=float)
    if np.ndim(result) == 0:
        with np.errstate(divide="ignore"):
            return float(np.exp(np.mean(np.log(a / b)))) if bool(valid) else float("nan")
    with np.errstate(divide="ignore"):
        result[valid] = np.exp(np.mean(np.log(a[:, valid] / b[:, valid]), axis=0))
    return result


def interval(values):
    a = np.asarray(values, float)
    good = a[np.isfinite(a)]
    return {"ci95": np.quantile(good, [0.025, 0.975]).tolist() if len(good) else None,
            "valid_replicates": int(len(good)), "invalid_replicates": int(len(a)-len(good))}


def acceptance(batch, scorer, coverage):
    n = len(batch["pairs"])
    flags = np.zeros(n, dtype=bool)
    flags[batch["orders"][scorer][:quota_count(n, coverage)]] = True
    return flags


def distribution(values):
    x = np.asarray(values, float)
    if not len(x):
        return {"n": 0, "mean": None, "median": None, "q25": None, "q75": None}
    return {"n": int(len(x)), "mean": float(x.mean()), "median": float(np.median(x)),
            "q25": float(np.quantile(x, .25)), "q75": float(np.quantile(x, .75))}


def confounding_diagnostics(batch):
    """Descriptive, not additional success tests or post-hoc endpoint selection."""
    result = {"quartiles_defined_from": "TRAIN_OOF_effect_score", "cutpoints": batch["cuts"].tolist(),
              "scorers": {}, "seen_strata": {}}
    quartile = np.searchsorted(batch["cuts"], batch["values"]["effect_score"], side="right")
    for scorer in SCORERS:
        accepted = acceptance(batch, scorer, .5)
        result["scorers"][scorer] = {
            "accepted_true_effect": distribution(batch["values"]["true_effect"][accepted]),
            "rejected_true_effect": distribution(batch["values"]["true_effect"][~accepted]),
            "predicted_effect_quartiles": []}
        for q in range(4):
            mask = quartile == q
            result["scorers"][scorer]["predicted_effect_quartiles"].append({
                "quartile": q+1, "n": int(mask.sum()), "accepted": int((mask & accepted).sum()),
                "coverage": float(accepted[mask].mean()) if mask.any() else None})
    # Within-stratum losses use the overall-batch selected set: report its actual
    # coverage and do not imply matched coverage inside a stratum.
    for seen in sorted(set(batch["seen"])):
        mask = np.array([x == seen for x in batch["seen"]])
        result["seen_strata"][seen] = {"n": int(mask.sum()), "scorers": {}}
        for scorer in SCORERS:
            accepted = acceptance(batch, scorer, .5) & mask
            result["seen_strata"][seen]["scorers"][scorer] = {
                "accepted": int(accepted.sum()), "coverage": float(accepted.sum()/mask.sum()),
                "mse": float(batch["values"]["loss_mse"][accepted].mean()) if accepted.any() else None}
    return result


def subset_batch(batch, mask):
    """Same frozen scores, with a method-common evaluable pool."""
    ids = np.flatnonzero(mask)
    old_to_new = np.full(len(mask), -1, dtype=int)
    old_to_new[ids] = np.arange(len(ids))
    out = dict(batch)
    out["pairs"] = [batch["pairs"][i] for i in ids]
    out["seen"] = [batch["seen"][i] for i in ids]
    out["values"] = {f: x[ids] for f, x in batch["values"].items()}
    out["orders"] = {s: old_to_new[order[mask[order]]] for s, order in batch["orders"].items()}
    return out


def quartile_control(group, batch_counts):
    """Within TRAIN-OOF predicted-effect quartiles, then weight by stratum n.

    Each stratum uses floor(0.5*n_h) for BOTH methods. A nonempty stratum with
    quota zero makes that seed/control replicate undefined; it is not silently
    dropped. An empty stratum carries zero population weight.
    """
    point, boot = {s: [] for s in ("S", "B")}, {s: [] for s in ("S", "B")}
    per_seed = []
    for b in group:
        seed = b["seed"]
        m = batch_counts[seed]
        total = m.sum(axis=1)
        q = np.searchsorted(b["cuts"], b["values"]["effect_score"], side="right")
        strata = []
        point_numerators = {s: 0. for s in point}
        boot_numerators = {s: np.zeros(len(m)) for s in point}
        valid = total > 0
        point_valid = True
        for h in range(4):
            mask = q == h
            n = int(mask.sum())
            quota = quota_count(n, .5)
            record = {"quartile": h+1, "n": n, "quota": quota}
            if not n:
                record["status"] = "empty_stratum_zero_weight"
                strata.append(record)
                continue
            view = subset_batch(b, mask)
            multiplicities = m[:, mask]
            nh = multiplicities.sum(axis=1)
            valid &= (nh == 0) | (np.floor(.5*nh) > 0)
            point_valid &= quota > 0
            for scorer in point:
                rp = selected_risk(view, scorer, "loss_mse", .5)[0]
                rb = selected_risk(view, scorer, "loss_mse", .5, multiplicities)
                record["risk_"+scorer] = float(rp)
                point_numerators[scorer] += n*rp
                # Empty bootstrap strata carry no weight; nonempty undefined
                # strata invalidate the replicate via `valid` above.
                boot_numerators[scorer] += np.where(nh > 0, nh*np.nan_to_num(rb, nan=0.), 0.)
            strata.append(record)
        row = {"seed": seed, "strata": strata, "estimable": bool(point_valid)}
        for scorer in point:
            value = point_numerators[scorer]/len(b["pairs"]) if point_valid else float("nan")
            row["weighted_risk_"+scorer] = value
            point[scorer].append(value)
            boot[scorer].append(np.divide(boot_numerators[scorer], total,
                out=np.full(len(m), np.nan), where=valid))
        per_seed.append(row)
    ratio = geometric_ratio(np.asarray(point["S"]), np.asarray(point["B"]))
    boot_ratio = geometric_ratio(np.asarray(boot["S"]), np.asarray(boot["B"]))
    ci = interval(1-boot_ratio)
    return {"axis": "effect_score", "cutpoints": "fixed TRAIN-OOF quartiles per seed",
        "quota": "floor(0.5*n_h) independently within each stratum",
        "weighting": "stratum population n_h; not number accepted",
        "loss": "loss_mse", "comparison": "S/B", "ratio": ratio, "gain": 1-ratio,
        "gain_ci95": ci["ci95"], "valid_bootstrap_replicates": ci["valid_replicates"],
        "invalid_bootstrap_replicates": ci["invalid_replicates"], "seed_results": per_seed}, boot_ratio


def analyze_dataset(dataset, batches, replicates, master_seed):
    group = sorted((b for (d, _), b in batches.items() if d == dataset), key=lambda b: b["seed"])
    pairs = sorted(set(p for b in group for p in b["pairs"]))
    pair_index = {p: i for i, p in enumerate(pairs)}
    # One cluster multiplicity shared by ALL seeds, methods and coverages.
    boot_seed = int(hashlib.sha256((PROTOCOL_ID+"|"+dataset+"|bootstrap|"+str(master_seed)).encode()).hexdigest()[:8], 16)
    rng = np.random.RandomState(boot_seed)
    counts = rng.multinomial(len(pairs), np.full(len(pairs), 1/len(pairs)), size=replicates)
    batch_counts = {b["seed"]: counts[:, [pair_index[p] for p in b["pairs"]]] for b in group}
    common_losses = [f for f in LOSS_FIELDS if all(f in b["values"] for b in group)]
    result = {"dataset": dataset, "seeds": [b["seed"] for b in group],
              "query_seed_units": sum(len(b["pairs"]) for b in group),
              "unique_pair_clusters": len(pairs), "bootstrap_seed": boot_seed,
              "bootstrap_replicates_requested": replicates, "available_losses": common_losses,
              "results": [], "seed_results": [], "diagnostics": {}, "evaluable_sets": {},
              "secondary_pool_rule": "per loss, common finite pool for all scorers; rerank frozen scores within that pool"}
    bootstrap_arrays = {}
    for b in group:
        result["diagnostics"][str(b["seed"])] = confounding_diagnostics(b)
    pools, field_views = {}, {}
    for field in common_losses:
        masks = [np.isfinite(b["values"][field]) for b in group]
        # Equal masks share one matrix calculation across loss columns.
        key = tuple(mask.tobytes() for mask in masks)
        if key not in pools:
            pools[key] = {"fields": [], "views": [subset_batch(b, mask) for b, mask in zip(group, masks)],
                          "counts": [batch_counts[b["seed"]][:, mask] for b, mask in zip(group, masks)]}
        pools[key]["fields"].append(field)
        field_views[field] = pools[key]["views"]
        result["evaluable_sets"][field] = [{"seed": b["seed"], "original_n": len(mask),
            "evaluable_n": int(mask.sum()), "excluded_n": int((~mask).sum()),
            "excluded_pairs": [p for p, good in zip(b["pairs"], mask) if not good]}
            for b, mask in zip(group, masks)]
    for coverage in COVERAGES:
        risks = {f: {} for f in common_losses}
        boot_risks = {f: {} for f in common_losses}
        for pool in pools.values():
            for scorer in SCORERS:
                point_by_seed, boot_by_seed = [], []
                for b, multiplicities in zip(pool["views"], pool["counts"]):
                    order = b["orders"][scorer]
                    matrix = np.column_stack([b["values"][f] for f in pool["fields"]])[order]
                    k = quota_count(len(order), coverage)
                    point_by_seed.append(matrix[:k].mean(axis=0) if k else np.full(len(pool["fields"]), np.nan))
                    weights, quota = quota_weights(multiplicities[:, order], coverage)
                    sums = np.einsum("bi,ij->bj", weights, matrix, optimize=False)
                    boot_by_seed.append(np.divide(sums, quota[:, None], out=np.full_like(sums, np.nan), where=quota[:, None] > 0))
                point_by_seed, boot_by_seed = np.asarray(point_by_seed), np.asarray(boot_by_seed)
                for j, field in enumerate(pool["fields"]):
                    risks[field][scorer] = point_by_seed[:, j]
                    boot_risks[field][scorer] = boot_by_seed[:, :, j]
        for field in common_losses:
            views = field_views[field]
            risk, boot_risk = risks[field], boot_risks[field]
            all_risk = np.array([b["values"][field].mean() if len(b["pairs"]) else np.nan for b in views])
            for i, b in enumerate(views):
                n = len(b["pairs"])
                result["seed_results"].append({"dataset": dataset, "seed": b["seed"],
                    "coverage_requested": coverage, "n": n, "original_n": len(group[i]["pairs"]),
                    "accepted": quota_count(n, coverage),
                    "coverage_actual": quota_count(n, coverage)/n if n else None,
                    "loss": field, "full_risk": float(all_risk[i]),
                    **{"risk_"+s: float(risk[s][i]) for s in SCORERS}})
            for comparator in ("B", "effect", "variance"):
                q = geometric_ratio(risk["S"], risk[comparator])
                boot_q = geometric_ratio(boot_risk["S"], boot_risk[comparator])
                gain_ci = interval(1-boot_q)
                row = {"dataset": dataset, "coverage": coverage, "loss": field,
                    "comparison": "S/"+comparator, "ratio": q, "gain": 1-q,
                    "gain_ci95": gain_ci["ci95"], "ratio_ci95": interval(boot_q)["ci95"],
                    "valid_bootstrap_replicates": gain_ci["valid_replicates"],
                    "invalid_bootstrap_replicates": gain_ci["invalid_replicates"],
                    "seed_ratios": {str(b["seed"]): float(risk["S"][i]/risk[comparator][i])
                        if risk[comparator][i] > 0 else None for i, b in enumerate(group)}}
                result["results"].append(row)
                key = "{}__c{}__{}__S_vs_{}".format(dataset, str(coverage).replace(".", "p"), field, comparator)
                bootstrap_arrays[key] = boot_q
            for scorer in SCORERS:
                result["results"].append({"dataset": dataset, "coverage": coverage, "loss": field,
                    "comparison": scorer+"/random_expectation", "ratio": geometric_ratio(risk[scorer], all_risk),
                    "interpretation": "random fixed-quota expected mean loss equals full mean loss; descriptive reference"})
    controlled, controlled_boot = quartile_control(group, batch_counts)
    result["within_effect_quartile_control"] = controlled
    bootstrap_arrays[dataset+"__within_effect_quartile_control__S_vs_B"] = controlled_boot
    return result, bootstrap_arrays


def get_result(results, dataset, comparator="B", field="loss_mse"):
    for row in results.get(dataset, {}).get("results", []):
        if row["coverage"] == .5 and row["loss"] == field and row["comparison"] == "S/"+comparator:
            return row
    return None


def gate_summary(results, allow_partial, replicates):
    out = {"scope": "retrospective risk-increment screen; not independent confirmation",
           "minimum_practical_gain": .1, "prediction_gate": "PENDING_EXTERNAL_COMPARATOR_ASSESSMENT",
           "confounding_sensitivity_review": "REQUIRED_MANUAL_REVIEW_OF_LOCKED_DIAGNOSTICS",
           "full_experiment_gate": "PENDING"}
    primary = get_result(results, "Norman")
    if allow_partial or not primary:
        out["risk_screen"] = "DESCRIPTIVE_ONLY_PARTIAL_INPUT"
        return out
    if replicates != 5000:
        out["risk_screen"] = "DESCRIPTIVE_ONLY_NONPROTOCOL_BOOTSTRAP_COUNT"
        return out
    if primary["invalid_bootstrap_replicates"] or primary["gain_ci95"] is None:
        out["risk_screen"] = "YELLOW_UNDEFINED_BOOTSTRAP_REQUIRES_REVIEW"
        return out
    point = primary["gain"]
    lo, hi = primary["gain_ci95"]
    if hi < .1:
        out["risk_screen"] = "RED_PRACTICAL_GAIN_RULED_OUT"
        out["reason"] = "Upper confidence limit below the prespecified 10% practical threshold; not a claim of zero effect"
        return out
    if not (point >= .1 and lo > 0):
        out["risk_screen"] = "YELLOW_STATISTICALLY_OR_PRACTICALLY_UNCERTAIN"
        return out
    wessels = get_result(results, "Wessels")
    other = [get_result(results, "Norman", c) for c in ("effect", "variance")]
    out["wessels_positive_direction"] = bool(wessels and math.isfinite(wessels["gain"]) and wessels["gain"] > 0)
    out["no_point_disadvantage_against_simple_scores"] = bool(all(r and math.isfinite(r["ratio"]) and r["ratio"] <= 1 for r in other))
    if out["wessels_positive_direction"] and out["no_point_disadvantage_against_simple_scores"]:
        out["risk_screen"] = "GREEN_STATISTICAL_SCREEN_SENSITIVITY_REVIEW_PENDING"
    else:
        out["risk_screen"] = "YELLOW_PRIMARY_POSITIVE_BUT_REPLICATION_OR_CONTROL_CHECK_NOT_MET"
    return out


def json_safe(value):
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def write_outputs(outdir, report, arrays, batches):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir/"evaluation.json").write_text(json.dumps(json_safe(report), ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    np.savez_compressed(str(outdir/"bootstrap_ratios.npz"), **arrays)
    with open(outdir/"selection_manifest.csv", "w", newline="", encoding="utf-8") as f:
        fields = ["dataset", "seed", "pair", "coverage_requested", "scorer", "accepted", "tie_hash"]
        writer = csv.DictWriter(f, fields); writer.writeheader()
        for (dataset, seed), batch in sorted(batches.items()):
            for coverage in COVERAGES:
                for scorer in SCORERS:
                    flags = acceptance(batch, scorer, coverage)
                    for pair, accept in zip(batch["pairs"], flags):
                        writer.writerow(dict(dataset=dataset, seed=seed, pair=pair, coverage_requested=coverage,
                                             scorer=scorer, accepted=int(accept), tie_hash=tie_hash(dataset, pair)))
    with open(outdir/"coverage_results.csv", "w", newline="", encoding="utf-8") as f:
        fields = ["dataset", "coverage", "loss", "comparison", "ratio", "gain", "gain_ci_low", "gain_ci_high", "valid_bootstrap_replicates"]
        writer = csv.DictWriter(f, fields); writer.writeheader()
        for result in report["datasets"].values():
            for row in result["results"]:
                record = {k: row.get(k) for k in fields}
                ci = row.get("gain_ci95")
                record["gain_ci_low"], record["gain_ci_high"] = ci if ci else (None, None)
                writer.writerow(json_safe(record))


def self_test():
    # Compare weighted selection against literal duplicated rows, including a
    # duplicated boundary and score ties. No real data or artifacts are read.
    rng = np.random.RandomState(16)
    n = 9
    scores = np.array([0, 0, 1, 4, 2, 7, 8, 4, 2.], float)
    order = np.array(sorted(range(n), key=lambda i: (scores[i], str(i))))
    loss = np.arange(n, dtype=float)+1
    counts = rng.multinomial(n, np.full(n, 1/n), size=30)
    for coverage in COVERAGES:
        weights, quotas = quota_weights(counts[:, order], coverage)
        for row, quota, selected in zip(counts, quotas, weights):
            expanded = np.repeat(order, row[order])[:int(quota)]
            assert np.isclose(selected @ loss[order], loss[expanded].sum())
            assert np.isclose(selected.sum(), quota)
    # Exact point quota, including odd n; shared pair multiplicities across two
    # different seed subsets; identity scores must yield exactly zero gain.
    pairs = ["A+B", "A+C", "A+D", "B+C", "B+D"]
    counts = rng.multinomial(5, np.full(5, .2), size=100)
    batch = {"pairs": pairs, "orders": {"S": np.arange(5), "B": np.arange(5)},
             "values": {"loss_mse": np.arange(1, 6, dtype=float)}}
    rs = selected_risk(batch, "S", "loss_mse", .5, counts)
    rb = selected_risk(batch, "B", "loss_mse", .5, counts)
    assert np.array_equal(rs, rb)
    assert np.allclose(geometric_ratio(np.array([rs, rs]), np.array([rb, rb])), 1)
    assert geometric_ratio(np.array([0., 0.]), np.array([1., 2.])) == 0
    _, quota = quota_weights(np.ones((1, 5)), .5)
    assert quota[0] == 2
    assert quota_count(90, .7) == 63
    print("SYNTHETIC SELF-TEST PASS: exact quotas, duplicate boundary, ties, paired ratio")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--allow-partial", action="store_true", help="Descriptive smoke inputs only; disables a formal gate")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test(); return
    if args.root is None or args.out is None:
        parser.error("--root and --out are required")
    if args.bootstrap < 1:
        parser.error("--bootstrap must be positive")
    batches, inventory = load_batches(args.root, args.allow_partial)
    report = {"protocol_id": PROTOCOL_ID, "design": "retrospective_development_and_diagnostic",
              "primary_dataset": "Norman", "replication_dataset": "Wessels",
              "primary_comparison": "S/B at batch-wise floor(0.5*n), nontarget MSE",
              "aggregation": "exp(mean_seed(log(accepted_MSE_S/accepted_MSE_B)))",
              "bootstrap": "canonical pair clusters shared across all seeds; rerank and reselect exact quota",
              "non_target_loss_check": "required input metadata assertion; arrays evaluated upstream",
              "score_freeze_check": "caller responsibility; labelled CSV alone cannot verify timing",
              "input_inventory": inventory, "evaluator_sha256": sha256(__file__), "datasets": {}}
    arrays = {}
    for dataset in sorted(set(d for d, _ in batches)):
        result, boot = analyze_dataset(dataset, batches, args.bootstrap, args.seed)
        report["datasets"][dataset] = result; arrays.update(boot)
    report["gate"] = gate_summary(report["datasets"], args.allow_partial, args.bootstrap)
    write_outputs(args.out, report, arrays, batches)
    print(json.dumps(json_safe(report["gate"]), ensure_ascii=False, indent=2))
    print("EVALUATION REPORT: {}".format((args.out/"evaluation.json").resolve()))


if __name__ == "__main__":
    main()
