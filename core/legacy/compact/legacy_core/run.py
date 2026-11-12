#!/usr/bin/env python3
"""Replay exact E17/E18/E19/E24 score inputs without restoring raw RNA arrays."""
from pathlib import Path
import argparse
import gzip
import hashlib
import importlib.util
import json
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
TOL = 1e-12
sys.dont_write_bytecode = True


def read_csv(path):
    return pd.read_csv(path, float_precision="round_trip")


def compare_tables(actual, expected, keys):
    a = actual.set_index(keys).sort_index()
    b = expected.set_index(keys).sort_index()
    assert a.index.equals(b.index), (keys, "row identities differ")
    assert set(a.columns) == set(b.columns), "columns differ"
    errors = {}
    for col in b.columns:
        if pd.api.types.is_numeric_dtype(b[col]):
            av, bv = a[col].to_numpy(float), b[col].to_numpy(float)
            assert np.array_equal(np.isnan(av), np.isnan(bv)), col
            valid = np.isfinite(av) & np.isfinite(bv)
            errors[col] = float(np.max(np.abs(av[valid] - bv[valid]))) if valid.any() else 0.
            assert errors[col] <= TOL, (col, errors[col])
        else:
            assert a[col].fillna("<NULL>").equals(b[col].fillna("<NULL>")), col
    return {"rows": len(actual), "max_abs_difference": max(errors.values(), default=0.), "columns": errors}


def verify_source_bytes():
    manifest = json.loads((ROOT / "SOURCE_MAP.json").read_text())
    for row in manifest["files"]:
        raw = (ROOT / row["compact_path"]).read_bytes()
        if row["transformation"] == "gzip_lossless":
            raw = gzip.decompress(raw)
        assert len(raw) == row["source_bytes"], row["compact_path"]
        assert hashlib.sha256(raw).hexdigest() == row["source_sha256"], row["compact_path"]
    return {"status": "PASS", "source_exact_files": len(manifest["files"])}


def paired(label, out):
    stem = ROOT / "data/ablations" / label
    f = read_csv(stem / "paired_by_condition.csv.gz")
    ref = read_csv(stem / "paired_summary.csv.gz")
    old, new = ("v13", "v14") if label == "released_vs_prior" else ("v1", "v2")
    rows = []
    for scope in ref.scope:
        sub = f if scope == "all" else f[f.dataset == scope] if scope in set(f.dataset) else f[f.subgroup == scope]
        row = {"scope": scope, "condition_seed_units": len(sub)}
        for metric in ["all_gene_mse", "top100_mse", "top100_pcc"]:
            for version in [old, new]:
                col = metric + "_" + version
                row[col] = float(sub[col].mean())
        for metric in ["all_gene_mse", "top100_mse"]:
            row["relative_" + metric + "_gain"] = float(1 - sub[metric + "_" + new].sum() / sub[metric + "_" + old].sum())
        row["top100_pcc_delta"] = float(sub.top100_pcc_delta.mean())
        rows.append(row)
    result = pd.DataFrame(rows)[ref.columns]
    check = compare_tables(result, ref, ["scope"])
    result.to_csv(out / (label + "_paired_summary.csv"), index=False)
    check["condition_seed_units"] = len(f)
    check["subgroups"] = f.subgroup.value_counts().to_dict()
    check["all"] = result[result.scope == "all"].iloc[0].to_dict()
    return check


def rebuild_operation_scores(original, official):
    blocks = []
    for (_, _, panel), block in original.groupby(["DataSet", "op", "DEG"], sort=False):
        b = block.copy()
        rank_columns = []
        for metric, name, ascending in [("cor", "Rank_pcc", False), ("edistance", "Rank_edistance", True), ("mse", "Rank_mse", True), ("sym_kldiv", "Rank_sym", True)]:
            b[name] = b[metric].rank(ascending=ascending, method="min")
            rank_columns.append(name)
        for metric, name in [("edistance", "edistance_score"), ("mse", "mse_score"), ("sym_kldiv", "sym_score")]:
            b[name] = official.inverse_minmax(b[metric])
        score_columns = ["cor", "edistance_score", "mse_score", "sym_score"]
        if panel == 100:
            b["DEG_score"] = b.common_degs
            b["was_score"] = official.inverse_minmax(b.wasserstein)
            b["Rank_deg"] = b.common_degs.rank(ascending=False, method="min")
            b["Rank_was"] = b.wasserstein.rank(ascending=True, method="min")
            rank_columns += ["Rank_deg", "Rank_was"]
            score_columns += ["DEG_score", "was_score"]
        b["ac_score"] = b[score_columns].mean(axis=1)
        b["ac_rank"] = b[rank_columns].mean(axis=1)
        b["Rank"] = b.ac_rank.rank(ascending=True, method="min")
        blocks.append(b)
    return pd.concat(blocks, ignore_index=True)


def official_tables(full):
    cell = full.groupby(["DataSet", "DEG", "method"], as_index=False).agg(mean_operation_rank=("Rank", "mean"), mean_metric_rank=("ac_rank", "mean"), mean_composite_score=("ac_score", "mean"), operations=("op", "nunique"))
    cell["dataset_panel_position"] = cell.groupby(["DataSet", "DEG"]).mean_operation_rank.rank(ascending=True, method="min")
    eq = cell.groupby("method", as_index=False).agg(mean_rank_4datasets_x_2panels=("mean_operation_rank", "mean"), mean_dataset_panel_position=("dataset_panel_position", "mean"), cells_won=("dataset_panel_position", lambda v: int(np.sum(v == 1))), cells=("dataset_panel_position", "size"))
    eq["overall_position"] = eq.mean_rank_4datasets_x_2panels.rank(ascending=True, method="min")
    eq100 = cell[cell.DEG == 100].groupby("method", as_index=False).agg(mean_rank_4datasets_top100=("mean_operation_rank", "mean"), mean_dataset_position_top100=("dataset_panel_position", "mean"), datasets_won_top100=("dataset_panel_position", lambda v: int(np.sum(v == 1))))
    eq100["top100_position"] = eq100.mean_rank_4datasets_top100.rank(ascending=True, method="min")
    top = full[full.DEG == 100].groupby("method", as_index=False).agg(mean_operation_rank=("Rank", "mean"), mean_metric_rank=("ac_rank", "mean"), mean_composite_score=("ac_score", "mean"), operation_seed_units=("op", "size"))
    top["official_position"] = top.mean_operation_rank.rank(ascending=True, method="min")
    both = full.groupby("method", as_index=False).agg(mean_operation_rank=("Rank", "mean"), mean_metric_rank=("ac_rank", "mean"), mean_composite_score=("ac_score", "mean"), operation_panel_units=("op", "size"))
    both["robustness_position"] = both.mean_operation_rank.rank(ascending=True, method="min")
    return {"dataset_panel_leaderboard": (cell, ["DataSet", "DEG", "method"]), "formal_equal_weight_leaderboard": (eq, ["method"]), "dataset_equal_weight_top100_sensitivity": (eq100, ["method"]), "formal_top100_leaderboard": (top, ["method"]), "formal_both_panels_leaderboard": (both, ["method"])}


def replay_official(label, out, scorer):
    stem = ROOT / "data/official" / label
    dst = out / label
    dst.mkdir()
    f = read_csv(stem / "published_plus_witness_per_operation.csv.gz")
    rebuilt = rebuild_operation_scores(f, scorer)
    checks = {"per_operation": compare_tables(rebuilt, f, ["DataSet", "op", "method", "DEG"])}
    tables = official_tables(rebuilt)
    for name, (table, keys) in tables.items():
        checks[name] = compare_tables(table, read_csv(stem / (name + ".csv.gz")), keys)
        table.to_csv(dst / (name + ".csv"), index=False)
    scorer.write_metric_leaderboard(rebuilt, dst)
    checks["per_metric"] = compare_tables(read_csv(dst / "official_top100_per_metric_leaderboard.csv"), read_csv(stem / "official_top100_per_metric_leaderboard.csv.gz"), ["metric", "method"])
    # Exact original 20,000-replicate, dataset-stratified perturbation-cluster recipe.
    scorer.write_cluster_bootstrap(rebuilt, dst)
    boot = json.loads((dst / "paired_cluster_bootstrap.json").read_text())
    expected = json.loads((stem / "paired_cluster_bootstrap.json").read_text())
    checks["bootstrap"] = {}
    for method, vals in boot["comparisons"].items():
        differences = []
        for key, value in vals.items():
            differences.append(float(np.max(np.abs(np.asarray(value) - np.asarray(expected["comparisons"][method][key])))))
        err = max(differences)
        assert err <= TOL, (label, method, err)
        checks["bootstrap"][method] = {"max_abs_difference": err, **vals}
    checks["WitnessCell_top100"] = tables["formal_top100_leaderboard"][0].query("method == 'WitnessCell'").iloc[0].to_dict()
    return checks


def raw_rows():
    stem = ROOT / "data/official/released"
    raw = pd.concat([read_csv(p) for p in sorted((stem / "raw").glob("*.csv.gz"))], ignore_index=True)
    raw["op"] = raw.perturb.astype(str) + "_" + raw.seed.astype(str)
    raw["metric"] = raw.metric.replace({"pearson_distance": "cor"})
    table = read_csv(stem / "published_plus_witness_per_operation.csv.gz")
    table = table[table.method == "WitnessCell"]
    fields = ["cor", "mse", "edistance", "sym_kldiv", "wasserstein", "common_degs"]
    melted = table.melt(id_vars=["DataSet", "op", "method", "DEG"], value_vars=fields, var_name="metric", value_name="performance").dropna(subset=["performance"])
    return compare_tables(raw[melted.columns], melted, ["DataSet", "op", "method", "DEG", "metric"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="New output directory outside this module")
    args = parser.parse_args()
    out = args.out.resolve()
    if out == ROOT or ROOT in out.parents:
        raise ValueError("Output must be outside the frozen compact module")
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    report = {"scope": "Exact score-row replay, not raw-cell preprocessing or model refitting", "source_integrity": verify_source_bytes()}
    spec = importlib.util.spec_from_file_location("official_score", ROOT / "source/aggregate_official_full.py")
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    report["ablations"] = {label: paired(label, out) for label in ["endpoint", "dual_head", "released_vs_prior", "released_vs_backbone"]}
    report["released_raw_metric_rows"] = raw_rows()
    report["official"] = {label: replay_official(label, out, scorer) for label in ["endpoint", "dual_head", "released"]}
    cert = json.loads((ROOT / "data/ARRAY_CERTIFICATES.json").read_text())
    assert len(cert["records"]) == 36
    assert all(r["no_truth_or_object_arrays"] for r in cert["records"])
    inactive = [r for r in cert["v14_vs_prior"] if not r["active"]]
    assert len(inactive) == 4 and all(r["inactive_exact_fallback"] and r["max_abs_prediction_delta"] == 0 for r in inactive)
    report["source_array_audit"] = {"scope": cert["scope"], "archived_records": 36, "inactive_exact_fallback_splits": 4}
    report["status"] = "PASS"
    report["elapsed_seconds"] = time.time() - started
    (out / "VALIDATION.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "out": str(out), "elapsed_seconds": report["elapsed_seconds"], "top100_released": report["official"]["released"]["WitnessCell_top100"]}, indent=2))


if __name__ == "__main__":
    main()
