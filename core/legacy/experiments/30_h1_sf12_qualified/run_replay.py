"""Portable numerical replay; never trains or selects a new model.

Partner scores are reconstructed from measured targets and frozen vectors.
Full-response scores use archived per-edge losses and case_R quadratics.
SF12 and case_R coefficients are independently re-solved from saved normal
equations. Expression feature construction remains an upstream operation.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent


def replay(output):
    recorded = json.loads((ROOT / "run/SCORES.json").read_text())
    spec = json.loads((ROOT / "run/PROTOCOL.json").read_text())["spec"]
    ev = np.load(ROOT / "run/EVALUATION_STATS.npz", allow_pickle=False)
    vec = np.load(ROOT / "run/COMPACT_VECTORS.npz", allow_pickle=False)
    sf12 = np.load(ROOT / "run/SF12_SHARED_TRAIN_STATS.npz", allow_pickle=False)
    methods = [str(x) for x in vec["names"]]
    lookup = {int(e): i for i, e in enumerate(vec["unique_ids"])}
    c = vec["contrast"]
    assert np.isclose(c @ c, 1)
    incidence_contrasts = np.zeros((len(spec["tasks"]), len(vec["pairs"])))
    errors, coef_errors, quad_errors, sf12_errors, rows = [], [], [], [], []
    partner = []
    full = ev["per_edge_SSE"]
    for i, task in enumerate(spec["tasks"]):
        name = task["name"]
        assert list(vec["query_ids"][i]) == task["query_ids"]
        incidence_contrasts[i, task["query_ids"]] = c
        assert set(task["train_ids"]) == set(range(21)) - set(task["query_ids"])
        for episode in task["episodes"]:
            assert set(episode["val_ids"]) <= set(task["train_ids"])
            assert set(episode["train_ids"]) == set(task["train_ids"]) - set(episode["val_ids"])
        truth = vec["unique_truth"][[lookup[int(e)] for e in vec["query_ids"][i]]]
        truth_contrast = c @ truth
        ps = np.sum((truth_contrast[None] - vec["partner_predictions"][i]) ** 2, axis=1)
        partner.append(ps)
        # Solve the unnormalized equivalent of the archived standardized ridge.
        prefix = name + "_"
        gram, rhs = sf12[prefix + "M"], sf12[prefix + "b"]
        gram = (gram + gram.T) / 2
        scale12 = np.sqrt(np.maximum(np.diag(gram), 0))
        scale12 = np.maximum(scale12, max(float(scale12.max(initial=0)) * 1e-10, 1e-15))
        theta12 = np.linalg.solve(gram + spec["old12"]["rho"] * np.diag(scale12 * scale12), rhs)
        trace12 = json.loads((ROOT / "run" / name / "TRACE.json").read_text())["outer_old12"]
        sf12_errors.extend([float(np.max(abs(theta12 - sf12[prefix + "theta"]))), float(np.max(abs(theta12 - np.asarray(trace12["theta"]))))])
        assert len(theta12) == 12
        assert list(sf12[prefix + "train_ids"]) == task["train_ids"]
        assert list(sf12[prefix + "query_ids"]) == task["query_ids"]
        train_sse = float(sf12[prefix + "train_projected_residual_SSE"] - 2 * theta12 @ rhs + theta12 @ gram @ theta12)
        assert train_sse >= -1e-9
        z = np.load(ROOT / "run" / name / "TRAIN_STATS.npz", allow_pickle=False)
        count = len(task["episodes"])
        m = sum(z[f"e{k}_M"] for k in range(count))
        b = sum(z[f"e{k}_b"] for k in range(count))
        assert np.allclose(m, z["M"], atol=1e-15, rtol=1e-12)
        assert np.allclose(b, z["b"], atol=1e-15, rtol=1e-12)
        scale = np.sqrt(np.maximum(np.diag(m), 0))
        scale = np.maximum(scale, max(float(scale.max(initial=0)) * 1e-10, 1e-15))
        theta = np.linalg.solve(m + spec["case_R"]["rho"] * np.diag(scale * scale), b)
        coef_errors.append(float(np.max(abs(theta - z["theta"]))))
        for suffix, key in [("F", "full_SSE"), ("C", "partner_SSE")]:
            old = recorded["tasks"][name]["methods"]["old12"][key]
            expected = recorded["tasks"][name]["methods"]["case_R"][key]
            recomputed = old - 2 * float(theta @ ev[name + "_b" + suffix]) + float(theta @ ev[name + "_M" + suffix] @ theta)
            quad_errors.append(abs(recomputed - expected))
        for j, method in enumerate(methods):
            ref = recorded["tasks"][name]["methods"][method]
            fs = float(full[i, j].sum())
            pg = float(100 * (1 - ps[j] / ps[0]))
            fg = float(100 * (1 - fs / full[i, 0].sum()))
            errors.extend([abs(fs - ref["full_SSE"]), abs(float(ps[j]) - ref["partner_SSE"]), abs(pg - ref["partner_gain_vs_M0_pct"]), abs(fg - ref["full_gain_vs_M0_pct"])])
            rows.append(dict(task=name, method=method, full_SSE=fs, partner_SSE=float(ps[j]), full_gain_pct=fg, partner_gain_pct=pg))
    partner = np.stack(partner)
    aggregate = {}
    selections = {"all5": list(range(5)), "exclude_Q4": [0, 1, 2, 4]}
    for label, idx in selections.items():
        ps, fs = partner[idx].sum(0), full[idx].sum(axis=(0, 2))
        edge_rows = []
        for edge in vec["unique_ids"]:
            entries = [full[i, :, int(np.flatnonzero(vec["query_ids"][i] == edge)[0])] for i in idx if edge in vec["query_ids"][i]]
            if entries:
                edge_rows.append(np.mean(entries, axis=0))
        edge = np.mean(edge_rows, axis=0)
        aggregate[label] = {}
        for j, method in enumerate(methods):
            got = dict(pooled_partner_SSE=float(ps[j]), pooled_full_SSE=float(fs[j]), partner_gain_vs_M0_pct=float(100 * (1 - ps[j] / ps[0])), partner_gain_vs_old12_pct=float(100 * (1 - ps[j] / ps[1])), pooled_full_gain_vs_M0_pct=float(100 * (1 - fs[j] / fs[0])), edge_balanced_full_SSE=float(edge[j]), edge_balanced_full_gain_vs_M0_pct=float(100 * (1 - edge[j] / edge[0])), edge_balanced_full_gain_vs_old12_pct=float(100 * (1 - edge[j] / edge[1])))
            expected = recorded["aggregate"][label][method]
            errors.extend(abs(got[k] - expected[k]) for k in got)
            aggregate[label][method] = got
    leave_one = {}
    for i, task in enumerate(spec["tasks"]):
        ps = partner[[j for j in range(5) if j != i]].sum(0)
        leave_one[task["name"]] = {m: float(100 * (1 - ps[j] / ps[0])) for j, m in enumerate(methods)}
        errors.extend(abs(leave_one[task["name"]][m] - recorded["leave_one_task_out_partner_gain_vs_M0_pct"][task["name"]][m]) for m in methods)
    rank = int(np.linalg.matrix_rank(incidence_contrasts))
    assert rank == 3 and len(vec["unique_ids"]) == 8
    assert max(errors) < 1e-9 and max(coef_errors) < 1e-12 and max(quad_errors) < 1e-9 and max(sf12_errors) < 1e-12
    report = dict(status="PASS", score_max_abs=max(errors), SF12_coefficient_max_abs=max(sf12_errors), SF12_coefficient_solves=5, case_R_coefficient_max_abs=max(coef_errors), case_R_quadratic_max_abs=max(quad_errors), contrast_identity_rank=rank, unique_query_edges=len(vec["unique_ids"]), aggregate=aggregate, leave_one_task_out_partner_gain_vs_M0_pct=leave_one, all_task_scores=rows, evidence="Historical-panel development extension; overlapping tasks; no stable new-partner superiority claim.", scope="Partner vectors and targets, archived full per-edge losses, exact SF12 and case_R normal-equation replay. SF12 sufficient statistics were derived from the verified original fixed-setting intermediate inputs and matched all archived coefficients and partner corrections. No raw-cell processing, new tuning or feature construction during compact replay.")
    output = Path(output).resolve()
    if output == ROOT / "run" or (ROOT / "run") in output.parents:
        raise ValueError("Output must not overwrite archived run/ data")
    output.mkdir(parents=True, exist_ok=True)
    (output / "REPLAY_REPORT.json").write_text(json.dumps(report, indent=2) + "\n")
    with (output / "ALL_TASK_SCORES.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="replay_output")
    args = parser.parse_args()
    result = replay(args.output)
    print(json.dumps({k: result[k] for k in ["status", "score_max_abs", "SF12_coefficient_max_abs", "case_R_coefficient_max_abs", "case_R_quadratic_max_abs", "contrast_identity_rank", "unique_query_edges"]}, indent=2))
