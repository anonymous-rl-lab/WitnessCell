"""Recover the compact SF12 sufficient statistics from original intermediates.

No model settings are changed. The five exact fixed-setting fits must reproduce
the archived coefficients and frozen partner corrections before any output is
accepted. This command is optional; run.py needs no external inputs.
"""
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset", required=True, type=Path, help="Original TRAINING_BATCH_ONLY.npz")
    parser.add_argument("--knowledge", type=Path, default=ROOT / "data/KNOWLEDGE_NODE.npz")
    parser.add_argument("--out", required=True, type=Path, help="Fresh directory outside the module")
    args = parser.parse_args()
    output = args.out.resolve()
    if output == ROOT or ROOT in output.parents:
        parser.error("Output must be outside the module")
    if output.exists() and any(output.iterdir()):
        parser.error("Output must be absent or empty")
    manifest = json.loads((ROOT / "provenance/RECONSTRUCTION_INPUTS.json").read_text())
    expected = {r["filename"]: r["sha256"] for r in manifest["inputs"]}
    inputs = {"TRAINING_BATCH_ONLY.npz": args.asset, "KNOWLEDGE_NODE.npz": args.knowledge}
    for name, path in inputs.items():
        if not path.is_file() or sha(path) != expected[name]:
            parser.error(f"Missing or mismatched input: {name}")
    sys.path.insert(0, str(ROOT / "source"))
    import legacy_model as model
    from fold_factory import FoldFactory
    spec = json.loads((ROOT / "run/PROTOCOL.json").read_text())["spec"]
    with np.load(args.knowledge, allow_pickle=False) as z:
        knowledge = z["knowledge_node"].copy()
        knowledge_nodes = z["nodes"].copy()
    with np.load(ROOT / "run/COMPACT_VECTORS.npz", allow_pickle=False) as z:
        predicted = z["partner_predictions"].copy()
        contrast = z["contrast"].copy()
    saved, records, capture = {}, [], {}
    original_solve = model.solve

    def tracked_solve(m, b):
        capture["M"], capture["b"] = m.copy(), b.copy()
        return original_solve(m, b)

    start = time.monotonic()
    model.solve = tracked_solve
    try:
        for i, task in enumerate(spec["tasks"]):
            fac = FoldFactory(asset_path=args.asset, knowledge_node=knowledge)
            assert np.array_equal(fac.asset["nodes"], knowledge_nodes)
            query = np.asarray(task["query_ids"])
            fac.asset["batch_values"][np.isin(fac.asset["obs_condition"], len(fac.nodes) + query)] = np.nan
            fold = fac.get_fold(task["train_ids"], task["query_ids"])
            delta, meta = model.old12(fold)
            archived = json.loads((ROOT / "run" / task["name"] / "TRACE.json").read_text())["outer_old12"]
            ce = float(np.max(abs(np.asarray(meta["theta"]) - np.asarray(archived["theta"]))))
            pe = float(np.max(abs(contrast @ delta - (predicted[i, 1] - predicted[i, 0]))))
            assert ce < 1e-12 and pe < 1e-12, (ce, pe)
            assert time.monotonic() - start < 60, "Exceeded bounded derivation budget"
            prefix = task["name"] + "_"
            values = dict(M=capture["M"], b=capture["b"], theta=np.asarray(meta["theta"]), train_ids=np.asarray(task["train_ids"]), query_ids=query, P_rank=np.asarray(fold.prank), train_projected_residual_SSE=np.asarray(np.sum((fold.P @ fold.R) ** 2)))
            saved.update({prefix + k: v for k, v in values.items()})
            records.append(dict(task=task["name"], train_ids=task["train_ids"], query_ids=task["query_ids"], training_pair_count=len(task["train_ids"]), permitted_single_count=len(fac.nodes), query_expression_masked=True, original_trace_sha256=sha(ROOT / "run" / task["name"] / "TRACE.json"), theta_max_abs=ce, partner_correction_max_abs=pe, rank=meta["rank"], effective_df=meta["effective_df"]))
    finally:
        model.solve = original_solve
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "SF12_SHARED_TRAIN_STATS.npz", **saved)
    record = dict(status="PASS", seconds=time.monotonic() - start, tasks=records, scope="Exact fixed-settings SF12 outer-support sufficient-statistic derivation from original hash-verified intermediate measurements; no hyperparameter search and no raw-cell preprocessing.", input_hashes={n: sha(p) for n, p in inputs.items()}, code_hashes={str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), ROOT / "source/legacy_model.py", ROOT / "source/fold_factory.py"]})
    (output / "SF12_SUFFICIENT_STATISTICS_DERIVATION.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
