"""Optional exact reconstruction from hash-verified external intermediate inputs.

The compact replay does not need this command. This runner refits the archived
calibration, M0 and SF12 at fixed settings, retaining all five predefined tasks.
It is not raw-cell preprocessing and performs no hyperparameter search.
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
from pathlib import Path
import argparse
import hashlib
import json
import sys
sys.dont_write_bytecode = True
import time
import numpy as np

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset", type=Path, required=True, help="TRAINING_BATCH_ONLY.npz")
    parser.add_argument("--halves", type=Path, required=True, help="TRAINING_RAW_HALVES.npz")
    parser.add_argument("--knowledge", type=Path, default=ROOT / "data/KNOWLEDGE_NODE.npz")
    parser.add_argument("--output", type=Path, required=True, help="Fresh output directory; never archived run/")
    parser.add_argument("--check-inputs-only", action="store_true")
    args = parser.parse_args()
    args.output = args.output.resolve()
    if args.output == ROOT / "run" or (ROOT / "run") in args.output.parents:
        parser.error("Choose a fresh output directory outside the archived run/.")
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("Output must be absent or empty; archived results will not be overwritten.")
    manifest = json.loads((ROOT / "provenance/RECONSTRUCTION_INPUTS.json").read_text())
    expected = {r["filename"]: r["sha256"] for r in manifest["inputs"]}
    paths = {"TRAINING_BATCH_ONLY.npz": args.asset.resolve(), "TRAINING_RAW_HALVES.npz": args.halves.resolve(), "KNOWLEDGE_NODE.npz": args.knowledge.resolve(), "ORIGINAL_SPLIT.json": ROOT / "source/ORIGINAL_SPLIT.json"}
    checks = []
    for name, path in paths.items():
        if not path.is_file():
            parser.error(f"Missing input: {path}")
        got = sha(path)
        if got != expected[name]:
            parser.error(f"Input checksum mismatch: {name}")
        checks.append(dict(input=name, sha256=got, bytes=path.stat().st_size))
    with np.load(args.asset, allow_pickle=False) as a, np.load(args.knowledge, allow_pickle=False) as k:
        assert np.array_equal(a["nodes"], k["nodes"])
    if args.check_inputs_only:
        print(json.dumps(dict(status="PASS", stage="input_verification_only", inputs=checks, trained=False), indent=2))
        return
    sys.path.insert(0, str(ROOT))
    import run_five as historical
    historical.ASSET = paths["TRAINING_BATCH_ONLY.npz"]
    historical.HALVES = paths["TRAINING_RAW_HALVES.npz"]
    historical.KNOWLEDGE = paths["KNOWLEDGE_NODE.npz"]
    historical.OUT = args.output
    args.output.mkdir(parents=True, exist_ok=True)
    historical.freeze()
    spec = historical.check_frozen()
    start = time.monotonic()
    q4 = spec["tasks"][3]
    pred, _, stats, _ = historical.compute(q4, spec)
    poisoned, _, _, _ = historical.compute(q4, spec, poison_counts=True)
    with np.load(ROOT / "run/COMPACT_VECTORS.npz", allow_pickle=False) as z:
        expected_contrast = z["partner_predictions"][3]
        contrast = z["contrast"]
    with np.load(ROOT / "run/Q4/TRAIN_STATS.npz", allow_pickle=False) as z:
        expected_theta = z["theta"]
    difference = float(np.max(abs(np.einsum("e,meg->mg", contrast, pred) - expected_contrast)))
    theta_difference = float(np.max(abs(stats["theta"] - expected_theta)))
    poison_difference = float(np.max(abs(pred - poisoned)))
    assert difference < 1e-10 and theta_difference < 1e-12 and poison_difference == 0
    assert time.monotonic() - start < 600
    historical.save(args.output / "SMOKE.json", dict(status="PASS", Q4_partner_prediction_max_abs=difference, Q4_theta_max_abs=theta_difference, withheld_counts_poison_max_abs=poison_difference, seconds=time.monotonic() - start, scope="Portable replay of frozen Q4 contrast and coefficient; withheld counts poison test. No new query scoring or selection."))
    historical.predict()  # All five predictions are frozen before score_five.
    import score_five
    score_five.main()
    historical.save(args.output / "PORTABLE_RECONSTRUCTION.json", dict(status="PASS", wrapper_sha256=sha(__file__), inputs=checks, parent_protocol_sha256=sha(ROOT / "run/PROTOCOL.json"), scope="Archived algorithm with explicit relocated intermediate inputs. Exact fixed-parameter refits; no search; not raw-cell preprocessing."))


if __name__ == "__main__":
    main()
