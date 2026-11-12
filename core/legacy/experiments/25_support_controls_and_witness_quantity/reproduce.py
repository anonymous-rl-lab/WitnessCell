#!/usr/bin/env python3
"""Replay the completed follow-up without downloading or reading raw cells.

The original analysis scripts, protocols, chronology and reference outputs are
read-only. This adapter relocates only verified input paths. All reconstruction
and repeated statistics are written to a new output directory.
"""
from __future__ import annotations

import os
for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
STEP1 = "experiments/23_joint_prediction_risk_step1/witnesscell_step1_run"


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                                   allow_nan=False) + "\n", encoding="utf-8")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_inputs(repo):
    result = []
    for entry in read(HERE / "CANONICAL_INPUT_MAP.json")["files"]:
        base = repo if entry["scope"] == "repository" else HERE
        path = base / entry["path"]
        require(path.is_file(), "Missing canonical asset: " + str(path))
        require(sha(path) == entry["sha256"], "Canonical input hash mismatch: " + str(path))
        result.append({"scope": entry["scope"], "path": entry["path"],
                       "sha256": entry["sha256"], "status": "PASS"})
    return result


def array_check(path, entry):
    """Require every original .npy byte, including metadata and numerical bits."""
    expected = {row["member"]: row for row in entry["members"]}
    with zipfile.ZipFile(path) as archive:
        require(set(archive.namelist()) == set(expected), "NPZ member names changed: " + str(path))
        for name, row in expected.items():
            content = archive.read(name)
            require(len(content) == row["size_bytes"], "NPY size changed: " + str(path) + "/" + name)
            require(hashlib.sha256(content).hexdigest() == row["npy_sha256"],
                    "NPY bit identity failed: " + str(path) + "/" + name)
    return {"path": entry["path"], "status": "PASS", "npy_members_bit_exact": len(expected),
            "npz_archive_byte_exact": sha(path) == entry["sha256"]}


def omitted_index():
    return {entry["path"]: entry for entry in read(HERE / "RECONSTRUCTED_ASSETS.json")["files"]}


def witness_statistics(root):
    # Keep the historical lock and all input hashes. Only rerun output files.
    stats = root / "statistics"
    stats.mkdir(exist_ok=True)
    shutil.copyfile(HERE / "wessels_witness/statistics/analysis_freeze.json",
                    stats / "analysis_freeze.json")
    module = load_module("wc25_original_bootstrap", HERE / "wessels_witness/bootstrap_witness_quantity.py")
    module.run(SimpleNamespace(root=root, analysis_plan=HERE / "analysis_plan.json"))
    comparator = load_module("wc25_statistics_comparator", HERE / "reproduction/reproduce_seen.py")
    audit = {"maximum_absolute_numeric_difference": 0.0, "numeric_values_compared": 0}
    comparator.compare_values(read(HERE / "wessels_witness/statistics/summary.json"),
                              read(stats / "summary.json"), "witness_summary", audit)
    audit.update(status="PASS", all_reported_summary_fields_compared=True,
                 bootstrap_replicates=5000, absolute_tolerance=1e-12, relative_tolerance=1e-12)
    audit["bootstrap_arrays"] = array_check(stats / "bootstrap_draws.npz",
                                          omitted_index()["wessels_witness/statistics/bootstrap_draws.npz"])
    return audit


def stage_witness_text(output):
    destination = output / "wessels_witness"
    # This is a fresh reproduction workspace, not a new historical experiment.
    shutil.copytree(HERE / "wessels_witness", destination,
                    ignore=lambda directory, names: [n for n in names if n == "statistics" or n == "__pycache__"])
    return destination


def witness_worker(repo, output, statistics_only):
    start = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    root = stage_witness_text(output)
    report = {"status": "PASS", "scope": "statistics" if statistics_only else "endpoint-to-statistics",
              "new_scientific_experiment": False, "original_inputs_modified": False}
    if not statistics_only:
        report["canonical_input_checks"] = verify_inputs(repo)
        protocol = read(root / "protocol.json")
        config_hash = read(root / "protocol_freeze.json")["protocol_sha256"]
        require(sha(root / "protocol.json") == config_hash, "Protocol hash mismatch")
        require(sha(root / "subsets.json") == protocol["subsets_sha256"], "Subset hash mismatch")
        source = load_module("wc25_original_witness", HERE / "wessels_witness/witness_quantity.py")
        args = SimpleNamespace(repo=repo, step1_root=repo / STEP1, out=root)
        reconstructed = omitted_index()
        checks = []
        log = {"note": "Fresh replay timing/metadata are retained here. Historical metadata are copied back only after hash verification of regenerated numeric files.",
               "replayed": []}
        for seed in (1, 2, 3):
            state, metadata = source.restore_state(args, seed, config_hash)
            relative = f"wessels_witness/seed{seed}/frozen_state.npz"
            path = root / f"seed{seed}/frozen_state.npz"
            checked = array_check(path, reconstructed[relative])
            checks.append(checked)
            # The original unmodified verifier validates the original archive
            # hashes. A different compression implementation may change archive
            # bytes despite identical arrays; fail explicitly rather than
            # rewriting or weakening the original historical hash locks.
            require(checked["npz_archive_byte_exact"],
                    "Arrays are exact but historical archive encoding differs. Use the pinned NumPy/Python environment for the original hash-locked verifier: " + relative)
            log["replayed"].append({"kind": "state", "seed": seed, "metadata": metadata})
            shutil.copyfile(HERE / f"wessels_witness/seed{seed}/restoration.json",
                            root / f"seed{seed}/restoration.json")
            for order, fraction in ((1, 1.0), (1, .25), (1, .5), (2, .25), (2, .5), (3, .25), (3, .5)):
                metadata = source.predict_one(args, seed, order, fraction, state, config_hash)
                name = source.prediction_name(order, fraction)
                relative = f"wessels_witness/seed{seed}/{name}"
                path = root / f"seed{seed}/{name}"
                checked = array_check(path, reconstructed[relative])
                checks.append(checked)
                require(checked["npz_archive_byte_exact"], "Historical NPZ encoding differs: " + relative)
                log["replayed"].append({"kind": "prediction", "seed": seed,
                                         "order": order, "fraction": fraction, "metadata": metadata})
                shutil.copyfile(HERE / f"wessels_witness/seed{seed}/{name.replace('.npz', '.json')}",
                                path.with_suffix(".json"))
            dest = root / f"evaluation_inputs/Wessels_seed{seed}.evaluation_only.npz"
            # These are the same canonical intermediate truth/weight arrays, not
            # a second payload inside the compact release.
            shutil.copyfile(repo / STEP1 / f"results/Wessels_seed{seed}/evaluation_only.npz", dest)
        write(output / "fresh_reconstruction_metadata.json", log)
        verifier = load_module("wc25_original_verifier", HERE / "wessels_witness/replay_packaged.py")
        report["prediction_and_loss_replay"] = verifier.verify(root)
        report["endpoint_training_reconstructed_by_adapter"] = True
        report["endpoint_training_note"] = "The nested original verifier concerns numerical state only; before calling it, this adapter refitted the three original endpoint states using canonical final-v14 code and condition aggregates, and checked all original state array bytes. Canonical code and GO inputs are required for this added replay stage."
        report["reconstructed_array_checks"] = checks
    report["statistics_replay"] = witness_statistics(root)
    report["elapsed_seconds"] = time.monotonic() - start
    write(output / "reproduction_status.json", report)
    return report


def compare_seen_artifacts(output):
    copied = output / "seen"
    checks = []
    # Risk scores, model coefficients and all fixed-coverage accepted identities
    # are regenerated, not merely the aggregate confidence endpoints.
    for dataset in ("Norman", "Wessels"):
        for seed in (1, 2, 3):
            for name in ("new_frozen_scores.csv", "new_frozen_decisions.csv", "new_risk_models.json", "new_scored_queries.csv"):
                rel = f"results/{dataset}_seed{seed}/{name}"
                require(sha(HERE / "seen" / rel) == sha(copied / rel),
                        "Replayed seen artifact differs: " + rel)
                checks.append({"path": rel, "byte_exact": True})
    arrays = array_check(copied / "evaluation/bootstrap_ratios.npz",
                         omitted_index()["seen/evaluation/bootstrap_ratios.npz"])
    return {"status": "PASS", "risk_scores_decisions_models_and_losses": checks,
            "bootstrap_arrays": arrays}


def run(args):
    output = args.output.resolve()
    require(not output.exists(), "Output directory already exists; choose a new path")
    require(output != HERE and HERE not in output.parents, "Output must be outside this frozen supplement")
    if args.mode in ("all", "witness"):
        require(args.repo is not None, "--repo is required to reconstruct witness vectors")
    output.mkdir(parents=True)
    started = time.monotonic()
    result = {"status": "PASS", "mode": args.mode, "new_scientific_experiment": False,
              "source_inputs_read_only": True, "max_seconds": args.max_seconds,
              "steps": {}, "python_version": sys.version}
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    def command(name, argv):
        remaining = args.max_seconds - (time.monotonic() - started)
        require(remaining > 0, "Reproduction time limit exhausted")
        with (output / (name + ".log")).open("w", encoding="utf-8") as stream:
            subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT,
                           env=environment, check=True, timeout=remaining)
    try:
        if args.mode in ("all", "statistics", "seen"):
            destination = output / "seen_reproduction"
            command("seen", [sys.executable, str(HERE / "reproduction/reproduce_seen.py"),
                             "--package-root", str(HERE), "--output", str(destination)])
            result["steps"]["seen_statistics"] = read(destination / "reproduction_status.json")
            result["steps"]["seen_frozen_artifacts"] = compare_seen_artifacts(destination)
        if args.mode in ("all", "statistics", "witness"):
            destination = output / "witness_reproduction"
            argv = [sys.executable, str(Path(__file__).resolve()), "--worker-witness",
                    "--output", str(destination), "--mode", args.mode]
            if args.repo:
                argv += ["--repo", str(args.repo.resolve())]
            command("witness", argv)
            result["steps"]["witness"] = read(destination / "reproduction_status.json")
        result["elapsed_seconds"] = time.monotonic() - started
        write(output / "reproduction_status.json", result)
        return result
    except Exception as error:
        result.update(status="FAIL", error_type=type(error).__name__, error=str(error),
                      elapsed_seconds=time.monotonic() - started)
        write(output / "reproduction_status.json", result)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, help="Extracted canonical repository, grn_repo_v20")
    parser.add_argument("--output", type=Path, required=True, help="New directory outside the frozen supplement")
    parser.add_argument("--mode", choices=("all", "statistics", "seen", "witness"), default="all")
    parser.add_argument("--max-seconds", type=float, default=600,
                        help="Hard total subprocess wall limit; default 600 seconds")
    parser.add_argument("--worker-witness", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_witness:
        report = witness_worker(args.repo.resolve() if args.repo else None,
                                args.output.resolve(), args.mode == "statistics")
    else:
        report = run(args)
    print(json.dumps({"status": report["status"], "elapsed_seconds": report["elapsed_seconds"],
                      "output": str(args.output.resolve())}), flush=True)


if __name__ == "__main__":
    main()
