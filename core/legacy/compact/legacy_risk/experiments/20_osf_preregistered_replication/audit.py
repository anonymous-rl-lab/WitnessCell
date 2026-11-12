#!/usr/bin/env python3
"""Fail-closed audit for the OSF xrsmd computational-replication assets."""
from __future__ import annotations

from io import BytesIO
from pathlib import Path
import csv
import hashlib
import json
import tarfile
import zipfile

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
EXPECTED_SOURCE_SHA = "955d446a0cfe65a6f4f191d74c7a7d16a6083488b9beb40a5e155ab94d528bd7"
EXPECTED_PREREG_SHA = "8d2bcc1c06a42fee039a515ee12588d516eb04c9979697181eb35ebd4c4a8c63"
EXPECTED_RESULT_SHA = "511860f645158a9af536a3944e03d4cc440674e23b55d826bf615e02d96bea0c"
EXPECTED_FALLBACKS = {
    ("Schmidt", 3),
    ("Wessels", 1),
    ("Wessels", 2),
    ("Wessels", 3),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tar_bytes(archive: tarfile.TarFile, member: str) -> bytes:
    stream = archive.extractfile(member)
    assert stream is not None, member
    return stream.read()


def check_protocol_manifest(protocol_root: Path) -> None:
    manifest = protocol_root / "FREEZE_MANIFEST.sha256"
    for raw_line in manifest.read_text().splitlines():
        if not raw_line.strip():
            continue
        expected, relative = raw_line.split(maxsplit=1)
        path = protocol_root / relative.removeprefix("*").removeprefix("./")
        assert path.is_file(), path
        assert sha256(path) == expected, path


def main() -> None:
    prereg_zip = ROOT / "protocol/WitnessCell_OSF_Computational_Preregistration_v2.zip"
    protocol_root = ROOT / "protocol/WitnessCell_OSF_Computational_Preregistration_v2"
    result_archive = ROOT / "result_archive/WitnessCell_preregistered_replication_20260810T105615Z.tar.gz"

    assert sha256(prereg_zip) == EXPECTED_PREREG_SHA
    assert sha256(result_archive) == EXPECTED_RESULT_SHA
    with zipfile.ZipFile(prereg_zip) as archive:
        assert archive.testzip() is None
    check_protocol_manifest(protocol_root)

    expected_verdict = json.loads((ROOT / "evidence/CONFIRMATORY_VERDICT.json").read_text())
    expected_protocol = json.loads((ROOT / "evidence/FROZEN_PROTOCOL.json").read_text())
    assert expected_protocol["protocol_version"] == "2.0.0"
    assert expected_protocol["source"]["sha256"] == EXPECTED_SOURCE_SHA
    assert expected_protocol["wet_lab_component"] is False
    assert expected_protocol["new_biological_data_collection"] is False

    active: set[tuple[str, int]] = set()
    inactive: set[tuple[str, int]] = set()
    comparisons = 0
    global_max_delta = 0.0

    with tarfile.open(result_archive, "r:gz") as archive:
        names = archive.getnames()
        prefixes = {name.split("/", 1)[0] for name in names if "/" in name}
        assert len(prefixes) == 1, prefixes
        prefix = next(iter(prefixes))

        embedded_verdict = json.loads(tar_bytes(archive, f"{prefix}/CONFIRMATORY_VERDICT.json"))
        embedded_protocol = json.loads(tar_bytes(archive, f"{prefix}/FROZEN_PROTOCOL.json"))
        assert embedded_verdict == expected_verdict
        assert embedded_protocol == expected_protocol
        assert tar_bytes(archive, f"{prefix}/environment.txt") == (ROOT / "evidence/environment.txt").read_bytes()
        assert tar_bytes(archive, f"{prefix}/run_started_utc.txt") == (ROOT / "evidence/run_started_utc.txt").read_bytes()
        assert tar_bytes(archive, f"{prefix}/run_completed_utc.txt") == (ROOT / "evidence/run_completed_utc.txt").read_bytes()

        deploy_members = sorted(
            name for name in names
            if name.endswith("/deploy_predictions.npz") and "/gate19_regenerated/" in name
        )
        assert len(deploy_members) == 12, len(deploy_members)

        for member in deploy_members:
            parts = member.split("/")
            gate_index = parts.index("gate19_regenerated")
            dataset = parts[gate_index + 1]
            seed_name = parts[gate_index + 2]
            seed = int(seed_name.removeprefix("seed"))
            frozen_path = (
                REPO
                / "experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo"
                / dataset
                / seed_name
                / "deploy_predictions.npz"
            )
            assert frozen_path.is_file(), frozen_path
            with np.load(BytesIO(tar_bytes(archive, member)), allow_pickle=False) as generated, np.load(
                frozen_path, allow_pickle=False
            ) as frozen:
                assert len(generated.files) == 13
                assert set(generated.files) == set(frozen.files)
                assert "truth" not in generated.files and "truth_variance" not in generated.files
                for key in generated.files:
                    left, right = generated[key], frozen[key]
                    assert left.shape == right.shape
                    assert left.dtype == right.dtype
                    assert left.dtype.kind != "O"
                    assert np.array_equal(left, right), (dataset, seed, key)
                    if left.dtype.kind not in "US" and left.size:
                        assert np.isfinite(left).all(), (dataset, seed, key)
                        global_max_delta = max(
                            global_max_delta,
                            float(np.max(np.abs(left - right))),
                        )
                    comparisons += 1

            manifest_member = member.rsplit("/", 1)[0] + "/manifest.json"
            manifest = json.loads(tar_bytes(archive, manifest_member))
            key = (dataset, seed)
            if bool(manifest["identity_head"]["correction_active"]):
                active.add(key)
            else:
                inactive.add(key)

    assert comparisons == 156
    assert global_max_delta == 0.0
    assert len(active) == 8
    assert inactive == EXPECTED_FALLBACKS

    with (ROOT / "evidence/paired_summary.csv").open(newline="") as stream:
        summary = next(csv.DictReader(stream))
    assert summary["scope"] == "all"
    assert int(summary["condition_seed_units"]) == 654
    assert float(summary["relative_all_gene_mse_gain"]) == 0.001079217249737252
    assert float(summary["relative_top100_mse_gain"]) == 0.0030453806218846413
    assert float(summary["top100_pcc_delta"]) == 0.0013006984799760163

    assert expected_verdict["status"] == "PASS_PREREGISTERED_COMPUTATIONAL_REPLICATION"
    for criterion in (
        "C1_execution_schema",
        "C2_numerical_identity",
        "C3_gate_fallback_identity",
        "C4_paired_mean_response",
        "C5_inherited_gate18_aggregation_integrity",
    ):
        assert expected_verdict[criterion] == "PASS"
    assert expected_verdict["new_biological_measurements"] is False
    assert expected_verdict["raw_h5ad_six_metric_rescoring_performed"] is False

    result = {
        "status": "PASS_GATE20_OSF_PREREGISTERED_REPLICATION_ASSET_AUDIT",
        "registration_id": "xrsmd",
        "protocol_version": "2.0.0",
        "source_archive_sha256": EXPECTED_SOURCE_SHA,
        "preregistration_package_sha256": EXPECTED_PREREG_SHA,
        "result_archive_sha256": EXPECTED_RESULT_SHA,
        "prediction_archives": 12,
        "arrays_compared": comparisons,
        "global_max_abs_delta": global_max_delta,
        "correction_active_splits": len(active),
        "exact_fallback_splits": len(inactive),
        "paired_units": 654,
        "criteria": {name: "PASS" for name in ("C1", "C2", "C3", "C4", "C5")},
        "new_biological_measurements": False,
        "raw_h5ad_gate19_v14_rescoring": False,
    }
    assert result == json.loads((ROOT / "RELEASE_AUDIT.json").read_text())
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

