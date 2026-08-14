#!/usr/bin/env python3
"""Build or verify the compact, append-only Experiment 22 release inventory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
E22 = HERE.parent
FORMAL = E22 / "results/formal_e22"
GENERATED = {
    "postrun/ASSET_CATALOG.csv",
    "postrun/RESULT_TABLE_CATALOG.csv",
    "postrun/RELEASE_INTEGRATION_AUDIT.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_manifest(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        digest, relative = line.split(None, 1)
        entries.append((digest, relative.removeprefix("./")))
    return entries


def verify_manifest(path: Path, root: Path, external_prefix: str | None = None) -> dict:
    entries = parse_manifest(path)
    verified = 0
    external = []
    errors = []
    for expected, relative in entries:
        target = root / relative
        if not target.exists() and external_prefix and relative.startswith(external_prefix):
            external.append(relative)
            continue
        if not target.is_file():
            errors.append({"path": relative, "error": "missing"})
            continue
        actual = sha256(target)
        if actual != expected:
            errors.append({"path": relative, "expected": expected, "actual": actual})
            continue
        verified += 1
    return {
        "entries": len(entries),
        "verified": verified,
        "external_not_duplicated": external,
        "errors": errors,
        "pass": not errors and verified + len(external) == len(entries),
    }


def classify(relative: str) -> str:
    if relative.startswith("results/formal_e22/"):
        return "formal_result"
    if relative.startswith("assets/formal_contracts/") or relative.startswith("assets/gate21_gene_contract"):
        return "derived_frozen_contract"
    if relative.startswith("amendments/"):
        return "execution_amendment"
    if relative.startswith("audit/"):
        return "prefreeze_audit"
    if relative.startswith(("src/", "tests/", "firewall/")) or relative.startswith("run_"):
        return "executable_code"
    if relative.startswith("assets/smoke/") or relative.startswith("results/smoke/"):
        return "development_smoke"
    if relative in {
        "PROTOCOL_v2.md",
        "PRE_FREEZE_CHECKLIST.md",
        "FROZEN_MANIFEST.sha256",
        "FREEZE_RECEIPT.json",
        "PREFREEZE_CORRECTIONS.md",
        "SUPERSESSION_NOTICE.md",
    }:
        return "frozen_protocol"
    if relative in {
        "SOURCE_LOCK.json",
        "COMPARATOR_SOURCE_LOCK.json",
        "COMPARATOR_SCOPE.json",
        "comparator_inventory.csv",
        "THIRD_PARTY_NOTICE.md",
    }:
        return "source_and_comparator_lock"
    if relative in {"FORMAL_REPORT_CN.md", "POSTRUN_RECEIPT.json", "README.md", "VERSION"}:
        return "report_and_identity"
    if relative.startswith("postrun/"):
        return "release_audit"
    return "supporting_asset"


def build_asset_rows() -> list[dict]:
    rows = []
    for path in sorted(E22.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(E22).as_posix()
        if relative in GENERATED:
            continue
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        if relative.startswith((".venv/", ".cache/", ".pytest_cache/", "assets/cell_level/")):
            continue
        rows.append(
            {
                "relative_path": relative,
                "category": classify(relative),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return rows


def csv_shape(path: Path) -> tuple[int, int]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = sum(1 for _ in reader)
    return rows, len(header)


def build_table_rows() -> list[dict]:
    rows = []
    for path in sorted(FORMAL.rglob("*.csv")):
        n_rows, n_columns = csv_shape(path)
        rows.append(
            {
                "relative_path": path.relative_to(E22).as_posix(),
                "phase": path.relative_to(FORMAL).parts[0],
                "rows": n_rows,
                "columns": n_columns,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write canonical catalogs before release freeze")
    args = parser.parse_args()

    required = [
        "PROTOCOL_v2.md",
        "PRE_FREEZE_CHECKLIST.md",
        "FROZEN_MANIFEST.sha256",
        "FREEZE_RECEIPT.json",
        "SOURCE_LOCK.json",
        "src/metric_core.py",
        "src/phase_m_metric_validity.py",
        "src/phase_p_prediction_stress.py",
        "src/phase_d_decision_stress.py",
        "assets/CONTRACT_SEAL.json",
        "assets/gate21_gene_contract.npz",
        "results/formal_e22/FORMAL_VERDICT.json",
        "results/formal_e22/FORMAL_RESULTS_MANIFEST.sha256",
        "FORMAL_REPORT_CN.md",
        "POSTRUN_RECEIPT.json",
        "amendments/A1_phase_d_index_reset/AMENDMENT_A1_MANIFEST.sha256",
        "amendments/A2_phase_d_complete_case_alignment/AMENDMENT_A2_MANIFEST.sha256",
    ]
    missing_required = [relative for relative in required if not (E22 / relative).is_file()]

    # The scientific freeze seals the cell-asset audit ledger and all derived
    # contracts, not the multi-gigabyte public matrices themselves.
    frozen = verify_manifest(E22 / "FROZEN_MANIFEST.sha256", E22)
    formal = verify_manifest(
        FORMAL / "FORMAL_RESULTS_MANIFEST.sha256", FORMAL
    )
    amendment_a1 = verify_manifest(
        E22 / "amendments/A1_phase_d_index_reset/AMENDMENT_A1_MANIFEST.sha256",
        E22 / "amendments/A1_phase_d_index_reset",
    )
    amendment_a2 = verify_manifest(
        E22 / "amendments/A2_phase_d_complete_case_alignment/AMENDMENT_A2_MANIFEST.sha256",
        E22 / "amendments/A2_phase_d_complete_case_alignment",
    )

    cell_ledger = json.loads((E22 / "audit/cell_asset_ledger.json").read_text())
    external_expected = {
        Path(item[key]).relative_to("experiments/22_metric_calibration_stress").as_posix()
        for item in cell_ledger["datasets"]
        for key in ("compressed", "h5ad")
    }
    external_present = {
        relative for relative in external_expected if (E22 / relative).exists()
    }
    external_ledger_pass = (
        cell_ledger.get("pass") is True
        and len(cell_ledger["datasets"]) == 4
        and len(external_expected) == 8
        and not external_present
    )

    verdict = json.loads((FORMAL / "FORMAL_VERDICT.json").read_text())
    postrun = json.loads((E22 / "POSTRUN_RECEIPT.json").read_text())
    asset_rows = build_asset_rows()
    table_rows = build_table_rows()

    report = {
        "status": "PASS_GATE22_COMPACT_RELEASE_INTEGRATION",
        "experiment": "Experiment 22 metric-calibration stress test",
        "required_assets_missing": missing_required,
        "frozen_manifest": frozen,
        "formal_results_manifest": formal,
        "amendment_a1_manifest": amendment_a1,
        "amendment_a2_manifest": amendment_a2,
        "public_cell_assets": {
            "policy": "external_not_duplicated",
            "datasets": len(cell_ledger["datasets"]),
            "files": len(external_expected),
            "files_present_in_compact_release": len(external_present),
            "ledger_pass": external_ledger_pass,
            "source": cell_ledger["source"],
        },
        "release_assets": {
            "catalogued_files": len(asset_rows),
            "formal_csv_tables": len(table_rows),
            "formal_csv_rows": sum(int(row["rows"]) for row in table_rows),
        },
        "formal_verdict": {
            "status": verdict["status"],
            "metric_validity": verdict["METRIC_VALIDITY"],
            "prediction_uninformative": verdict["PRED_UNINFORMATIVE"],
            "prediction_linear": verdict["PRED_LINEAR"],
            "gate21_full": verdict["GATE21_WMSE"],
            "shadow_gate": verdict["SHADOW_GATE"],
        },
        "hash_chain": {
            "frozen_manifest_sha256": sha256(E22 / "FROZEN_MANIFEST.sha256"),
            "formal_results_manifest_sha256": sha256(FORMAL / "FORMAL_RESULTS_MANIFEST.sha256"),
            "postrun_receipt_parent_frozen_manifest_sha256": postrun[
                "parent_frozen_manifest_sha256"
            ],
            "postrun_receipt_formal_results_manifest_sha256": postrun[
                "formal_results_manifest_sha256"
            ],
        },
    }
    checks = [
        not missing_required,
        frozen["pass"],
        formal["pass"],
        amendment_a1["pass"],
        amendment_a2["pass"],
        external_ledger_pass,
        verdict["status"] == "COMPLETE_EXPERIMENT22_FORMAL_EXECUTION",
        postrun["status"] == "COMPLETE_EXPERIMENT22_FORMAL_EXECUTION",
        report["hash_chain"]["frozen_manifest_sha256"]
        == postrun["parent_frozen_manifest_sha256"],
        report["hash_chain"]["formal_results_manifest_sha256"]
        == postrun["formal_results_manifest_sha256"],
    ]
    if not all(checks):
        report["status"] = "FAIL_GATE22_COMPACT_RELEASE_INTEGRATION"

    if args.write:
        HERE.mkdir(parents=True, exist_ok=True)
        write_csv(HERE / "ASSET_CATALOG.csv", asset_rows)
        write_csv(HERE / "RESULT_TABLE_CATALOG.csv", table_rows)
        (HERE / "RELEASE_INTEGRATION_AUDIT.json").write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n"
        )
    else:
        stored = json.loads((HERE / "RELEASE_INTEGRATION_AUDIT.json").read_text())
        if stored != report:
            raise SystemExit("release integration audit differs from sealed catalog")
        with (HERE / "ASSET_CATALOG.csv").open(newline="", encoding="utf-8") as handle:
            stored_assets = list(csv.DictReader(handle))
        with (HERE / "RESULT_TABLE_CATALOG.csv").open(newline="", encoding="utf-8") as handle:
            stored_tables = list(csv.DictReader(handle))
        normalized_assets = [{key: str(value) for key, value in row.items()} for row in asset_rows]
        normalized_tables = [{key: str(value) for key, value in row.items()} for row in table_rows]
        if stored_assets != normalized_assets:
            raise SystemExit("asset catalog differs from current release")
        if stored_tables != normalized_tables:
            raise SystemExit("result-table catalog differs from current release")

    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
