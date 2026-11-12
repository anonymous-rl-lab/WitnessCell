#!/usr/bin/env python3
"""Load the archived renderer and supply missing corner labels in memory only.

The archived numerical CSV contains four ordered pairs but omits row_endpoint
and column_endpoint expected by the original plotting script. This adapter
derives those labels by set intersection. It never changes CSV values, fits,
predictions, the embedded archive, or the archived renderer source.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import runpy
from pathlib import Path


def infer_corner_labels(rows):
    case_reports = []
    for case in sorted({r["case"] for r in rows}):
        group = [r for r in rows if r["case"] == case]
        if len(group) != 4:
            raise ValueError("Expected exactly four corners: " + case)
        for column in ("gene", "rectangle", "fold"):
            if len({r[column] for r in group}) != 1:
                raise ValueError("Corners must share readout/rectangle/fit: " + case)
        lookup = {}
        for row in group:
            raw = row["corner"]
            corner = {"1": "AB", "2": "AD", "3": "CB", "4": "CD"}.get(raw, raw)
            if corner not in ("AB", "AD", "CB", "CD") or corner in lookup:
                raise ValueError("Unknown/duplicate corner: " + repr(raw))
            pair = row["pair"].split("+")
            if len(pair) != 2 or len(set(pair)) != 2:
                raise ValueError("Expected an unordered pair of distinct genes")
            expected_sign = 1 if corner in ("AB", "CD") else -1
            if float(row["contrast_sign"]) != expected_sign:
                raise ValueError("Four-corner signs disagree: " + case)
            lookup[corner] = (row, set(pair))
        intersections = {
            "A": lookup["AB"][1] & lookup["AD"][1],
            "B": lookup["AB"][1] & lookup["CB"][1],
            "C": lookup["CB"][1] & lookup["CD"][1],
            "D": lookup["AD"][1] & lookup["CD"][1],
        }
        if any(len(v) != 1 for v in intersections.values()):
            raise ValueError("Pair identities do not determine four-corner axes")
        labels = {key: next(iter(value)) for key, value in intersections.items()}
        if len(set(labels.values())) != 4:
            raise ValueError("Expected four distinct rectangle endpoints")
        for corner, (row, pair) in lookup.items():
            expected = {labels[corner[0]], labels[corner[1]]}
            if pair != expected:
                raise ValueError("Inconsistent rectangle edge identities")
            for column, label in (("row_endpoint", labels[corner[0]]),
                                  ("column_endpoint", labels[corner[1]])):
                if column in row and row[column] != label:
                    raise ValueError("Existing corner metadata conflicts")
                row[column] = label
        case_reports.append({"case": case, "fold": group[0]["fold"],
                             "rectangle": group[0]["rectangle"], "axes": labels})
    return case_reports


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", required=True, type=Path)
    args = parser.parse_args()
    payload = args.payload.resolve()
    script = payload / "code/build_display_figures.py"
    table = payload / "rna/selected_four_corners.csv"
    before = {"renderer": sha(script), "corner_csv": sha(table)}
    namespace = runpy.run_path(str(script), run_name="witnesscell_archived_renderer")
    rows = namespace["corners"]
    original = [dict(row) for row in rows]
    cases = infer_corner_labels(rows)
    # Original numerical and identity fields must remain character-for-character.
    if any(any(row[key] != value for key, value in old.items())
           for old, row in zip(original, rows)):
        raise ValueError("Compatibility adapter altered an archived data field")
    for function in ("rna_main", "rna_context", "jurkat_main", "jurkat_context"):
        namespace[function]()
    after = {"renderer": sha(script), "corner_csv": sha(table)}
    if before != after:
        raise ValueError("Renderer or numerical source CSV changed on disk")
    print(json.dumps({
        "passed": True,
        "adaptation": "Derive plotting row/column labels from four ordered pair identities in memory",
        "case_axes": cases,
        "archived_renderer_and_corner_csv_unchanged": True,
        "source_sha256": after,
        "new_model_fits": 0,
        "figures_created": 4,
        "formats": ["pdf", "svg", "png"],
    }, indent=2))


if __name__ == "__main__":
    main()
