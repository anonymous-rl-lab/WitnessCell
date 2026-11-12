#!/usr/bin/env python3
"""Recompute RNA M0/M1/uniform/R results from the inherited signed E28 rows.

No duplicated data, downloads, model fits or selection based on test labels.
This rescoring is not a replay of training or inner selection.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np

ARCHIVE = 'experiments/28_rna_jurkat_response_examples/WitnessCell_RNA_Jurkat_Response_Examples_20261019.zip'
ARCHIVE_SHA256 = '7f9f63bc812df611955a824b2573066d893c5c6d63b8bfc6294e550db562ec41'
MEMBER = 'WitnessCell_RNA_Jurkat_Response_Examples_20261019/rna/all_rectangle_readout_candidates.csv.gz'
EXPECTED = {'M1': 9.438239831528272, 'scalar': 9.698805654041065, 'R': 15.96223210421443}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--out', type=Path, required=True, help='New external result JSON path')
    args = parser.parse_args()
    repo, out = args.repo.resolve(), args.out.resolve()
    if (repo / 'VERSION').read_text().strip() != '25.0.0':
        raise ValueError('Expected grn_repo_v25 VERSION 25.0.0')
    if out.exists() or out.is_relative_to(repo):
        raise ValueError('--out must be a new file outside the release')
    archive = repo / ARCHIVE
    if archive.is_symlink() or hashlib.sha256(archive.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise ValueError('E28 archive identity mismatch')
    with zipfile.ZipFile(archive) as z:
        rows = list(csv.DictReader(io.StringIO(gzip.decompress(z.read(MEMBER)).decode('utf-8'))))
    identities = [(int(r['rectangle']), r['gene']) for r in rows]
    if len(rows) != 11 * 4995 or len(set(identities)) != len(rows):
        raise ValueError('Expected exactly 11 x 4,995 unique rectangle/readout rows')
    rectangles = sorted({q[0] for q in identities})
    genes = sorted({q[1] for q in identities})
    if rectangles != list(range(1, 12)) or len(genes) != 4995:
        raise ValueError('Rectangle or readout identities differ')
    values = {(int(r['rectangle']), r['gene']): r for r in rows}
    observed = np.array([[float(values[q, g]['C_observed']) for g in genes] for q in rectangles])
    predictions = {m: np.array([[float(values[q, g]['C_' + m]) for g in genes] for q in rectangles])
                   for m in ['M0', 'M1', 'scalar', 'R']}
    errors = {m: (v - observed) ** 2 for m, v in predictions.items()}
    mse = {m: float(v.mean()) for m, v in errors.items()}
    gain = {m: float(100 * (1 - mse[m] / mse['M0'])) for m in EXPECTED}
    for model, expected in EXPECTED.items():
        if abs(gain[model] - expected) > 1e-6:
            raise AssertionError((model, gain[model], expected))
    r_vs_uniform = float(100 * (1 - mse['R'] / mse['scalar']))
    if abs(r_vs_uniform - 6.9361502) > 1e-6:
        raise AssertionError(('R_vs_uniform', r_vs_uniform))
    per_rectangle = []
    for i, rectangle in enumerate(rectangles):
        er = {m: float(e[i].sum()) for m, e in errors.items()}
        per_rectangle.append({'rectangle': rectangle, 'SSE': er,
                              'gain_R_vs_uniform_pct': 100 * (1 - er['R'] / er['scalar'])})
    if sum(r['gain_R_vs_uniform_pct'] > 0 for r in per_rectangle) != 11:
        raise AssertionError('R must improve uniform in all eleven rectangles')
    result = {'status': 'PASS', 'repository_version': '25.0.0', 'source_archive': ARCHIVE,
              'source_sha256': ARCHIVE_SHA256, 'signed_rows': len(rows), 'rectangles': 11,
              'non_target_readouts': 4995, 'MSE': mse, 'gain_vs_M0_pct': gain,
              'R_vs_uniform_gain_pct': r_vs_uniform, 'per_rectangle': per_rectangle,
              'scope': 'Exact scoring from saved signed predictions; no training or inner-selection replay',
              'interpretation': 'Eleven overlapping contrasts are dependent; optimization repeats are not biological replicates.',
              'model_identity': 'Matched RNA M0/M1/R; distinct from H1 SF12 and case_R.',
              'new_model_fits': 0, 'external_downloads': 0}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['status', 'signed_rows', 'MSE', 'gain_vs_M0_pct', 'R_vs_uniform_gain_pct']}, indent=2))

if __name__ == '__main__':
    main()
