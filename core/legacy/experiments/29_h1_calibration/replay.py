#!/usr/bin/env python3
"""Replay the H1 calibration experiment from compact sufficient statistics.

Recomputes errors and selection, and independently solves constrained WLS
preprocessing operators. Does not reconstruct full-gene predictions or raw cells.
Original reference files are never overwritten.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def independent_operator(design, record):
    """Solve a condition-intercept + constrained batch-offset WLS model."""
    qi = design['observation_condition']; bi = design['observation_batch']
    counts = design['counts']; pi = design['pi']
    allowed = list(range(10)) + record['train_global']
    batches = list(range(9)) if record['arm'] != 'B2_cal7' else list(np.flatnonzero(pi))
    L = np.zeros((len(allowed), len(qi)))
    if record['arm'] in ('raw', 'B1_direct7'):
        for row, qid in enumerate(allowed):
            props = pi if record['arm'] == 'B1_direct7' else counts[qid] / counts[qid].sum()
            for b, proportion in enumerate(props):
                observed = np.flatnonzero((qi == qid) & (bi == b))
                if len(observed):
                    L[row, observed[0]] = proportion
        return L

    localpi = pi[batches]
    pivot = int(np.flatnonzero(localpi)[0])
    other = [i for i in range(len(batches)) if i != pivot]
    Z = np.zeros((len(batches), len(other)))
    for col, j in enumerate(other):
        Z[j, col] = 1
        Z[pivot, col] = -localpi[j] / localpi[pivot]
    obs = np.flatnonzero(np.isin(qi, allowed) & np.isin(bi, batches))
    D = np.zeros((len(obs), len(allowed) + len(other)))
    weights = []
    for row, o in enumerate(obs):
        a = allowed.index(int(qi[o])); b = batches.index(int(bi[o]))
        D[row, a] = 1
        D[row, len(allowed):] = Z[b]
        weights.append(counts[qi[o], bi[o]] / counts[qi[o], batches].sum())
    sw = np.sqrt(weights)
    op = np.linalg.pinv(D * sw[:, None]) * sw[None, :]
    A = np.zeros((len(batches), len(qi)))
    A[:, obs] = Z @ op[len(allowed):]
    for row, qid in enumerate(allowed):
        props = counts[qid, batches] / counts[qid, batches].sum()
        for b, proportion in zip(batches, props):
            observed = np.flatnonzero((qi == qid) & (bi == b))
            if len(observed):
                L[row, observed[0]] = proportion
        L[row] -= (props - localpi) @ A
    return L


def replay():
    ref = HERE / 'reference'
    provenance = json.loads((HERE / 'SOURCE_PROVENANCE.json').read_text())
    for row in provenance['original_files']:
        got = hashlib.sha256((HERE / row['path']).read_bytes()).hexdigest()
        check(got == row['sha256'], 'Original source differs: ' + row['path'])
    protocol = json.loads((ref / 'PROTOCOL_BEFORE_FIT.json').read_text())
    expected = json.loads((ref / 'SUMMARY.json').read_text())
    rows = json.loads((ref / 'SCORE_GRAM_ROWS.json').read_text())
    K = np.load(ref / 'SCORE_GRAM.npz')['gram']
    truth = {(r['global_id'], r['half']): i for i, r in enumerate(rows) if r['type'] == 'reference'}
    with (ref / 'PER_CONDITION_HALF_ERRORS.csv').open() as f:
        archived = {(int(r['fold']), r['arm'], r['condition'], int(r['half'])): float(r['mse']) for r in csv.DictReader(f)}
    scored = {}
    for i, row in enumerate(rows):
        if row['type'] == 'prediction':
            for half in (0, 1):
                j = truth[row['global_id'], half]
                scored[row['fold'], row['arm'], row['condition'], half] = float(K[i, i] + K[j, j] - 2 * K[i, j])
    check(scored.keys() == archived.keys(), 'Scoring row identities differ')
    score_diff = max(abs(v - archived[k]) for k, v in scored.items())

    selection = json.loads((ref / 'SELECTION.json').read_text())
    stats = np.load(ref / 'SELECTION_QUADRATICS.npz')
    gamma_diff = risk_diff = 0.
    for s in selection:
        prefix = f"fold{s['fold']}__{s['arm']}"
        t, dot, q = [stats[prefix + suffix] for suffix in ('__truth2', '__dot', '__raw2')]
        candidates = []
        for ai, _alpha in enumerate(protocol['baseline']['alphas']):
            for ni, eta in enumerate(protocol['baseline']['noises']):
                gamma = float(np.clip(dot[ai, ni] / q[ai, ni], 0, 1)) if q[ai, ni] > 1e-12 else 0.
                loss = float(t[ai] - 2 * gamma * dot[ai, ni] + gamma ** 2 * q[ai, ni])
                candidates.append((loss, -eta, ai, ni, gamma))
        best = min(candidates)
        check([best[2], best[3]] == s['index'], 'Selected candidate differs')
        gamma_diff = max(gamma_diff, abs(best[4] - s['selected']['gamma']))
        denominator = len(s['outertrain']) * expected['all_gene_count']
        risk_diff = max(risk_diff, abs(best[0] / denominator - s['oof_mse']))

    design = np.load(HERE / 'metadata/design.npz')
    fits = json.loads((ref / 'PREPROCESS_FITS.json').read_text())
    saved = np.load(ref / 'PREPROCESS_COEFFICIENTS.npz')
    names = {k: i for i, k in enumerate(design['condition_ids'].astype(str))}
    wls_diff = excluded_max = affine_error = 0.
    for fit in fits:
        L = saved[fit['key']]
        reconstructed = independent_operator(design, fit)
        wls_diff = max(wls_diff, float(np.max(abs(reconstructed - L))))
        excluded = [names[n] for n in fit['excluded_pairs']]
        excluded_max = max(excluded_max, float(np.max(abs(L[:, np.isin(design['observation_condition'], excluded)]))))
        affine_error = max(affine_error, float(np.max(abs(L.sum(1) - 1))))

    conditions = sorted({k[2] for k in scored})
    arms = protocol['arms']
    vals = {a: np.array([np.mean([v for k, v in scored.items() if k[1] == a and k[2] == c]) for c in conditions]) for a in arms}
    mse = {a: float(v.mean()) for a, v in vals.items()}
    gains = {a: 100 * (1 - v / mse['raw']) for a, v in mse.items()}
    loo = [100 * (1 - np.delete(vals['B2_cal9'], i).sum() / np.delete(vals['raw'], i).sum()) for i in range(len(conditions))]
    fold_gains = {}
    for fold in sorted({k[0] for k in scored}):
        fold_mse = {a: np.mean([v for k, v in scored.items() if k[0] == fold and k[1] == a]) for a in arms}
        fold_gains[str(fold)] = 100 * (1 - fold_mse['B2_cal9'] / fold_mse['raw'])
    half_gains = {str(h): 100 * (1 - np.mean([v for k, v in scored.items() if k[1] == 'B2_cal9' and k[3] == h]) / np.mean([v for k, v in scored.items() if k[1] == 'raw' and k[3] == h])) for h in (0, 1)}
    headline_diff = max(abs(mse[a] - expected['arms_mse'][a]) for a in arms)
    check(len(scored) == 168 and len(selection) == 20 and len(fits) == 80, 'Unexpected experiment coverage')
    check(score_diff < 1e-12 and gamma_diff < 1e-12 and risk_diff < 1e-12, 'Numerical replay failed')
    check(wls_diff < 1e-10 and affine_error < 1e-10 and excluded_max == 0, 'Preprocessing reconstruction/isolation failed')
    check(headline_diff < 1e-12, 'Headline values differ')
    return {
        'status': 'PASS',
        'scope': 'Compact full-RNA score and selection replay; independent preprocessing-operator fit. No raw-cell or full-expression retraining.',
        'input_files_hash_verified': len(provenance['original_files']),
        'condition_half_scores': len(scored), 'selected_models_verified': len(selection),
        'preprocessing_operators_refit': len(fits),
        'max_score_abs': score_diff, 'max_selection_gamma_abs': gamma_diff,
        'max_selection_risk_abs': risk_diff, 'independent_WLS_operator_max_abs': wls_diff,
        'held_condition_coefficient_max': excluded_max, 'affine_row_sum_error': affine_error,
        'max_headline_mse_abs': headline_diff,
        'MSE': mse, 'gain_vs_raw_pct': gains,
        'cal7_gain_vs_direct7_pct': 100 * (1 - mse['B2_cal7'] / mse['B1_direct7']),
        'cal9_gain_vs_direct7_pct': 100 * (1 - mse['B2_cal9'] / mse['B1_direct7']),
        'cal9_improved_conditions': int(np.sum(vals['B2_cal9'] < vals['raw'])),
        'condition_count': len(conditions), 'cal9_vs_raw_fold_gains_pct': fold_gains,
        'cal9_vs_raw_technical_half_gains_pct': half_gains,
        'leave_one_condition_cal9_vs_raw_range_pct': [min(loo), max(loo)]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional JSON report path. No reference data are modified.')
    args = parser.parse_args()
    result = replay()
    text = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
