"""Independent SVD loss checks and comparison to archived full-precision outputs."""
import os
for key in ['OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS']:
    os.environ[key] = '1'
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(); ap.add_argument('--expected', type=Path); a = ap.parse_args()
checks = []
scores = pd.read_csv(R / 'results/MATCHED_REFERENCE_SCORES.csv')
for cfg in json.loads((R / 'CONFIGURATIONS.json').read_text()):
    label = cfg['label']; src = R / 'inputs' / label
    tr = np.load(src / 'training.npz'); y = np.load(src / 'outcomes.npz')['y']
    old = np.load(src / 'original.npz'); new = np.load(R / 'fitted' / label / 'PREDICTIONS.npz')
    edges = tr['query_edges']; n = len(tr['singles'])
    X = np.zeros((len(edges), n)); X[np.arange(len(edges))[:, None], edges] = 1
    u, s, _ = np.linalg.svd(X, full_matrices=False)
    u = u[:, s > max(X.shape) * np.finfo(float).eps * s[0]]
    project = lambda v: v - u @ (u.T @ v)
    preds = {'native_B': old['M0'], 'native_rank1': old['rank1_selected'],
             'native_spectral': old['spectral_selected'], 'global_B': new['B']}
    for family in ['rank1', 'spectral']:
        preds['global_' + family] = new['B'] + new[family + '_D']
        preds['global_' + family + '_control'] = new[family + '_controlled']
    maximum = 0.
    for method, pred in preds.items():
        row = scores[(scores.configuration == label) & (scores.method == method) & (scores.repeat == 'pooled')].iloc[0]
        full = np.sum((y - pred) ** 2); partner = np.sum(project(y - pred) ** 2)
        maximum = max(maximum, abs(full - row.full_sse), abs(partner - row.partner_sse))
    assert maximum < 1e-9, (label, maximum)
    endpoint_error = identity_error = 0.
    for family in ['rank1', 'spectral']:
        before = preds['global_' + family]; after = preds['global_' + family + '_control']
        endpoint_error = max(endpoint_error, float(np.max(abs(u.T @ (after - before)))))
        dfull = np.sum((y - before) ** 2) - np.sum((y - after) ** 2)
        dpartner = np.sum(project(y - before) ** 2) - np.sum(project(y - after) ** 2)
        identity_error = max(identity_error, float(abs(dfull - dpartner)))
    assert endpoint_error < 1e-10 and identity_error < 1e-9
    checks.append(dict(configuration=label, independent_sse_max_abs=maximum,
                       endpoint_preservation_max_abs=endpoint_error,
                       control_sse_identity_max_abs=identity_error,
                       training_covers_all_endpoints=bool(len(np.unique(tr['edges'])) == n)))
comparisons = []
if a.expected:
    for p in sorted((R / 'results').glob('*.csv')):
        ref = a.expected / p.name
        if not ref.exists(): continue
        actual = pd.read_csv(p); expected = pd.read_csv(ref)
        pd.testing.assert_frame_equal(actual, expected, check_exact=False, atol=1e-9, rtol=1e-8)
        comparisons.append(p.name)
theory = json.loads((R / 'results/THEORY_CHECKS.json').read_text())
assert theory['status'] == 'PASS'
report = dict(status='PASS', numerical_checks=checks, compared_tables=comparisons,
              scope='Independent SVD scoring of v39 additions; not a new replay of the historical benchmark.')
(R / 'results/VALIDATION.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
