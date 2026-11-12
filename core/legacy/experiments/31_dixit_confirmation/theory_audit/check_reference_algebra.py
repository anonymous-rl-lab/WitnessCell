#!/usr/bin/env python3
"""Independent algebra checks for fixed-reference partner delivery.

Uses a deterministic synthetic kernel on exactly the ten-endpoint / 45-pair
complete graph. No measured outcomes, model selection, or result optimization.
Writes quantitative checks suitable for scientific audit/reproduction.
"""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


def incidence(edges, n):
    x = np.zeros((len(edges), n))
    for i, (a, b) in enumerate(edges):
        x[i, a] = x[i, b] = 1.0
    return x


def residualizer(x):
    return np.eye(len(x)) - x @ np.linalg.pinv(x)


def weighted_left_inverse(x, w):
    return np.linalg.pinv(x * np.sqrt(w)[:, None]) * np.sqrt(w)[None, :]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2026102019)
    n = 10
    edges = list(itertools.combinations(range(n), 2))
    x = incidence(edges, n)
    n_pairs = len(edges)
    endpoint = rng.normal(size=(n, 13))
    eg = endpoint @ endpoint.T / endpoint.shape[1]
    a, b = np.asarray(edges).T
    k = 0.5 * (eg[np.ix_(a, a)] * eg[np.ix_(b, b)]
               + eg[np.ix_(a, b)] * eg[np.ix_(b, a)])
    # One frozen training fit, many differently assembled request batches.
    tr = np.sort(rng.choice(n_pairs, 30, replace=False))
    q = np.asarray(sorted(set(range(n_pairs)) - set(tr)))
    xt = x[tr]
    pt = residualizer(xt)
    ktt = k[np.ix_(tr, tr)]
    projected = pt @ ktt @ pt
    scale = np.trace(projected) / len(tr)
    y = rng.normal(size=(len(tr), 3))
    c = np.linalg.solve(projected / scale + 0.7 * np.eye(len(tr)), y)
    g = (k[:, tr] - x @ np.linalg.pinv(xt) @ ktt) @ pt @ c / scale
    left = np.linalg.pinv(x)
    beta = left @ g
    ref = g - x @ beta
    compact = (k[:, tr] - x @ left @ k[:, tr]) @ pt @ c / scale
    weighted = []
    for mode, w in [('uniform', np.ones(n_pairs)),
                    ('fixed_positive', np.exp(rng.normal(0, .4, n_pairs)))]:
        aw = weighted_left_inverse(x, w)
        rw = g - x @ aw @ g
        p = residualizer(x[q])
        weighted.append({
            'reference': mode,
            'reference_weighted_orthogonality_max_abs': float(np.max(np.abs(x.T @ (w[:, None] * rw)))),
            'query_P_change_max_abs': float(np.max(np.abs(p @ (rw[q] - g[q])))),
            'left_inverse_identity_max_abs': float(np.max(np.abs(aw @ x - np.eye(n)))),
        })

    # Complete graph closed form uses only endpoint row sums and global sum.
    endpoint_sum = x.T @ g
    closed_beta = (endpoint_sum - np.sum(g, axis=0, keepdims=True) / (n - 1)) / (n - 2)

    # Query invariance: parameters are fixed; only request rows may change.
    queries = [q, q[::-1], q[:7], q[[2, 2, 5, 0]]]
    queries.extend(np.asarray([j]) for j in q)
    max_request_delta = 0.0
    max_old_batch_delta = 0.0
    for ix in queries:
        delivered = g[ix] - x[ix] @ beta
        max_request_delta = max(max_request_delta, float(np.max(np.abs(delivered - ref[ix]))))
        old = residualizer(x[ix]) @ g[ix]
        max_old_batch_delta = max(max_old_batch_delta, float(np.max(np.abs(old - ref[ix]))))

    # Exact endpoint gauge invariance and endpoint-only null.
    gauge = rng.normal(size=(n, g.shape[1]))
    shifted = g + x @ gauge
    ref_shifted = shifted - x @ left @ shifted
    additive = x @ gauge
    null = additive - x @ left @ additive

    # Two-sided kernel residualization: same PSD feature map on every pair.
    pf = residualizer(x)
    kref = pf @ k @ pf.T
    qt_formula = (k[np.ix_(q, tr)]
                  - x[q] @ left @ k[:, tr]
                  - k[q] @ left.T @ xt.T
                  + x[q] @ left @ k @ left.T @ xt.T)
    train_equivalence = pt @ kref[np.ix_(tr, tr)] @ pt - projected
    cross_equivalence = qt_formula @ pt @ c / scale - ref[q]
    pq = residualizer(x[q])
    # Full-error identity is independently evaluated with random baseline error.
    err = rng.normal(size=(len(q), g.shape[1]))
    full_direct = np.sum((err + ref[q]) ** 2)
    full_decomp = (np.sum((pq @ (err + g[q])) ** 2)
                   + np.sum(((np.eye(len(q)) - pq) @ (err + g[q]) - x[q] @ beta) ** 2))

    checks = {
        'n_endpoints': n,
        'n_reference_pairs': n_pairs,
        'reference_incidence_rank': int(np.linalg.matrix_rank(x)),
        'train_pairs': tr.tolist(),
        'query_pairs': q.tolist(),
        'reference_generation': 'all unordered distinct pairs from endpoint universe',
        'uses_measured_data': False,
        'direct_vs_compact_kernel_max_abs': float(np.max(np.abs(ref - compact))),
        'pinv_vs_complete_graph_closed_form_max_abs': float(np.max(np.abs(beta - closed_beta))),
        'reference_orthogonality_max_abs': float(np.max(np.abs(x.T @ ref))),
        'reference_XTX_formula_max_abs': float(np.max(np.abs(x.T @ x - ((n-2)*np.eye(n) + np.ones((n,n)))))),
        'same_model_request_change_max_abs': max_request_delta,
        'reference_vs_old_batch_max_change': max_old_batch_delta,
        'nonzero_singleton_reference_correction_max_abs': float(np.max(np.abs(ref[q]))),
        'endpoint_gauge_invariance_max_abs': float(np.max(np.abs(ref - ref_shifted))),
        'endpoint_only_null_max_abs': float(np.max(np.abs(null))),
        'query_P_invariance_max_abs': float(np.max(np.abs(pq @ ref[q] - pq @ g[q]))),
        'reference_kernel_min_eigenvalue': float(np.linalg.eigvalsh((kref+kref.T)/2).min()),
        'reference_kernel_symmetry_max_abs': float(np.max(np.abs(kref-kref.T))),
        'reference_kernel_expanded_formula_max_abs': float(np.max(np.abs(qt_formula-kref[np.ix_(q,tr)]))),
        'second_training_projection_equivalence_max_abs': float(np.max(np.abs(train_equivalence))),
        'reference_kernel_prediction_equivalence_max_abs': float(np.max(np.abs(cross_equivalence))),
        'full_error_decomposition_abs': float(abs(full_direct-full_decomp)),
        'weighted_reference_checks': weighted,
        'interpretation': 'P preservation is an identity at fixed model; full-error improvement requires measured heldout evidence.',
    }
    errors = [v for key,v in checks.items()
              if (key.endswith('_max_abs') or key.endswith('_abs'))
              and key != 'nonzero_singleton_reference_correction_max_abs']
    errors.extend(a[key] for a in weighted for key in a if key.endswith('_max_abs'))
    passed = (max(errors) < 1e-10
              and checks['reference_kernel_min_eigenvalue'] > -1e-10
              and checks['nonzero_singleton_reference_correction_max_abs'] > 1e-3
              and checks['reference_vs_old_batch_max_change'] > 1e-3)
    checks['status'] = 'PASS' if passed else 'FAIL'
    args.out.write_text(json.dumps(checks, indent=2)+'\n')
    print(json.dumps({'status': checks['status'], 'largest_identity_error': max(errors),
                      'output': str(args.out)}))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
