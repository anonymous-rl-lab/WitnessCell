"""Small, symmetric empirical-case transfer features; no fitted parameters.

The caller must refit calibration and M0 with every evaluation combination
excluded before constructing ``resid``. Removing a query only from retrieval
does not remove information it may already have contributed to preprocessing.

This module uses raw M0 residuals, not purified biological interactions. They
can include endpoint error, saturation misspecification and measurement noise.
The resulting branch must be assessed with the independent partner contrast.
"""
import numpy as np


FEATURE_NAMES = (
    'case_constant',
    'replacement_difference',
    'anchor_x_replacement_difference',
    'replacement_sum_x_difference',
    'anchor_x_replacement_sum_x_difference',
)


def _edges(values, n):
    e = np.asarray(values, dtype=int).reshape(-1, 2)
    if len(e):
        if np.any(e < 0) or np.any(e >= n) or np.any(e[:, 0] == e[:, 1]):
            raise ValueError('Invalid unordered pair identity')
        e = np.sort(e, axis=1)
        if len({tuple(x) for x in e}) != len(e):
            raise ValueError('Duplicate unordered pairs are not independent cases')
    return e


def _residual_mapping(rows, total, seed):
    """Map case identity to response row inside one anchor's neighbor list."""
    mapping = np.arange(total)
    ix = np.asarray(rows, dtype=int)
    if seed is None or len(ix) < 2:
        return mapping
    rng = np.random.default_rng(int(seed))
    for _ in range(10000):
        shuffled = rng.permutation(ix)
        if np.all(shuffled != ix):
            mapping[ix] = shuffled
            return mapping
    raise RuntimeError('Could not generate fixed-point-free residual assignment')


def feature(singles, train_edges, resid, query_edges, eligiblemask,
            perm_seed=None, singles_centered=False, knowledge_node=None,
            weighting='similarity'):
    """Return H[query, RNA, 5] and traceable source-case metadata.

    ``singles`` should be fold-local calibrated single means [node, RNA]. If
    supplying the archived exact centered singles instead, set
    ``singles_centered=True`` to preserve their values without re-centering.

    For anchor a, source partner b and new partner c, d=T_c-T_b and the five
    fixed multipliers are [1,d,T_a*d,(T_c+T_b)*d,T_a*(T_c+T_b)*d]. T uses the
    archived shared12 panel RMS scaling/floor and tanh. The default case weight
    is proportional to exp(.5*cos(centered_single_b,centered_single_c)
    +.5*cos(knowledge_b,knowledge_c)), with temperature fixed at one. Weights
    depend on the replacement identity but are shared across RNA. The two anchor
    predictions receive exactly 1/2 each; an unsupported anchor contributes zero
    without rescaling. ``weighting='uniform'`` is an additive sanity control,
    not a strong partner-prediction reference.

    With perm_seed, observed residual vectors are deranged WITHIN each anchor's
    eligible neighbor list, with an anchor-specific seed. Graph, per-anchor
    unweighted case mean, response distribution, fixed replacement weights,
    single features and quality mask stay fixed. The weighted case mean can
    change, which is part of breaking the response-to-partner correspondence.
    Degree-one anchors cannot support a nontrivial permutation and are reported.
    """
    singles = np.asarray(singles, dtype=float)
    resid = np.asarray(resid, dtype=float)
    if singles.ndim != 2:
        raise ValueError('singles must have shape [node, RNA]')
    n, g = singles.shape
    te, qe = _edges(train_edges, n), _edges(query_edges, n)
    eligible = np.asarray(eligiblemask, dtype=bool)
    if resid.shape != (len(te), g) or eligible.shape != (len(te),):
        raise ValueError('Case response, pair identity and eligibility disagree')
    if not np.isfinite(singles).all() or not np.isfinite(resid).all():
        raise ValueError('Nonfinite source measurements')
    if set(map(tuple, te)).intersection(map(tuple, qe)):
        raise ValueError('A query combination is present in the support case set')

    centered = singles if singles_centered else singles - singles.mean(axis=0)
    if weighting not in ('similarity', 'uniform'):
        raise ValueError('weighting must be similarity or uniform')
    if weighting == 'similarity':
        if knowledge_node is None:
            raise ValueError('Fixed gene knowledge is required for similarity weighting')
        knowledge = np.asarray(knowledge_node, dtype=float)
        if knowledge.ndim != 2 or knowledge.shape[0] != n or not np.isfinite(knowledge).all():
            raise ValueError('knowledge_node must contain finite [node, feature] embeddings')
        def cosine_gram(x):
            norms = np.linalg.norm(x, axis=1)
            normalized = np.divide(x, norms[:, None], out=np.zeros_like(x),
                                   where=norms[:, None] > 0)
            return np.clip(normalized @ normalized.T, -1., 1.)
        similarity = .5 * cosine_gram(centered) + .5 * cosine_gram(knowledge)
    else:
        similarity = np.zeros((n, n), dtype=float)
    scale = np.sqrt(np.mean(centered * centered, axis=0))
    positive = scale[scale > 1e-12]
    floor = max(1e-10, float(np.median(positive)) * 1e-3) if len(positive) else 1e-10
    scale = np.maximum(scale, floor)
    t = np.tanh(centered / scale)
    result = np.zeros((len(qe), g, len(FEATURE_NAMES)), dtype=float)
    traces = []
    degrees = np.zeros(n, dtype=int)
    for a, b in te[eligible]:
        degrees[a] += 1
        degrees[b] += 1
    mappings = {}
    for anchor in range(n):
        rows = np.flatnonzero(eligible & np.any(te == anchor, axis=1))
        seed = None if perm_seed is None else int(perm_seed) + 1000003 * anchor
        mappings[anchor] = _residual_mapping(rows, len(te), seed)

    for qi, (left, right) in enumerate(qe):
        trace = {'query': [int(left), int(right)], 'anchors': []}
        for anchor, replacement in ((left, right), (right, left)):
            rows = np.flatnonzero(eligible & np.any(te == anchor, axis=1))
            mapping = mappings[int(anchor)]
            side = np.zeros((g, len(FEATURE_NAMES)), dtype=float)
            partners = [int(te[row, 1] if te[row, 0] == anchor else te[row, 0])
                        for row in rows]
            logits = similarity[partners, replacement]
            weights = np.exp(logits - logits.max()) if len(rows) else np.empty(0)
            if len(rows):
                weights /= weights.sum()
            for row, old_partner, weight in zip(rows, partners, weights):
                d = t[replacement] - t[old_partner]
                summed = t[replacement] + t[old_partner]
                u = np.column_stack((np.ones(g), d, t[anchor] * d,
                                     summed * d, t[anchor] * summed * d))
                side += weight * resid[mapping[row], :, None] * u
            result[qi] += .5 * side
            trace['anchors'].append({
                'anchor': int(anchor), 'new_partner': int(replacement),
                'case_rows': rows.tolist(), 'old_partners': partners,
                'residual_source_rows': mapping[rows].tolist(),
                'case_weights': weights.tolist(),
                'anchor_weight': .5, 'missing_anchor': len(rows) == 0,
                'permuted_case_count': int(np.sum(mapping[rows] != rows)),
                'permutation_supported': len(rows) >= 2,
            })
        trace['feature_rms'] = np.sqrt(np.mean(result[qi] ** 2, axis=0)).tolist()
        traces.append(trace)

    return result, {
        'feature_names': list(FEATURE_NAMES),
        'eligible_case_count': int(eligible.sum()),
        'eligible_case_rows': np.flatnonzero(eligible).tolist(),
        'eligible_degree': degrees.tolist(),
        'anchor_residual_mappings': {
            str(a): {'case_rows': np.flatnonzero(eligible & np.any(te == a, axis=1)).tolist(),
                     'source_rows': mappings[a][eligible & np.any(te == a, axis=1)].tolist()}
            for a in range(n)},
        'perm_seed': None if perm_seed is None else int(perm_seed),
        'permutation_nontrivial': bool(perm_seed is not None and np.any(degrees >= 2)),
        'degree_one_anchor_count': int(np.sum(degrees == 1)),
        'permutation_scope': 'eligible response rows within each anchor; unweighted case mean preserved',
        'weighting': weighting,
        'weighting_temperature': 1.,
        'weighting_information': 'fold-local centered singles and fixed full gene knowledge; no combination labels',
        'scale_floor': floor,
        'query_support': traces,
        'zero_residual_produces_zero': bool(np.all(result == 0)) if not np.any(resid) else None,
        'source': 'raw fold-local M0 residual; not isolated biological interaction',
        'no_support_policy': 'zero per missing anchor, fixed 1/2 weights',
    }


features = feature


def smoke():
    """Deterministic implementation probe; not evidence of real prediction."""
    rng = np.random.default_rng(20261026021)
    s = rng.normal(size=(5, 13))
    knowledge = rng.normal(size=(5, 64))
    te = np.array([[0, 1], [0, 2], [1, 3], [2, 3], [0, 4]])
    qe = np.array([[0, 3], [1, 2]])
    r = rng.normal(size=(len(te), s.shape[1]))
    mask = np.array([True, True, True, True, False])
    kwargs = dict(knowledge_node=knowledge)
    h, info = feature(s, te, r, qe, mask, **kwargs)
    swapped, _ = feature(s, te[:, ::-1], r, qe[:, ::-1], mask, **kwargs)
    singles = np.concatenate([feature(s, te, r, q[None], mask, **kwargs)[0] for q in qe])
    zero, _ = feature(s, te, np.zeros_like(r), qe, mask, **kwargs)
    empty, _ = feature(s, te, r, qe, np.zeros(len(te), bool), **kwargs)
    wrong, wi = feature(s, te, r, qe, mask, perm_seed=20261026021, **kwargs)
    wrong2, _ = feature(s, te, r, qe, mask, perm_seed=20261026021, **kwargs)
    # Manual uniform-control first-column calculation for the two anchors.
    uniform, _ = feature(s, te, r, qe, mask, weighting='uniform')
    manual = []
    for a, c in qe:
        sides = []
        for anchor in (a, c):
            rows = np.flatnonzero(mask & np.any(te == anchor, axis=1))
            sides.append(r[rows].mean(0) if len(rows) else np.zeros(r.shape[1]))
        manual.append(.5 * (sides[0] + sides[1]))
    checks = dict(exchange=float(np.max(abs(h - swapped))),
                  single_batch=float(np.max(abs(h - singles))),
                  zero_response=float(np.max(abs(zero))),
                  no_support=float(np.max(abs(empty))),
                  uniform_case_reference=float(np.max(abs(uniform[:, :, 0] - manual))),
                  deterministic_wrong=float(np.max(abs(wrong - wrong2))),
                  wrong_weighted_case_change=float(np.max(abs(wrong[:, :, 0] - h[:, :, 0]))),
                  wrong_fixed_points=int(sum(
                      np.sum(np.array(v['case_rows']) == np.array(v['source_rows']))
                      for v in wi['anchor_residual_mappings'].values() if len(v['case_rows']) >= 2)),
                  wrong_changes_features=bool(np.max(abs(wrong - h)) > 1e-6))
    assert all(checks[k] < 1e-12 for k in ('exchange', 'single_batch', 'zero_response',
                                         'no_support', 'uniform_case_reference', 'deterministic_wrong'))
    assert checks['wrong_fixed_points'] == 0 and checks['wrong_changes_features']
    # The fixed replacement weighting must itself support a nonadditive output;
    # uniform anchor averaging is an exactly additive negative control.
    sr = rng.normal(size=(7, 13)); kr = rng.normal(size=(7, 64))
    tr = np.array([[0,4],[0,5],[1,4],[1,6],[2,4],[2,5],[3,5],[3,6]])
    qr = np.array([[0,2],[0,3],[1,2],[1,3]])
    rr = rng.normal(size=(len(tr), 13)); mr = np.ones(len(tr), bool)
    hr, _ = feature(sr, tr, rr, qr, mr, knowledge_node=kr)
    hu, _ = feature(sr, tr, rr, qr, mr, weighting='uniform')
    c = np.array([1.,-1.,-1.,1.])
    checks['scalar_case_rectangle_norm'] = float(np.linalg.norm(c @ hr[:,:,0]))
    checks['uniform_rectangle_norm'] = float(np.linalg.norm(c @ hu[:,:,0]))
    assert checks['scalar_case_rectangle_norm'] > 1e-6
    assert checks['uniform_rectangle_norm'] < 1e-12
    return dict(status='PASS', checks=checks,
                scope='synthetic indexing/linearity/symmetry probe; no model training')


if __name__ == '__main__':
    import json
    print(json.dumps(smoke(), indent=2))
