"""Training-only M0 comparators and exact combination-axis RNA diagnostics.

No function in fit_baselines reads a held-out test response.  Effects, not raw
counts, are accepted; all outputs remain in the original log-normalized RNA
units.  Endpoint projection acts on CONDITIONS and is unrelated to removing
RNA columns corresponding to the manipulated genes.
"""
from __future__ import annotations
import hashlib
import numpy as np

ALPHAS = (0., .02, .05, .1, .2, .5)
NOISES = (.01, .03, .1, .3, 1.)
MULTI_RIDGES = (.001, .01, .1, 1.)


def sat(z, alpha):
    return z / (1. + alpha * np.abs(z))


def endpoint_basis(incidence, rtol=1e-10):
    """Orthonormal rows spanning the left null space of condition incidence."""
    x = np.asarray(incidence, dtype=np.float64)
    u, s, _ = np.linalg.svd(x, full_matrices=True)
    tol = rtol * max(float(s[0]) if len(s) else 0., 1.)
    rank = int(np.sum(s > tol))
    return u[:, rank:].T


def _ids(x):
    return np.asarray(x, dtype=np.int64)


def _identity_folds(ids, incidence, nfolds=3):
    """Partition identities, never RNA measurements; no stratification by y."""
    ids = _ids(ids)
    ordered = sorted(ids.tolist(), key=lambda i: hashlib.sha256(
        ('20261027|' + ','.join(map(str, np.flatnonzero(incidence[i])))).encode()
    ).hexdigest())
    folds = [np.array(ordered[k::nfolds], dtype=np.int64) for k in range(nfolds)]
    return [f for f in folds if len(f)]


def _m0_tune(e, y, x, fit_ids, non_target):
    """Same alpha/noise/gamma grid and operator as legacy baseline_core.py.

    Only the identity-based internal split is new when this API is given a new
    training subset.  Original M0 saturation and endpoint shrinkage are intact.
    """
    folds = _identity_folds(fit_ids, x)
    additive = x @ e
    dot = np.zeros((len(ALPHAS), len(NOISES)))
    raw2 = np.zeros_like(dot)
    truth2 = np.zeros(len(ALPHAS))
    allowed = set(map(int, fit_ids))
    scored = 0
    for va in folds:
        tr = np.array(sorted(allowed - set(va)), dtype=np.int64)
        if not len(tr):
            raise ValueError('M0 requires at least two fitting conditions')
        xt, xv = x[tr], x[va]
        kt, kv = xt @ xt.T / 2., xv @ xt.T / 2.
        eig, v = np.linalg.eigh(kt)
        left = kv @ v
        operators = [(left / (eig + noise)) @ v.T for noise in NOISES]
        yt, yv = y[tr][:, non_target], y[va][:, non_target]
        for ai, alpha in enumerate(ALPHAS):
            residual = yt - sat(additive[tr][:, non_target], alpha)
            target = yv - sat(additive[va][:, non_target], alpha)
            truth2[ai] += np.sum(target * target)
            for ni, op in enumerate(operators):
                raw = op @ residual
                dot[ai, ni] += np.sum(raw * target)
                raw2[ai, ni] += np.sum(raw * raw)
        scored += yv.size
    gamma = np.divide(dot, raw2, out=np.zeros_like(dot), where=raw2 > 1e-12).clip(0., 1.)
    losses = (truth2[:, None] - 2 * gamma * dot + gamma * gamma * raw2) / scored
    ai, ni = min(np.ndindex(losses.shape), key=lambda ij: (
        losses[ij], ALPHAS[ij[0]], NOISES[ij[1]]))
    return dict(alpha=ALPHAS[ai], noise=NOISES[ni], gamma=float(gamma[ai, ni]),
                oof_mse=float(losses[ai, ni]), internal_folds=[v.tolist() for v in folds],
                candidate_losses=losses.tolist(), candidate_gammas=gamma.tolist(),
                selection='fit-only whole-condition cross-validation')


def _m0_state(e, y, x, fit_ids, pars):
    xt = x[fit_ids]
    residual = y[fit_ids] - sat(xt @ e, pars['alpha'])
    # Primal is exactly equivalent to the old XX'/2 kernel with noise I.
    coef = np.linalg.solve(xt.T @ xt + 2 * pars['noise'] * np.eye(x.shape[1]), xt.T @ residual)
    return dict(single_effects=e, coef=coef, alpha=pars['alpha'], gamma=pars['gamma'])


def predict_state(state, u):
    u = np.asarray(u, dtype=np.float64)
    return sat(u @ state['single_effects'], state['alpha']) + state['gamma'] * (u @ state['coef'])


def _multi_training(e, y, x, fit_ids, triple_u, triple_y):
    t = len(e)
    us = [np.eye(t), x[fit_ids]]
    ys = [e, y[fit_ids]]
    normalizers = [max(float(np.mean(e**2)), 1e-8),
                   max(float(np.mean(y[fit_ids]**2)), 1e-8)]
    weights = [np.full(t, 1. / (t * normalizers[0])),
               np.full(len(fit_ids), 1. / (len(fit_ids) * normalizers[1]))]
    if len(triple_u):
        us.append(triple_u)
        ys.append(triple_y)
        normalizers.append(max(float(np.mean(triple_y**2)), 1e-8))
        weights.append(np.full(len(triple_u), .1 / (len(triple_u) * normalizers[2])))
    return np.concatenate(us), np.concatenate(ys), np.concatenate(weights), normalizers


def _multi_state(e, y, x, fit_ids, triple_u, triple_y, alpha, ridge, gamma=1.):
    u, yy, w, normalizers = _multi_training(e, y, x, fit_ids, triple_u, triple_y)
    r = yy - sat(u @ e, alpha)
    normal = u.T @ (w[:, None] * u)
    # The network's training partner loss is also 0.5 * mean((C e)^2)/pden.
    # C @ pair_incidence == 0, so an endpoint correction has exactly zero
    # derivative for that term. It cannot change this closed-form endpoint fit.
    # Normalize the ridge grid to the weighted design, avoiding an accidental
    # loss of shrinkage when RNA-effect energy normalization is introduced.
    ridge_scale = max(float(np.trace(normal) / len(e)), 1e-12)
    coef = np.linalg.solve(normal + ridge * ridge_scale * np.eye(len(e)),
                           u.T @ (w[:, None] * r))
    return dict(single_effects=e, coef=coef, alpha=float(alpha), gamma=float(gamma),
                training_energy_normalizers=normalizers, ridge_scale=ridge_scale)



def _assert_no_triple_leak(triple_u, heldout_x):
    if not len(triple_u) or not len(heldout_x):
        return
    contains = (triple_u[:, None, :] >= heldout_x[None, :, :]).all(axis=2)
    if contains.any():
        raise ValueError('An auxiliary triple contains a held-out validation/test pair')


def fit_baselines(single_effects, pair_effects, incidence, fit_ids, val_ids, test_ids,
                  triple_u, triple_effects, non_target):
    """Fit original M0 and a multi-order information-matched extension.

    Returns predictions for fit/validation/test and numerical states supporting
    query-only deployment.  M0 uses pairs only, preserving its old operator.
    M0_multi uses the network's energy-normalized order losses with multipliers
    1/1/0.1; alpha/ridge/gamma use its normalized full+0.5*partner validation
    objective on whole validation pairs, never on test y.
    The latter is an explicit new comparator, not a claim of old-M0 equivalence.
    """
    e, y, x = (np.asarray(a, dtype=np.float64) for a in
               (single_effects, pair_effects, incidence))
    tr, va, te = map(_ids, (fit_ids, val_ids, test_ids))
    nt = np.asarray(non_target, dtype=bool)
    tu = np.asarray(triple_u, dtype=np.float64).reshape(-1, x.shape[1])
    ty = np.asarray(triple_effects, dtype=np.float64).reshape(-1, y.shape[1])
    if set(tr) & set(va) or set(tr) & set(te) or set(va) & set(te):
        raise ValueError('Fit, validation, and test identities overlap')
    if len(tu) != len(ty):
        raise ValueError('Triple inputs/targets differ in length')
    if not np.isfinite(y[tr]).all() or (len(va) and not np.isfinite(y[va]).all()):
        raise ValueError('Allowed train/validation labels must be finite')
    _assert_no_triple_leak(tu, x[np.r_[va, te]])
    if not len(va):
        raise ValueError('M0_multi selection requires an explicit whole-combination validation set')

    pars = _m0_tune(e, y, x, tr, nt)
    m0 = _m0_state(e, y, x, tr, pars)
    cval = endpoint_basis(x[va])
    additive_val = x[va] @ e
    # Exactly the same formula as Training.validate in train.py. All response
    # energy normalizers use the permitted validation labels, not the test.
    vden = max(float(np.mean(y[va][:, nt]**2)), 1e-8)
    vpden = max(float(np.mean((cval @ y[va][:, nt])**2)), 1e-8) if len(cval) else 1.
    def val_inner(a, b):
        a, b = a[:, nt], b[:, nt]
        result = np.mean(a * b) / vden
        if len(cval):
            result += .5 * np.mean((cval @ a) * (cval @ b)) / vpden
        return float(result)
    scores = []
    states = []
    for alpha in ALPHAS:
        base = sat(additive_val, alpha)
        target = y[va] - base
        for ridge in MULTI_RIDGES:
            state = _multi_state(e, y, x, tr, tu, ty, alpha, ridge)
            raw = x[va] @ state['coef']
            gamma = np.clip(val_inner(raw, target) / max(val_inner(raw, raw), 1e-30), 0., 1.)
            state['gamma'] = float(gamma)
            err = target - gamma * raw
            score = val_inner(err, err)
            scores.append(dict(alpha=alpha, ridge=ridge, gamma=float(gamma), validation_normalized_score=score))
            states.append(state)
    selected = min(range(len(scores)), key=lambda i: (
        scores[i]['validation_normalized_score'], -scores[i]['ridge'], scores[i]['alpha']))
    states_out = dict(M0=m0, M0_multi=states[selected])
    queries = dict(fit=tr, validation=va, test=te)
    predictions = {name: {split: predict_state(state, x[ids]) for split, ids in queries.items()}
                   for name, state in states_out.items()}
    predictions['additive'] = {split: x[ids] @ e for split, ids in queries.items()}
    return dict(predictions=predictions, states=states_out, parameters=dict(
        M0=pars, M0_multi=dict(**scores[selected], candidate_scores=scores,
            group_loss_multipliers=dict(single=1., pair=1., triple=.1),
            training_energy_normalizers=states[selected]['training_energy_normalizers'],
            ridge_scale=states[selected]['ridge_scale'],
            ridge_units='dimensionless multiplier of trace(weighted normal matrix)/n_TF',
            validation_energy_normalizers=dict(full=vden, partner=vpden),
            fit_pair_ids=tr.tolist(), validation_pair_ids=va.tolist(),
            n_auxiliary_triples=len(tu), selection='whole-validation-combination full_MSE/full_energy + 0.5*partner_MSE/partner_energy',
            triple_reference_caveat='Triple effects use the supplied reference; only five triple-guide controls exist.'),
        test_labels_used=False))


def score(pred, truth, query_incidence, non_target, reference=None):
    """Original-unit errors, signed direction, and exact SSE decomposition."""
    pred, truth = map(lambda a: np.asarray(a, dtype=np.float64), (pred, truth))
    nt = np.asarray(non_target, dtype=bool)
    if pred.shape != truth.shape:
        raise ValueError('Prediction/target shapes differ')
    c = endpoint_basis(query_incidence)
    err = pred - truth
    full = float(np.sum(err[:, nt] ** 2))
    part = float(np.sum((c @ err[:, nt]) ** 2))
    out = dict(full_SSE=full, full_MSE=full / max(err[:, nt].size, 1),
               partner_SSE=part, independent_partner_dimensions=len(c),
               endpoint_SSE=full-part,
               all_gene_SSE=float(np.sum(err**2)),
               manipulated_gene_columns_SSE=float(np.sum(err[:, ~nt]**2)),
               per_condition_non_target_SSE=np.sum(err[:, nt]**2, axis=1).tolist())
    if reference is not None:
        ref = np.asarray(reference, dtype=np.float64)
        r = truth[:, nt] - ref[:, nt]
        delta = pred[:, nt] - ref[:, nt]
        base_full = float(np.sum(r*r))
        pr, pd = c @ r, c @ delta
        base_part = float(np.sum(pr*pr))
        inner = float(np.sum(pr*pd))
        delta2 = float(np.sum(pd*pd))
        # Roundoff from projecting a purely endpoint correction is not a
        # partner direction; do not report huge oracle scalings for 1e-31 energy.
        direction_present = delta2 > 1e-20 * max(float(np.sum(delta*delta)), 1e-30)
        cosine = inner / np.sqrt(base_part*delta2) if direction_present and base_part*delta2 > 0 else 0.
        endpoint_change = (full-part) - (base_full-base_part)
        partner_change = part - base_part
        total_change = full - base_full
        out.update(full_gain_pct=100.*(1.-full/base_full) if base_full else None,
                   partner_gain_pct=100.*(1.-part/base_part) if base_part else None,
                   partner_direction_cosine=float(cosine),
                   partner_correction_energy=delta2,
                   oracle_nonnegative_scalar=float(max(inner/delta2, 0.)) if direction_present else 0.,
                   partner_direction_numerically_present=bool(direction_present),
                   oracle_nonnegative_direction_gain_pct=float(max(cosine, 0.)**2*100.),
                   endpoint_SSE_change=endpoint_change,
                   partner_SSE_change=partner_change,
                   total_SSE_change=total_change,
                   decomposition_numerical_error=float(total_change-endpoint_change-partner_change))
    return out
