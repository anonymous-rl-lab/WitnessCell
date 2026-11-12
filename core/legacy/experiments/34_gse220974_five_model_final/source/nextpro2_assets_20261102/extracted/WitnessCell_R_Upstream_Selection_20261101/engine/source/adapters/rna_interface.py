"""Training-only RNA interface for the four-model comparison.

The original M0 operator is unchanged. Each validation M0 is tuned and fitted
strictly inside that validation split's training subset, including its own
auxiliary identity folds. No outer-query outcomes are stored in these inputs.
"""
from pathlib import Path
import json, sys
import numpy as np

ROOT = Path(__file__).resolve().parent
ORIGINAL = ROOT / 'sf12_original'
sys.path.insert(0, str(ORIGINAL / 'original_code'))
import global_pair_model as gm


def load_npz(path):
    with np.load(path, allow_pickle=False) as f:
        return {k: f[k].copy() for k in f.files}


def outer_inputs(fold):
    common = load_npz(ORIGINAL / 'inputs/common.npz')
    data = load_npz(ORIGINAL / f'inputs/fold_{fold}.npz')
    design = json.loads((ORIGINAL / 'inputs/frozen_original_design.json').read_text())
    names = common['target_names'].astype(str).tolist()
    assert names == sorted(names) and len(names) == 14
    pairs = common['all_pair_names'].astype(str).tolist()
    lookup = {p: i for i, p in enumerate(pairs)}
    train = np.array([lookup[p] for p in data['train_names']], dtype=int)
    query = np.array([lookup[p] for p in data['query_names']], dtype=int)
    assert not set(train) & set(query)
    assert len(train) + len(query) == 55
    response = np.full((55, len(common['genes'])), np.nan)
    response[train] = data['train_y']
    singles = dict(zip(names, data['single_effects']))
    blocks = design['model_panel_geometry']['picked']
    full = gm.make_dataset(pairs, singles, response, blocks)
    mask = common['score_gene_mask'].astype(bool)
    score = gm.make_dataset(pairs, {n: v[mask] for n, v in singles.items()}, response[:, mask], blocks)
    inner = [train[data['inner_fold_id'] == j] for j in range(3)]
    assert np.array_equal(full.x[train], data['train_u'])
    assert np.array_equal(full.x[query], data['query_u'])
    assert mask.sum() == 4995 and np.isnan(full.y[query]).all()
    return common, data, full, score, train, query, inner


def edges_from_u(u):
    assert np.all(np.sum(u != 0, axis=1) == 2)
    return np.array([np.flatnonzero(row) for row in u], dtype=np.int64)


def make_case(common, full, score, train, query, tune_folds, with_truth=False):
    """Build one fit/validation case; no outcome outside train is used in M0."""
    train, query = np.asarray(train, int), np.asarray(query, int)
    # Erase all non-fitting responses even from the numerical Dataset handed to
    # the fitting routine, making the data boundary testable rather than implied.
    yfit = np.full_like(full.y, np.nan); yfit[train] = full.y[train]
    yscore = np.full_like(score.y, np.nan); yscore[train] = score.y[train]
    fit_full = gm.Dataset(full.pairs, full.nodes, yfit, full.additive, full.x, full.pair_kernel, full.blocks)
    fit_score = gm.Dataset(score.pairs, score.nodes, yscore, score.additive, score.x, score.pair_kernel, score.blocks)
    pars = gm.tune_baseline(fit_score, train, tune_folds)
    bt = gm.baseline_predict(fit_full, train, train, pars)
    bq = gm.baseline_predict(fit_full, train, query, pars)
    result = dict(n_nodes=len(common['target_names']),
                  train_edges=edges_from_u(full.x[train]),
                  train_residual=full.y[train] - bt,
                  train_y=full.y[train].copy(), train_base=bt,
                  query_edges=edges_from_u(full.x[query]), query_base=bq,
                  train_names=np.array([full.pairs[i] for i in train]),
                  query_names=np.array([full.pairs[i] for i in query]),
                  train_global_indices=train, query_global_indices=query,
                  baseline_parameters=pars,
                  baseline_tune_folds=[f.tolist() for f in tune_folds])
    if with_truth:
        assert np.isfinite(full.y[query]).all(), 'Validation must belong to allowed outer-training responses'
        result.update(val_edges=result['query_edges'], val_base=bq,
                      val_y=full.y[query].copy())
    return result


def prepare_fold(fold):
    common, data, full, score, tr, q, inner_ids = outer_inputs(fold)
    final = make_case(common, full, score, tr, q, inner_ids)
    inner = []
    for val in inner_ids:
        fit = np.array(sorted(set(tr) - set(val)), int)
        auxiliary = gm.preserved_folds(fit, full.pairs, full.blocks)
        case = make_case(common, full, score, fit, val, auxiliary, with_truth=True)
        assert not set(case['train_global_indices']) & set(case['query_global_indices'])
        inner.append(case)
    return dict(fold=int(fold), common=common, final=final, inner=inner,
                n_nodes=len(common['target_names']), genes=common['genes'],
                selection_mask=common['score_gene_mask'].copy(),
                score_gene_mask=common['score_gene_mask'].copy(),
                target_names=common['target_names'], target_indices=common['target_indices'],
                control=common['control'],
                query_outcomes_available=False)


if __name__ == '__main__':
    import time
    t = time.monotonic()
    d = prepare_fold(0)
    print(json.dumps(dict(fold=d['fold'], nodes=d['n_nodes'], genes=len(d['genes']),
                         selection_genes=int(d['selection_mask'].sum()),
                         train=d['final']['train_residual'].shape,
                         query=d['final']['query_base'].shape,
                         inner=[dict(train=x['train_residual'].shape, val=x['val_y'].shape)
                                for x in d['inner']], seconds=time.monotonic()-t,
                         query_outcomes_available=d['query_outcomes_available']), indent=2))
