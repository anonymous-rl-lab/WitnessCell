"""Exact archived SF12 adapter; new coefficients, no query-response access.

The archived rho and amplitude are frozen. The four RNA channels and two
GenePT axes are byte-for-byte implemented in legacy_shared_core.py. This
module adds shape/identity checks and full-RNA target-coordinate bookkeeping.
The historical model learns from singles and AB; it does not natively use ABC.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import numpy as np
import legacy_shared_core as legacy

RHO = legacy.RHO
AMPLITUDE = legacy.AMPLITUDE


def pair_edges(u):
    u = np.asarray(u, dtype=np.float64)
    if u.ndim != 2 or not np.isfinite(u).all():
        raise ValueError('Expected finite pair-incidence matrix')
    if not np.all((u == 0) | (u == 1)) or not np.all(u.sum(1) == 2):
        raise ValueError('Archived SF12 accepts exactly two distinct targets')
    return np.stack([np.flatnonzero(row) for row in u]).astype(np.int64)


def fit_predict(single_effects, knowledge, train_u, train_y, train_base,
                query_u, query_base, non_target, *, knowledge_nodes=None,
                nodes=None):
    """Return full-coordinate predictions and sufficient fitted state.

    Inputs contain training responses only. query_base is already frozen from
    permitted training labels. Every RNA coordinate not selected by non_target
    stays equal to supplied M0, matching historical Dixit SF12 evaluation.
    """
    e = np.asarray(single_effects, dtype=np.float64)
    k = np.asarray(knowledge)
    ut, uq = np.asarray(train_u, float), np.asarray(query_u, float)
    yt, bt, bq = map(lambda a: np.asarray(a, dtype=np.float64),
                     (train_y, train_base, query_base))
    nt = np.asarray(non_target, dtype=bool)
    if e.ndim != 2 or k.ndim != 2 or len(k) != len(e):
        raise ValueError('Single/knowledge target dimension mismatch')
    if nt.shape != (e.shape[1],) or not nt.any():
        raise ValueError('RNA mask mismatch')
    if yt.shape != bt.shape or yt.shape != (len(ut), e.shape[1]):
        raise ValueError('Training expression shape mismatch')
    if bq.shape != (len(uq), e.shape[1]) or ut.shape[1] != len(e) or uq.shape[1] != len(e):
        raise ValueError('Query expression or target dimension mismatch')
    if not all(np.isfinite(a).all() for a in (e, k, ut, uq, yt, bt, bq)):
        raise ValueError('Nonfinite allowed input')
    if (knowledge_nodes is None) != (nodes is None):
        raise ValueError('Supply both identity arrays or neither')
    if nodes is not None and not np.array_equal(np.asarray(nodes), np.asarray(knowledge_nodes)):
        raise ValueError('Knowledge identity order mismatch')
    et, eq = pair_edges(ut), pair_edges(uq)
    # The historical caller centers the allowed single panel before features().
    # features() itself centers knowledge, not RNA: both operations are needed.
    centered = e[:, nt] - e[:, nt].mean(axis=0, keepdims=True)
    ft, fq = legacy.features(centered, k, et, eq)
    state = legacy.fit_fixed_m0((yt - bt)[:, nt], ut, ft,
                                rho=RHO, amplitude=AMPLITUDE)
    pred_nt, delta_nt = legacy.predict(state, bq[:, nt], uq, fq)
    fitted_nt, fit_delta_nt = legacy.predict(state, bt[:, nt], ut, ft)
    prediction = bq.copy(); prediction[:, nt] = pred_nt
    fit_prediction = bt.copy(); fit_prediction[:, nt] = fitted_nt
    delta = np.zeros_like(bq); delta[:, nt] = delta_nt
    fit_delta = np.zeros_like(bt); fit_delta[:, nt] = fit_delta_nt
    info = dict(method='old12_fixed', rho=RHO, amplitude=AMPLITUDE,
                coefficients_refit=True, hyperparameters_refit=False,
                trained_coefficients=12,
                training_double_conditions=len(ut), training_genes=int(nt.sum()),
                independent_edge_contrasts=state['independent_edge_contrasts'],
                projected_feature_rank=state['projected_feature_rank'],
                single_centering='per-RNA mean over all allowed measured singles',
                knowledge='fixed archived GenePT; two axes from allowed target identities',
                combination_feature_scaling='computed only on this task training pairs',
                supervision='single-derived features and AB responses; ABC unused by archived SF12',
                query_response_received=False,
                legacy_source_sha256=hashlib.sha256(
                    (Path(__file__).parent/'legacy_shared_core.py').read_bytes()).hexdigest())
    return dict(prediction=prediction, delta=delta,
                fit_prediction=fit_prediction, fit_delta=fit_delta,
                state=state, metadata=info)


def run_task(name):
    """Fit one frozen task and save predictions before any evaluation read."""
    import json
    from data_prepare import load_task, baseline, read_npz, save_npz, save_json, sha
    root=Path(__file__).resolve().parent
    common,train=load_task(name,include_evaluation=False)
    knowledge=read_npz(root/'data/knowledge.npz')
    bt=baseline(common,train,train['u'])
    bq=baseline(common,train,common['query_u'])
    fitted=fit_predict(common['single_effects'],knowledge['embedding'],
                       train['u'],train['y'],bt,common['query_u'],bq,
                       common['non_target'],knowledge_nodes=knowledge['nodes'],
                       nodes=common['nodes'])
    output=root/'baselines'/name;output.mkdir(parents=True,exist_ok=True)
    save_npz(output/'predictions.npz',M0_fit=bt,M0_query=bq,
             old12_fit=fitted['fit_prediction'],old12_query=fitted['prediction'],
             theta=fitted['state']['theta'],shift=fitted['state']['shift'],
             gram=fitted['state']['gram'],rhs=fitted['state']['rhs'])
    report=dict(task=name,**fitted['metadata'],
                prediction_sha256=sha(output/'predictions.npz'),
                evaluation_response_opened=False,
                code_sha256={p.name:sha(p) for p in [Path(__file__),root/'legacy_shared_core.py']},
                input_sha256={str(p.relative_to(root)):sha(p) for p in
                   [root/'data/common.npz',root/f'data/{name}_train.npz',root/'data/knowledge.npz']})
    save_json(output/'FROZEN.json',report)
    print(json.dumps(dict(task=name,method='old12_fixed',frozen=True,
                         rank=report['projected_feature_rank'],
                         evaluation_response_opened=False)),flush=True)
    return report


if __name__=='__main__':
    import argparse
    from data_prepare import TASKS, save_json
    ap=argparse.ArgumentParser();ap.add_argument('--task',choices=('all',*TASKS),default='all')
    args=ap.parse_args()
    chosen=TASKS if args.task=='all' else (args.task,)
    summaries=[run_task(name) for name in chosen]
    if args.task=='all':
        save_json(Path(__file__).resolve().parent/'OLD12_ALL_FROZEN.json',
                  dict(tasks=summaries,evaluation_response_opened=False,
                       hyperparameter_selection='archived fixed rho/amplitude; coefficient refits only'))
