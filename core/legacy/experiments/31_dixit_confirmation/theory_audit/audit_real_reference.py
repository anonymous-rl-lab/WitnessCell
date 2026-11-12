#!/usr/bin/env python3
"""Independent real-output reconstruction and conditional sensitivity.

Reconstructs the fixed-reference predictor from saved training states using a
direct linear solve and the complete-graph closed form. Imports no production
prediction code. Evaluation labels are opened only after reconstructed
predictions are written. Sensitivity is conditional on frozen predictions.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def design(pairs, nodes):
    idx = {g:i for i,g in enumerate(nodes)}
    edges = np.array([[idx[g] for g in p.split('+')] for p in pairs])
    x = np.zeros((len(pairs), len(nodes)))
    x[np.arange(len(pairs)), edges[:,0]] = 1
    x[np.arange(len(pairs)), edges[:,1]] = 1
    return x, edges


def pkernel(a, b, gram):
    i,j = a.T
    u,v = b.T
    return .5*(gram[np.ix_(i,u)]*gram[np.ix_(j,v)]
               +gram[np.ix_(i,v)]*gram[np.ix_(j,u)])


def p_residual(x, y):
    return y-x@np.linalg.lstsq(x,y,rcond=None)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model-dir', type=Path, required=True)
    ap.add_argument('--kernel-dir', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    tick = time.monotonic()
    mp = np.load(a.model_dir/'PREDICTIONS.npz')
    kp = np.load(a.kernel_dir/'PREDICTIONS.npz')
    pairs = mp['pairs'].tolist()
    nobs = len(pairs)
    folds = mp['foldids']
    names = ['scalar_saturation','global_selected'] + [f'global_fixed_{v}' for v in (.1,1.,10.)]
    reconstructed = {name:np.full(nobs,np.nan) for name in names}
    checks=[]
    sources=[a.model_dir/'PREDICTIONS.npz',a.kernel_dir/'PREDICTIONS.npz']
    for fold in sorted(np.unique(folds)):
        sp=a.model_dir/f'FITTED_STATE_{fold}.npz'
        fp=a.model_dir/f'FOLD_{fold}.json'
        sources += [sp, fp]
        s=np.load(sp)
        info=json.loads(fp.read_text())
        nodes=s['nodes'].tolist()
        ref=['+'.join(z) for z in itertools.combinations(nodes,2)]
        xr,er=design(ref,nodes)
        xt,et=design(s['fit_pairs'].tolist(),nodes)
        q=np.array(info['query'])
        qp=[pairs[i] for i in q]
        xq,eq=design(qp,nodes)
        e=s['e']; gram=e@e.T/e.shape[1]
        ktt=pkernel(et,et,gram)
        kft=pkernel(er,et,gram)
        kqt=pkernel(eq,et,gram)
        pt=np.eye(len(xt))-xt@np.linalg.pinv(xt)
        kt=pt@ktt@pt
        scale=np.trace(kt)/len(xt)
        n=len(nodes)
        # Compute reference feature centering without enumerating pair rows.
        u,v=et.T
        colsum=gram.sum(axis=0)
        endpoint_rows=.5*(gram[:,u]*colsum[v]+gram[:,v]*colsum[u])-gram[:,u]*gram[:,v]
        global_sum=.5*(colsum[u]*colsum[v]-np.sum(gram[:,u]*gram[:,v],axis=0))
        feature_beta=(endpoint_rows-global_sum[None,:]/(n-1))/(n-2)
        sigma=s['sigma']; pi=s['pi']; alpha=float(s['alpha'])
        def base(x):
            z=np.mean((x@e)[:,pi]/sigma[pi],axis=1)
            return z/(1+alpha*np.abs(z))
        bq=base(xq)
        reconstructed['scalar_saturation'][q]=bq
        selected=float(info['choices']['global'])
        invrows=[]
        for lam in (.1,1.,10.):
            # Direct solve; production uses an eigendecomposition.
            coef=np.linalg.solve(kt/scale+lam*np.eye(len(xt)),s['resid'])
            rawf=(kft-xr@np.linalg.pinv(xt)@ktt)@pt@coef/scale
            rawq=(kqt-xq@np.linalg.pinv(xt)@ktt)@pt@coef/scale
            # Complete graph closed form; production uses pseudoinverse.
            beta=(xr.T@rawf-np.sum(rawf,axis=0,keepdims=True)/(n-1))/(n-2)
            correction=rawq-xq@beta
            prediction=bq+correction.ravel()
            reconstructed[f'global_fixed_{lam}'][q]=prediction
            if lam==selected:
                reconstructed['global_selected'][q]=prediction
            compact=(kqt-xq@np.linalg.pinv(xr)@kft)@pt@coef/scale
            no_enumeration=(kqt-xq@feature_beta)@pt@coef/scale
            Pq=np.eye(len(q))-xq@np.linalg.pinv(xq)
            soft=(np.eye(len(q))-xq@np.linalg.inv(xr.T@xr)@xq.T)@rawq
            # No label arrays or realized query pool in the pointwise function.
            def request(request_pairs):
                xx,ee=design(request_pairs,nodes)
                kk=pkernel(ee,et,gram)
                raw=(kk-xx@np.linalg.pinv(xt)@ktt)@pt@coef/scale
                return base(xx)+(raw-xx@beta).ravel()
            single=np.array([request([p])[0] for p in qp])
            chunk=np.concatenate([request(qp[j:j+4]) for j in range(0,len(qp),4)])
            perm=np.arange(len(qp))[::-1]
            shuffled=request([qp[j] for j in perm])[np.argsort(perm)]
            allreq=request(ref)[[ref.index(p) for p in qp]]
            dup=request(qp+[qp[0]])[:-1]
            swapped=request(['+'.join(p.split('+')[::-1]) for p in qp])
            crossname=f'ref_trainP_fixed_{lam:g}'
            row={
                'fold':int(fold),'lambda':lam,'selected':lam==selected,
                'saved_predictions_max_abs':float(np.max(abs(prediction-mp[f'global_fixed_{lam}'][q]))),
                'independent_PSD_kernel_max_abs':float(np.max(abs(prediction-kp[crossname][q]))),
                'compact_kernel_identity_max_abs':float(np.max(abs(compact-correction))),
                'no_reference_pair_enumeration_max_abs':float(np.max(abs(no_enumeration-correction))),
                'reference_feature_closed_form_max_abs':float(np.max(abs(feature_beta-np.linalg.pinv(xr)@kft))),
                'complete_panel_soft_projection_max_abs':float(np.max(abs(soft-correction))),
                'reference_orthogonality_max_abs':float(np.max(abs(xr.T@(rawf-xr@beta)))),
                'training_orthogonality_max_abs':float(np.max(abs(xt.T@rawf[[ref.index(p) for p in s['fit_pairs'].tolist()]]))),
                'P_query_preservation_max_abs':float(np.max(abs(Pq@correction-Pq@rawq))),
                'singleton_max_abs':float(np.max(abs(single-prediction))),
                'chunks4_max_abs':float(np.max(abs(chunk-prediction))),
                'order_max_abs':float(np.max(abs(shuffled-prediction))),
                'all_universe_request_max_abs':float(np.max(abs(allreq-prediction))),
                'duplicate_request_max_abs':float(np.max(abs(dup-prediction))),
                'endpoint_swap_max_abs':float(np.max(abs(swapped-prediction))),
                'nonzero_singleton_corrections':int(np.count_nonzero(abs(correction)>1e-12)),
            }
            row['PASS']=all(v<1e-11 for k,v in row.items() if k.endswith('_max_abs'))
            checks.append(row)
    assert all(np.isfinite(v).all() for v in reconstructed.values())
    np.savez_compressed(a.out/'INDEPENDENT_PREDICTIONS.npz',pairs=np.array(pairs),foldids=folds,**reconstructed)
    write_json(a.out/'RECONSTRUCTION_FROZEN.json',{
        'prediction_sha256':sha(a.out/'INDEPENDENT_PREDICTIONS.npz'),
        'evaluation_labels_opened':False,
        'reconstruction':'saved training residuals; direct linear solve; complete-graph closed-form reference endpoint removal',
    })
    # Only now access evaluation outcomes; they do not enter any reconstruction.
    el=np.load(a.model_dir/'EVALUATION_LABELS.npz')
    truth=el['truth']; x=el['x']
    assert el['pairs'].tolist()==pairs and np.array_equal(el['foldids'],folds)
    metrics=[]
    foldmetrics=[]
    for name,pred in reconstructed.items():
        total=np.zeros(3)
        for f in sorted(np.unique(folds)):
            q=np.flatnonzero(folds==f)
            err=pred[q]-truth[q]
            pe=p_residual(x[q],err)
            h=err-pe
            ss=np.array([err@err,pe@pe,h@h])
            total+=ss
            foldmetrics.append({'method':name,'fold':int(f),'full_sse':float(ss[0]),'P_sse':float(ss[1]),'endpoint_sse':float(ss[2])})
        metrics.append({'method':name,'full_sse':float(total[0]),'P_sse':float(total[1]),'endpoint_sse':float(total[2])})
    base=next(v for v in metrics if v['method']=='scalar_saturation')
    foldbase={v['fold']:v for v in foldmetrics if v['method']=='scalar_saturation'}
    for row in metrics:
        for metric in ('full_sse','P_sse','endpoint_sse'):
            row[metric+'_gain_vs_saturation']=1-row[metric]/base[metric]
    for row in foldmetrics:
        for metric in ('full_sse','P_sse','endpoint_sse'):
            row[metric+'_gain_vs_saturation']=1-row[metric]/foldbase[row['fold']][metric]

    # Fixed-prediction conditional node multipliers; all assumptions explicit.
    pred=reconstructed['global_selected']
    baseline=reconstructed['scalar_saturation']
    _,edge=design(pairs,nodes)
    rng=np.random.default_rng(2026102021)
    gains=np.empty((2000,2))
    def weighted_sse(w,keep):
        ans=np.zeros((2,2))
        dims=0
        for f in sorted(np.unique(folds)):
            q=np.flatnonzero((folds==f)&keep)
            if not len(q):continue
            rw=np.sqrt(w[q]); xx=x[q]*rw[:,None]
            dims+=len(q)-np.linalg.matrix_rank(xx)
            for j,pp in enumerate((baseline,pred)):
                ee=(pp[q]-truth[q])*rw
                pe=p_residual(xx,ee)
                ans[j]=ans[j]+[ee@ee,pe@pe]
        return ans,int(dims)
    keep=np.ones(nobs,dtype=bool)
    for j in range(len(gains)):
        nw=rng.exponential(1,len(nodes))
        w=nw[edge[:,0]]*nw[edge[:,1]]
        ss,dims=weighted_sse(w,keep)
        gains[j]=1-ss[1]/ss[0]
    np.savez_compressed(a.out/'CONDITIONAL_NODE_MULTIPLIER.npz',gain=gains)
    deletions=[]
    for j,node in enumerate(nodes):
        keep=np.all(edge!=j,axis=1)
        ss,dims=weighted_sse(np.ones(nobs),keep)
        deletions.append({'excluded_endpoint':node,'remaining_pairs':int(keep.sum()),'remaining_P_dimensions':dims,
                          'full_gain_vs_saturation':float(1-ss[1,0]/ss[0,0]),
                          'P_gain_vs_saturation':float(1-ss[1,1]/ss[0,1])})
    conditional={
        'policy':'global_selected',
        'draws':len(gains),'seed':2026102021,
        'distribution':'independent Exp(1) endpoint weights; pair weight is product of its two endpoint weights',
        'columns':['full_gain_vs_saturation','P_gain_vs_saturation'],
        'central_95_percent_range':np.quantile(gains,[.025,.975],axis=0).tolist(),
        'fraction_positive':np.mean(gains>0,axis=0).tolist(),
        'fixed_predictions':True,'weighted_evaluation_projection_recomputed':True,
        'refitting_or_selection_uncertainty_included':False,
        'biological_confidence_interval':False,
        'endpoint_deletion':deletions,
        'deletion_semantics':'remove evaluation records only; predictor/reference unchanged; recompute remaining evaluation projection',
    }
    write_csv(a.out/'REAL_IDENTITY_CHECKS.csv',checks)
    write_csv(a.out/'INDEPENDENT_RESULTS.csv',metrics)
    write_csv(a.out/'INDEPENDENT_FOLD_RESULTS.csv',foldmetrics)
    write_json(a.out/'CONDITIONAL_SENSITIVITY.json',conditional)
    write_csv(a.out/'CONDITIONAL_ENDPOINT_DELETION.csv',deletions)
    errors=[v for row in checks for k,v in row.items() if k.endswith('_max_abs')]
    receipt={
        'status':'PASS' if all(v['PASS'] for v in checks) else 'FAIL',
        'largest_identity_error':max(errors),
        'real_fold_lambda_combinations':len(checks),
        'all_selected_query_corrections_nonzero':all(row['nonzero_singleton_corrections']==15 for row in checks if row['selected']),
        'seconds':time.monotonic()-tick,
        'source_sha256':{str(p):sha(p) for p in sources},
        'boundary':'invariance is conditional on the fitted model and fixed measured-single endpoint universe; P preservation is algebraic',
        'development':True,
    }
    write_json(a.out/'REAL_AUDIT_RECEIPT.json',receipt)
    print(json.dumps({'audit':receipt['status'],'largest_error':max(errors),
                     'selected':next(v for v in metrics if v['method']=='global_selected'),
                     'conditional_range':conditional['central_95_percent_range'],
                     'deletion_positive_full':sum(z['full_gain_vs_saturation']>0 for z in deletions),
                     'deletion_positive_P':sum(z['P_gain_vs_saturation']>0 for z in deletions),
                     'seconds':receipt['seconds']}))
    if receipt['status']!='PASS':raise SystemExit(1)


if __name__=='__main__':
    main()
