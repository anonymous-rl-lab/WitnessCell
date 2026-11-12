"""Global complete-baseline plus endpoint-orthogonal symmetric-product predictor.

All fitting functions explicitly receive allowed training indices. Readout and
query-response arrays are never used to construct the endpoint/pair kernels.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import hashlib
from dataclasses import dataclass
import numpy as np

ALPHAS=(0.,.02,.05,.1,.2,.5)
NOISES=(.01,.03,.1,.3,1.)
LAMBDAS=(.1,1.,10.)


def sat(z,alpha): return z/(1.+alpha*np.abs(z))

def digest(edges):
    return hashlib.sha256(('20261017|'+'|'.join(sorted(edges))).encode()).hexdigest()


def preserved_folds(allowed, pairs, blocks, k=3):
    """Identity-only folds; complete original blocks remain intact."""
    allowed=set(map(int,allowed)); lookup={p:i for i,p in enumerate(pairs)}
    available=[];used=set()
    for b in blocks:
        ids=[lookup[p] for p in b]
        if set(ids)<=allowed:
            available.append((digest(b),ids));used.update(ids)
    groups=[ids for _,ids in sorted(available)]
    singles=sorted(allowed-used,key=lambda i:digest([pairs[i]]))
    groups.extend([[i] for i in singles])
    folds=[[] for _ in range(k)]
    for group in groups:
        f=min(range(k),key=lambda j:(len(folds[j]),j));folds[f].extend(group)
    return [np.array(sorted(f),int) for f in folds if f]


def incidence(pairs,nodes):
    lut={g:i for i,g in enumerate(nodes)}
    x=np.zeros((len(pairs),len(nodes)))
    for i,p in enumerate(pairs):
        for g in p.split('+'):x[i,lut[g]]=1.
    return x


def symmetric_kernel(edge_indices, endpoint_gram):
    a,b=edge_indices.T
    return .5*(endpoint_gram[np.ix_(a,a)]*endpoint_gram[np.ix_(b,b)]+endpoint_gram[np.ix_(a,b)]*endpoint_gram[np.ix_(b,a)])


@dataclass
class Dataset:
    pairs:list
    nodes:list
    y:np.ndarray
    additive:np.ndarray
    x:np.ndarray
    pair_kernel:np.ndarray
    blocks:list


def make_dataset(pairs,single_effects,y,blocks):
    nodes=sorted(single_effects)
    e=np.stack([single_effects[g] for g in nodes]).astype(float)
    x=incidence(pairs,nodes)
    lut={g:i for i,g in enumerate(nodes)}
    edge=np.array([[lut[g] for g in p.split('+')] for p in pairs])
    gram=e@e.T/e.shape[1]
    return Dataset(list(pairs),nodes,np.asarray(y,float),x@e,x,symmetric_kernel(edge,gram),blocks)


def tune_baseline(d,tr,folds):
    """Tune alpha/noise and closed-form gamma using only training-internal OOF."""
    tr=np.asarray(tr,int); allowed=set(tr)
    if set(np.concatenate(folds))!=allowed or sum(map(len,folds))!=len(tr):
        raise ValueError('Baseline folds do not partition allowed training')
    dot=np.zeros((len(ALPHAS),len(NOISES)));raw2=np.zeros_like(dot);truth2=np.zeros(len(ALPHAS))
    scored=0
    for va in folds:
        fit=np.array(sorted(allowed-set(va)),int)
        if not len(fit):raise ValueError('Empty baseline training fold')
        xt,xv=d.x[fit],d.x[va]
        kt=xt@xt.T/2.; kv=xv@xt.T/2.
        eig,V=np.linalg.eigh(kt)
        left=kv@V
        operators=[(left/(eig+noise))@V.T for noise in NOISES]
        yt,yv=d.y[fit],d.y[va]
        for ai,alpha in enumerate(ALPHAS):
            r=yt-sat(d.additive[fit],alpha)
            target=yv-sat(d.additive[va],alpha)
            truth2[ai]+=np.sum(target*target)
            for ni,op in enumerate(operators):
                raw=op@r
                dot[ai,ni]+=np.sum(raw*target)
                raw2[ai,ni]+=np.sum(raw*raw)
        scored+=yv.size
    gammas=np.divide(dot,raw2,out=np.zeros_like(dot),where=raw2>1e-12).clip(0,1)
    losses=(truth2[:,None]-2*gammas*dot+gammas*gammas*raw2)/scored
    ai,ni=min(np.ndindex(losses.shape),key=lambda t:(losses[t],ALPHAS[t[0]],NOISES[t[1]]))
    sat_ai=int(np.argmin(truth2))
    return {'alpha':ALPHAS[ai],'noise':NOISES[ni],'gamma':float(gammas[ai,ni]),
            'oof_mse':float(losses[ai,ni]),'saturation_only_alpha':ALPHAS[sat_ai],
            'saturation_only_oof_mse':float(truth2[sat_ai]/scored),
            'candidate_losses':losses.tolist(),'candidate_gammas':gammas.tolist(),
            'internal_fold_sizes':[len(f) for f in folds]}


def baseline_predict(d,tr,query,pars):
    tr=np.asarray(tr,int);q=np.asarray(query,int)
    xt,xq=d.x[tr],d.x[q]
    r=d.y[tr]-sat(d.additive[tr],pars['alpha'])
    # Same finite-sample operator and gamma bounds as the formal v14 backbone.
    raw=np.linalg.solve(xt@xt.T/2.+pars['noise']*np.eye(len(tr)),xt@xq.T/2.).T@r
    return sat(d.additive[q],pars['alpha'])+pars['gamma']*raw


def projected_blocks(xtr,xq,ktt,kqt,kqq):
    plus=np.linalg.pinv(xtr)
    P=np.eye(len(xtr))-xtr@plus
    P=(P+P.T)/2
    A=xq@plus
    tt=P@ktt@P;tt=(tt+tt.T)/2
    qt=(kqt-A@ktt)@P
    qq=kqq-A@kqt.T-kqt@A.T+A@ktt@A.T;qq=(qq+qq.T)/2
    trace=float(np.trace(tt));raw_trace=float(np.trace(ktt))
    zero=trace<=1e-12*max(abs(raw_trace),1e-300)
    scale=trace/len(xtr) if not zero else 0.
    if zero:tt=np.zeros_like(tt);qt=np.zeros_like(qt);qq=np.zeros_like(qq)
    else:tt/=scale;qt/=scale;qq/=scale
    audit={'training_trace_raw':raw_trace,'training_trace_projected':trace,'normalization_scale':scale,
           'zero_trace':bool(zero),'projection_annihilation':float(np.max(np.abs(P@xtr))),
           'training_orthogonality':float(np.max(np.abs(xtr.T@tt))),
           'query_self_min':float(np.min(np.diag(qq))) if len(qq) else None}
    return tt,qt,qq,audit


def pair_predictions(d,tr,query,base_pars,lambdas=LAMBDAS):
    tr=np.asarray(tr,int);q=np.asarray(query,int)
    ktt,kqt,kqq,audit=projected_blocks(d.x[tr],d.x[q],d.pair_kernel[np.ix_(tr,tr)],d.pair_kernel[np.ix_(q,tr)],d.pair_kernel[np.ix_(q,q)])
    base=baseline_predict(d,tr,q,base_pars)
    resid=d.y[tr]-baseline_predict(d,tr,tr,base_pars)
    eig,V=np.linalg.eigh(ktt)
    left=kqt@V
    preds={float(lam):base+((left/(eig+lam))@V.T)@resid for lam in lambdas}
    audit['kernel_min_eigenvalue']=float(eig.min())
    return base,preds,audit


def block_contrast(values,block_ids,x):
    # The +/- pattern is derived exclusively from the four-edge incidence.
    xb=x[block_ids]
    for sign in ((1,-1,-1,1),(1,-1,1,-1),(1,1,-1,-1)):
        sign=np.array(sign,float)
        if np.max(np.abs(sign@xb))<1e-10:return sign@values[block_ids],sign
    raise ValueError('Not a complete four-edge rectangle')


def available_blocks(d,allowed):
    allowed=set(map(int,allowed)); lut={p:i for i,p in enumerate(d.pairs)}
    return [[lut[p] for p in b] for b in d.blocks if {lut[p] for p in b}<=allowed]


def fit_outer(d,tr,query,inner_folds):
    tr=np.asarray(tr,int);query=np.asarray(query,int)
    pars=tune_baseline(d,tr,inner_folds)
    ss={lam:0. for lam in LAMBDAS};cnt=0; ssfull={lam:0. for lam in LAMBDAS};nf=0
    details=[];lut={p:i for i,p in enumerate(d.pairs)}
    for va in inner_folds:
        fit=np.array(sorted(set(tr)-set(va)),int)
        auxiliary=preserved_folds(fit,d.pairs,d.blocks)
        innerpars=tune_baseline(d,fit,auxiliary)
        _,pred,audit=pair_predictions(d,fit,va,innerpars)
        blocks=available_blocks(d,va);local={j:i for i,j in enumerate(va)}
        for b in blocks:
            true,sign=block_contrast(d.y,b,d.x)
            loc=[local[j] for j in b]
            for lam in LAMBDAS:
                delta=sign@pred[lam][loc]-true
                ss[lam]+=float(np.mean(delta*delta))
            cnt+=1
        for lam in LAMBDAS:ssfull[lam]+=float(np.sum((pred[lam]-d.y[va])**2))
        nf+=d.y[va].size
        details.append({'validation_pairs':[d.pairs[i] for i in va], 'training_pairs':[d.pairs[i] for i in fit],
                        'M0':innerpars,'pair_kernel':audit,'rectangle_blocks':len(blocks),
                        'auxiliary_M0_folds':[[d.pairs[i] for i in f] for f in auxiliary]})
    scores={lam:(ss[lam]/cnt if cnt else ssfull[lam]/nf) for lam in LAMBDAS}
    lam=min(LAMBDAS,key=lambda v:(scores[v],-v))
    m0,m1,audit=pair_predictions(d,tr,query,pars,lambdas=(lam,))
    saturation=sat(d.additive[query],pars['saturation_only_alpha'])
    result={'M0':m0,'M1':m1[lam],'additive':d.additive[query].copy(),'saturation_only':saturation}
    meta={'M0':pars,'lambda':lam,'lambda_scores':scores,'lambda_selection_metric':'rectangle_C_MSE' if cnt else 'fallback_full_response_MSE',
          'lambda_selection_blocks':cnt,'pair_kernel':audit,'inner_details':details}
    return result,meta
