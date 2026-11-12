"""Fixed-prediction, fixed-projector endpoint-reweighting sensitivity.

Not a confidence guarantee over experiments or populations. Exponential weights
are shared for the same physical genes across cell contexts. Neither predictors
nor the unweighted query projector are refitted during the sensitivity analysis.
"""
import numpy as np
from shen_adapter import projection_info

def endpoint_sensitivity(contexts,draws=10000,seed=20261107):
    if isinstance(draws,bool) or not isinstance(draws,int) or draws<=0:
        raise ValueError('draws must be positive integer')
    union=sorted({g for x in contexts.values() for g in x['genes']});lookup={g:i for i,g in enumerate(union)}
    stats={}; samples={k:[] for k in contexts}
    for name,x in contexts.items():
        y,base,pred=(np.asarray(x[k],float) for k in ('y','base','pred'))
        if y.shape!=base.shape or y.shape!=pred.shape or y.shape!=(len(x['pairs']),2):raise ValueError('shape')
        if not all(np.isfinite(v).all() for v in (y,base,pred)):raise ValueError('nonfinite')
        _,U,_=projection_info(x['genes'],x['pairs'])
        e=y-base;f=y-pred;pe=e-U@(U.T@e);pf=f-U@(U.T@f)
        ix=np.array([[lookup[a],lookup[b]] for a,b in x['pairs']],int)
        stats[name]=(ix,(pe*pe).sum(1),(pf*pf).sum(1))
    rng=np.random.default_rng(seed)
    for offset in range(0,draws,256):
        w=rng.exponential(size=(min(256,draws-offset),len(union)))
        for name,(ix,a,b) in stats.items():
            pairw=w[:,ix[:,0]]*w[:,ix[:,1]];den=pairw@a;num=pairw@b
            z=np.full(len(den),np.nan);ok=den>1e-24;z[ok]=1-num[ok]/den[ok]
            samples[name].append(z)
    samples={k:np.concatenate(v) for k,v in samples.items()}
    def summarize(z):
        z=z[np.isfinite(z)]
        return dict(valid_draws=int(len(z)),central95=None if not len(z) else np.quantile(z,[.025,.975]).tolist())
    result={k:summarize(v) for k,v in samples.items()}
    matrix=np.array(list(samples.values()));valid=np.isfinite(matrix).all(0)
    return {'method':'Exp(1) endpoint products; shared gene weights across contexts; fixed P and predictions',
            'draws':draws,'seed':seed,'per_context':result,
            'equal_context_mean':summarize(matrix[:,valid].mean(0)),
            'conditional_sensitivity_only':True,'population_confidence_claim':False}
