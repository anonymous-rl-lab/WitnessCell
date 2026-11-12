"""Shape-portable PairPick-TargetRank ranking layer.

This preserves the 20261101 algorithm and hyperparameter grid. It does NOT
build/retrain R-G, graph nuisance fits, GEARS, Scouter or PRESAGE. Those
upstream fits must be regenerated using only the allowed training subset
on each new dataset. No evaluation labels should be provided.

Input mapping d:
 X[P,T] binary incidence (exactly two distinct intervention identities/row),
 single_effects[T,G], objectives[T,G], single_utility[T,T],
 RG[P,T], Y[P,T] with outer-query rows NaN,
 target_names[T], pair_names[P], train[], query[].
 cases: train/query ids, train_blocks, train/query_RG, train/query_graph,
        train/query_raw; one case per inner validation block.
 outer: same structure, with train=d.train and query=d.query.

An intervention identity is (gene,mode), not merely gene, in mixed a/i data.
Objectives and all gene selection must be established using permitted inputs.
"""
from __future__ import annotations
import numpy as np

def validate_inputs(d, cases, outer):
    x = np.asarray(d['X'])
    if x.ndim != 2 or x.shape[0] < 1 or x.shape[1] < 2:
        raise ValueError('X must have shape [pairs, intervention identities]')
    p, t = x.shape
    if not np.isfinite(x).all() or not np.isin(x, [0,1]).all() or not np.all(x.sum(1)==2):
        raise ValueError('Each pair must contain exactly two distinct binary interventions')
    if len(d['target_names'])!=t or len(set(map(str,d['target_names'])))!=t:
        raise ValueError('target_names must be unique and mode-aware')
    if len(d['pair_names'])!=p or len(set(map(str,d['pair_names'])))!=p:
        raise ValueError('pair_names must be unique')
    e = np.asarray(d['single_effects']); v=np.asarray(d['objectives'])
    if e.ndim!=2 or e.shape[0]!=t or v.shape!=e.shape or e.shape[1]==0:
        raise ValueError('single_effects/objectives must share shape [T,G]')
    if not np.isfinite(e).all() or not np.isfinite(v).all() or not np.any(e):
        raise ValueError('Finite nonzero single-effect signal is required')
    if np.asarray(d['single_utility']).shape!=(t,t):
        raise ValueError('single_utility must be [T,T]')
    if not np.allclose(d['single_utility'],e@v.T,atol=1e-9,rtol=1e-7):
        raise ValueError('single_utility is inconsistent with single_effects and objectives')
    if np.asarray(d['Y']).shape!=(p,t) or np.asarray(d['RG']).shape!=(p,t):
        raise ValueError('Y and RG must be [P,T]')
    tr=np.asarray(d['train']); q=np.asarray(d['query'])
    for ids in (tr,q):
        if ids.ndim!=1 or ids.dtype.kind not in 'iu' or len(np.unique(ids))!=len(ids):
            raise ValueError('Split indices must be unique one-dimensional integers')
        if not len(ids) or np.min(ids)<0 or np.max(ids)>=p:
            raise ValueError('Empty or out-of-range split')
    if np.intersect1d(tr,q).size:
        raise ValueError('Outer train/query overlap')
    if not np.isfinite(d['Y'][tr]).all() or not np.isnan(d['Y'][q]).all():
        raise ValueError('Training labels must be finite; outer labels must be NaN')
    if not np.isfinite(d['RG'][np.r_[tr,q]]).all():
        raise ValueError('Required RG predictions must be finite')
    valid_ids=[]
    for c in [*cases,outer]:
        ci,cq=np.asarray(c['train']),np.asarray(c['query'])
        if np.intersect1d(ci,cq).size or not set(map(int,ci))<=set(map(int,tr)):
            raise ValueError('Invalid nested fit/query boundary')
        if len(c['train_blocks'])!=len(ci):
            raise ValueError('train_blocks shape mismatch')
        for side,ids in [('train',ci),('query',cq)]:
            for suffix in ['RG','graph','raw']:
                a=np.asarray(c[side+'_'+suffix])
                if a.shape!=(len(ids),t) or not np.isfinite(a).all():
                    raise ValueError(side+'_'+suffix+' must be finite [n,T]')
    for c in cases:
        if not set(map(int,c['query']))<=set(map(int,tr)):
            raise ValueError('Inner validation must stay in outer training')
        if set(map(int,c['train'])) | set(map(int,c['query'])) != set(map(int,tr)):
            raise ValueError('Each inner split must partition outer training')
        valid_ids.extend(map(int,c['query']))
    if sorted(valid_ids)!=sorted(map(int,tr)):
        raise ValueError('Inner validation blocks must cover training once')
    if set(map(int,outer['train']))!=set(map(int,tr)) or not np.array_equal(outer['query'],q):
        raise ValueError('Outer case identity mismatch')

def metrics(d,pred,blocks):
 us=[];cs=[]
 for va in blocks:
  for ai in range(d['X'].shape[1]):
   ids=va[d['X'][va,ai]>0]
   if len(ids)<2:continue
   # pair names alphabetical equivalent to partner alphabetical for fixed A
   order=sorted(ids,key=lambda i:d['pair_names'][i]);ps=pred[order,ai];truth=d['Y'][order,ai]
   pick=np.flatnonzero(ps>=ps.max()-1e-12)[0];us.append(truth[pick]);err=ps-truth;cs.append(np.mean((err-err.mean())**2))
 if not us:raise ValueError('No eligible within-block fixed-A choice groups')
 return dict(utility=float(np.mean(us)),center_mse=float(np.mean(cs)),groups=len(us))

def fixed_features(d):
 # e: singles by 4994 non-target RNA; target objectives follow the exact archived definition.
 e=d['single_effects'];names=d['target_names'];n=len(names);pairs=np.array([np.flatnonzero(x) for x in d['X']]);a,b=pairs.T
 # V was not copied into minimal inputs originally: recover via exact full-single directions and archived package.
 vp=None  # objective is an explicit training-only input
 v=d['objectives'];s=np.sqrt(np.mean(e*e,axis=0));s=np.maximum(s,np.sqrt(np.mean(e*e)));z=e/s
 add=e[a]+e[b]
 fields=[add]
 for alpha in (.02,.1,.5,1.):fields.append(add/(1+alpha*np.abs(add)))
 fields.extend([np.minimum(e[a],e[b]),np.maximum(e[a],e[b]),np.sign(e[a]*e[b])*np.minimum(abs(e[a]),abs(e[b])),
  s*(np.tanh(z[a]+z[b])-np.tanh(z[a])-np.tanh(z[b])),e[a]*e[b]/s])
 f=np.stack([y@v.T for y in fields],axis=-1)
 # Pair-conditioned scalar with explicit anchor/partner contribution preserves sign.
 su=d['single_utility']; u=d['X']@su
 own=np.diag(su)[None,:]*np.ones((d['X'].shape[0],1));partner=u-own
 extra=np.stack([partner,partner*own,partner*abs(partner),partner*own*own],axis=-1)
 return np.concatenate([f,extra],axis=-1)

def features(d,c,fix,side,mode):
 ids=c['train'] if side=='train' else c['query'];rg=c[side+'_RG'];fs=[rg[...,None]]
 if mode in ('simple','physical','context','graph'):
  fs.append((d['X'][ids]@d['single_utility'])[...,None])
 if mode in ('physical','context','graph'):fs.append(fix[ids,:,:10])
 if mode in ('context','graph'):
  fs.append(fix[ids,:,10:]);own=np.diag(d['single_utility'])[None,:];fs.extend([(rg*own)[...,None],(rg*abs(rg))[...,None]])
 if mode=='graph':
  gr=c[side+'_graph'];raw=c[side+'_raw'];su=d['single_utility']
  fs.extend([(gr@su)[...,None],(raw@su)[...,None]])
 return np.concatenate(fs,axis=-1)

def contrasts(d,c,X):
 ids=c['train'];groups=c['train_blocks'];XX=[];yy=[]
 for g in np.unique(groups):
  for a in range(d['X'].shape[1]):
   use=np.flatnonzero((groups==g)&(d['X'][ids,a]>0));m=len(use)
   if m<2:continue
   xx=X[use,a];y=d['Y'][ids[use],a]-c['train_RG'][use,a]
   XX.append((xx-xx.mean(0))/np.sqrt(m));yy.append((y-y.mean())/np.sqrt(m))
 if not XX:raise ValueError('No within-block contrasts: dataset/split lacks ranking support')
 return np.concatenate(XX),np.concatenate(yy)

def fitpredict(d,c,fix,mode,lam):
 f=features(d,c,fix,'train',mode);q=features(d,c,fix,'query',mode);X,y=contrasts(d,c,f)
 # Standardize using contrast training only; scale is unchanged for a new pool.
 scale=np.sqrt(np.mean(X*X,axis=0));scale=np.maximum(scale,1e-8);X=X/scale
 co=np.linalg.solve(X.T@X+lam*len(X)*np.eye(X.shape[1]),X.T@y)
 return (q/scale)@co,dict(coefficient=co.tolist(),scale=scale.tolist(),ncontrast_rows=len(y))

def old_branch(d,c):
 X=c['train_graph'];Y=d['Y'][c['train']]-c['train_RG'];gg=c['train_blocks'];rows=[]
 for lam in (.1,1.,10.):
  cv=np.zeros_like(Y)
  for g in np.unique(gg):
   fit=gg!=g;val=~fit;co=np.linalg.solve(X[fit].T@X[fit]+lam*np.eye(X.shape[1]),X[fit].T@Y[fit]);cv[val]=X[val]@co
  for frac in (0.,.25,.5,1.):
   losses=[]
   for g in np.unique(gg):
    for a in range(d['X'].shape[1]):
     sel=(gg==g)&(d['X'][c['train'],a]>0)
     if sel.sum()<2:continue
     err=frac*cv[sel,a]-Y[sel,a];losses.append(np.mean((err-err.mean())**2))
   rows.append(dict(lam=lam,fraction=frac,loss=float(np.mean(losses))))
 best=min(r['loss'] for r in rows);tol=64*np.finfo(float).eps*max(best,np.finfo(float).tiny)
 sel=min((r for r in rows if r['loss']<=best+tol),key=lambda r:(r['fraction'],-r['lam']))
 coef=np.linalg.solve(X.T@X+sel['lam']*np.eye(X.shape[1]),X.T@Y)
 return sel['fraction']*(c['query_graph']@coef),dict(selected=sel,coefficient=coef.tolist())

def train_ranker(d, cases, outer):
    """Fit the archived directional and graph ranking heads; no outer truth use."""
    validate_inputs(d,cases,outer)
    fix=fixed_features(d)
    blocks=[c['query'] for c in cases]
    rows=[dict(mode='R_G',lam=0.,fraction=0.,**metrics(d,d['RG'],blocks))]
    for mode in ('simple','physical','context','graph'):
        for lam in (.01,.1,1.,10.,100.):
            cor=np.zeros_like(d['RG'])
            for c in cases:
                cor[c['query']]=fitpredict(d,c,fix,mode,lam)[0]
            for fraction in (.25,.5,1.):
                rows.append(dict(mode=mode,lam=lam,fraction=fraction,
                                 **metrics(d,d['RG']+fraction*cor,blocks)))
    chosen=min(rows,key=lambda r:(-round(r['utility'],13),r['center_mse'],r['fraction']))
    if chosen['fraction']:
        direct,state=fitpredict(d,outer,fix,chosen['mode'],chosen['lam'])
    else:
        direct=np.zeros((len(d['query']),d['X'].shape[1]))
        state=dict(coefficient=[],scale=[],ncontrast_rows=0)
    graph,old=old_branch(d,outer)
    result=d['RG'][d['query']]+chosen['fraction']*direct+graph
    if not np.isfinite(result).all():
        raise FloatingPointError('Nonfinite final scores')
    record=dict(selected=chosen,candidates=rows,directional=state,legacy=old,
                training_ids=d['train'].tolist(),query_ids=d['query'].tolist())
    return result,record
