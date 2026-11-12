#!/usr/bin/env python3
"""Independent SVD recomputation; does not import the control implementation."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd,json
R=Path(__file__).resolve().parent
ia=np.load(R/'inputs/INNER_GROWTH_REPLAY.npz');oa=np.load(R/'inputs/CALIBRATION_OUTCOMES.npz')
d=pd.read_csv(R/'results/DECISIONS.csv').set_index('feature_id');s=pd.read_csv(R/'results/FOUR_ARM_SCORES.csv').set_index(['feature_id','policy'])
def orth(edges,n,v):
    x=np.zeros((len(edges),n));x[np.arange(len(edges))[:,None],edges]=1
    u,z,_=np.linalg.svd(x,full_matrices=False);rank=int((z>np.finfo(float).eps*max(x.shape)*z[0]).sum());u=u[:,:rank]
    return v-u@(u.T@v)
records=[]
for key in d.index:
    inner=lambda field:ia[key+'__'+field]
    n=int(inner('nodes'));base=inner('M0');pred=base if inner('selected_zero') else inner('regularized');a=b=0.
    for k in np.unique(inner('folds')):
        m=inner('folds')==k;r=orth(inner('edges')[m],n,inner('y')[m]-base[m]);z=orth(inner('edges')[m],n,pred[m]-base[m]);a+=float((r*z).sum());b+=float((z*z).sum())
    alpha=float(np.clip(a/b,0,1)) if b>1e-14 else 0.
    outer=lambda field:oa[key+'__'+field]
    y,base,orig,e=outer('truth'),outer('M0'),outer('original'),outer('edges');delta=orig-base;dp=orth(e,n,delta)
    policies={'baseline':base,'original':orig,'whole_shrink':base+alpha*delta,'partner_only':orig-(1-alpha)*dp}
    errors=[]
    for p,v in policies.items():
        er=y-v;pr=orth(e,n,er);f=float((er*er).sum());q=float((pr*pr).sum());stored=s.loc[(key,p)]
        errors.extend([abs(stored.full_sse-f),abs(stored.partner_sse-q)])
    records.append(dict(feature_id=key,alpha_error=abs(alpha-d.loc[key,'alpha']),max_sse_error=max(errors)))
df=pd.DataFrame(records);assert df.alpha_error.max()<1e-10 and df.max_sse_error.max()<1e-10
rna=pd.read_csv(R/'results/RNA_FOUR_ARM_SCORES.csv');w=rna.pivot(index=['record_id','pool_id','model'],columns='policy',values='partner_sse')
assert np.max(abs(w.whole_shrink-w.partner_only))<1e-10
df.to_csv(R/'results/INDEPENDENT_SVD_CHECKS.csv',index=False)
out=dict(status='PASS',functional_configurations=len(df),policy_comparisons=len(df)*4,max_alpha_error=float(df.alpha_error.max()),max_sse_error=float(df.max_sse_error.max()),RNA_partner_equivalence_rows=len(w),method='Independent SVD projector and direct quadratic-loss recomputation; no control implementation imported')
(R/'results/INDEPENDENT_CHECK.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
