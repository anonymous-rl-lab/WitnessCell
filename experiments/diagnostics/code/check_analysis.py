import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
from han_attribution import Project,features,fit_baselines
R=Path(__file__).resolve().parents[1];H=R/'inputs/han';O=R/'results';checks={}
for name in ['SCORES.csv','RANDOM_POOLED.csv','RANDOM_SUMMARY.csv','HAN_ORACLE_FOLDS.csv','HAN_REPEATABILITY.csv']:
 a=pd.read_csv(O/name);b=pd.read_csv(R/'replay'/name);assert a.shape==b.shape
 num=a.select_dtypes(include='number').columns;err=float(np.max(abs(a[num].to_numpy()-b[num].to_numpy())))
 checks['prepared_replay_'+name]=err;assert err<1e-10
pn='full_raw';inp=np.load(H/f'data/{pn}_inputs.npz');Y=np.load(H/f'data/{pn}_outcomes.npz')['outcomes'];genes=inp['genes'];pairs=inp['pairs'];idx={g:i for i,g in enumerate(genes)};edges=np.array([[idx[a],idx[b]] for a,b in pairs]);ids=json.loads((H/'data/splits'/pn/'20261108/0/identities.json').read_text());cal=np.array(ids['train']);query=np.array(ids['query']);B,Sat,Loc=features(inp['singles'],edges)
first,params,pc,pq=fit_baselines(edges,pairs,len(genes),B,Y,Sat,Loc,cal,query,20261108)
changed=Y.copy();changed[query]=np.random.default_rng(510).normal(size=changed[query].shape)*100
second,params2,_,_=fit_baselines(edges,pairs,len(genes),B,changed,Sat,Loc,cal,query,20261108)
assert params==params2
checks['query_outcome_replacement_max_prediction_change']=max(float(np.max(abs(first[k][1]-second[k][1]))) for k in first)
assert checks['query_outcome_replacement_max_prediction_change']==0
# Independent joint least squares for the selected saturation reference.
err=0.
for r in range(Y.shape[1]):
 x=np.column_stack([pc.X,Sat[cal,r,1]]);z=np.column_stack([pq.X,Sat[query,r,1]])
 coef=np.linalg.lstsq(x,Y[cal,r]-B[cal,r],rcond=None)[0];pred=B[query,r]+z@coef
 err=max(err,float(np.max(abs(pred-first['saturation'][1][:,r]))))
checks['joint_lstsq_vs_projected_saturation_max_abs']=err;assert err<1e-8
# Independent QR residualization for a learned direction and original two-block scalar.
F=np.load(H/'FROZEN_SF12_FEATURES.npz');ci=json.loads((H/'SOURCE_COEFFICIENTS.json').read_text());f=F[pn+'__HeLa'];theta=np.array(ci['HeLa']['theta']);D=np.einsum('nrf,f->nr',f,theta)
order=sorted(cal,key=lambda i:hashlib.sha256((str(20261108+1000)+'|'+'|'.join(sorted(pairs[i]))).encode()).hexdigest());A=C=0.
for block in np.array_split(order,2):
 pp=Project(edges[block],len(genes));u,sv,_=np.linalg.svd(pp.X,full_matrices=False);u=u[:,sv>sv.max()*1e-12];r=Y[block]-B[block];d=D[block];r-=u@(u.T@r);d-=u@(u.T@d);A+=float(np.sum(r*d));C+=float(np.sum(d*d))
alpha=float(np.clip(A/C,0,1));stored=json.loads((H/'results'/pn/'20261108/0/PARAMETERS.json').read_text())['alpha'][0]
checks['independent_svd_alpha_difference']=abs(alpha-stored);assert abs(alpha-stored)<1e-10
scores=pd.read_csv(O/'SCORES.csv');delta=(scores.full_sse-scores.base_full_sse)-(scores.partner_sse-scores.base_partner_sse)
checks['complete_partner_sse_identity_max_abs']=float(abs(delta).max());assert abs(delta).max()<1e-8
checks['status']='PASS';checks['scope']='Native and prepared new-control replay; query-outcome exclusion; independent baseline and projection arithmetic. No source-model retraining.'
(R/'replay/VALIDATION.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2))
