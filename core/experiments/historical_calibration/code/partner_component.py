"""Secondary fixed-query-batch decomposition; training-only alpha, no label use in predictor."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd,json
R=Path(__file__).resolve().parents[1];O=R/'results'
dec=pd.read_csv(O/'CALIBRATION_DECISIONS.csv');arc=np.load(R/'inputs/CALIBRATION_OUTCOMES.npz');old=pd.read_csv(O/'CALIBRATION_SCORES.csv')
rows=[];checks=[]
for r in dec.itertuples():
 key=r.feature_id;e=arc[key+'__edges'];n=int(arc[key+'__nodes']);base=arc[key+'__M0'];orig=arc[key+'__original']
 x=np.zeros((len(e),n));x[np.arange(len(e))[:,None],e]=1
 delta=orig-base;endpoint=x@np.linalg.lstsq(x,delta,rcond=None)[0];partner=delta-endpoint
 # Query identities enter the projection. Query outcomes do NOT enter this predictor.
 pred=base+endpoint+r.alpha*partner
 y=arc[key+'__truth'];res=y-pred;pres=res-x@np.linalg.lstsq(x,res,rcond=None)[0]
 ref=old[(old.feature_id==key)&(old.policy=='inner_shrink')].iloc[0]
 original=old[(old.feature_id==key)&(old.policy=='original')].iloc[0]
 p_sse=float(np.sum(pres**2));f_sse=float(np.sum(res**2));pg=100*(1-p_sse/ref.partner_base_sse)
 same=abs(pg-ref.partner_gain_pp);orth=float(np.max(abs(x.T@partner)))
 identity=abs((original.full_sse-f_sse)-(original.partner_sse-p_sse))
 checks.append(dict(feature_id=key,partner_gain_difference_vs_whole_shrink_pp=same,orthogonality_error=orth,full_partner_improvement_identity_error=identity))
 assert same<1e-7 and orth<1e-8 and identity<1e-8
 rows.append(dict(feature_id=key,panel=r.panel,group=ref.group,alpha=r.alpha,partner_base_sse=ref.partner_base_sse,partner_sse=p_sse,full_base_sse=ref.full_base_sse,full_sse=f_sse,partner_gain_pp=pg,full_gain_pp=100*(1-f_sse/ref.full_base_sse)))
a=pd.DataFrame(rows);a.to_csv(O/'PARTNER_COMPONENT_SCORES.csv',index=False)
g=a.groupby('group').agg(base=('partner_base_sse','sum'),loss=('partner_sse','sum'),full_base=('full_base_sse','sum'),full_loss=('full_sse','sum'),configs=('feature_id','size')).reset_index();g['partner_gain_pp']=100*(1-g.loss/g.base);g['full_gain_pp']=100*(1-g.full_loss/g.full_base);g.to_csv(O/'PARTNER_COMPONENT_BY_CONDITION.csv',index=False)
(O/'PARTNER_COMPONENT_CHECKS.json').write_text(json.dumps(checks,indent=2))
(O/'PARTNER_COMPONENT_SCOPE.json').write_text(json.dumps({'status':'Retrospective secondary structural transformation, added after observing whole-correction calibration','rule':'M0 + (I-P_query) delta + alpha P_query delta','alpha':'Identical archived training-only bounded calibration; no re-optimization','query_labels_used_to_predict':False,'query_identities_used':True,'batch_dependent':True,'scope':'Fixed query panels. Not an inductive, query-batch-invariant deployed predictor. Single-query contrast dimension can be zero.','scientific_meaning':'Retains original endpoint contribution while altering only the partner projection; full SSE improvement equals partner SSE improvement by orthogonality.'},indent=2))
print(g.to_string(index=False))
