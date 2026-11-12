"""Independent SVD scoring and archived decision checks; no original fitting code used."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
from common import weights,split,metric
R=Path(__file__).resolve().parents[1];O=R/'results';checks=[]
def ck(name,ok,details=None):
 checks.append({'check':name,'passed':bool(ok),'details':details});assert ok,(name,details)
def P(edges,n,y):
 x=np.zeros((len(edges),n));x[np.arange(len(edges))[:,None],edges]=1
 u,s,_=np.linalg.svd(x,full_matrices=False);keep=s>max(x.shape)*np.finfo(float).eps*s[0];u=u[:,keep]
 return y-u@(u.T@y)
d=pd.read_csv(O/'ENRICHED_INPUT.csv');d['zero_to_reference']=d.zero_to_reference.fillna(False).astype(bool)
ck('eligible_rows_configs_panels_studies', (len(d),d.feature_id.nunique(),d.panel.nunique(),d.study.nunique())==(788,64,8,5))
err=np.max(abs(d.gain_pp-100*(1-d.partner_sse/d.base_partner_sse)));ck('all_outer_gain_definitions',err<1e-7,float(err))
ck('zero_corrections_have_zero_gain',np.max(abs(d.loc[d.zero_to_reference,'gain_pp']))<1e-7)
ck('twelve_coefficients_not_matrix_rank',d.loc[d.model=='shared12_selected','model_rank'].isna().all())
ck('no_missing_inner_scores',d[['inner_partner_gain','inner_full_gain']].notna().all().all())
w=weights(d);ck('panel_weights_equal',np.max(abs(pd.Series(w).groupby(d.panel).sum()-1/8))<1e-12)
for cohort in d.cohort.unique():
 tr,te=split(d,'cohort',cohort);a=set(';'.join(d.loc[tr,'query_ids']).split(';'));b=set(';'.join(d.loc[te,'query_ids']).split(';'));ck('purged_query_'+cohort,not(a&b))
for study in d.study.unique():
 tr,te=split(d,'study',study);ck('study_exclusion_'+study,not set(d.loc[tr,'study'])&set(d.loc[te,'study']))
arc=np.load(R/'inputs/INNER_GROWTH_REPLAY.npz');dec=pd.read_csv(O/'CALIBRATION_DECISIONS.csv');maxerrs=[]
for r in dec.itertuples():
 key=r.feature_id;y=arc[key+'__y'];base=arc[key+'__M0'];pred=arc[key+'__regularized'];e=arc[key+'__edges'];n=int(arc[key+'__nodes']);fold=arc[key+'__folds']
 if bool(arc[key+'__selected_zero']):pred=base.copy()
 den=reward=cost=0.
 for f in sorted(set(fold)):
  ix=fold==f;res=P(e[ix],n,y[ix]-base[ix]);delta=P(e[ix],n,pred[ix]-base[ix]);den+=np.sum(res**2);reward+=2*np.sum(res*delta);cost+=np.sum(delta**2)
 vals=np.array([100*reward/den,100*cost/den,100*(reward-cost)/den,100*(1-np.sum((y-pred)**2)/np.sum((y-base)**2))]);saved=np.array([r.inner_reward,r.inner_cost,r.inner_partner_gain,r.inner_full_gain]);err=float(np.max(abs(vals-saved)));maxerrs.append(err)
 ck('inner_training_SVD_'+key,err<1e-7,err)
 alpha=np.clip(reward/(2*cost),0,1) if cost>1e-14*den/100 else 0.
 ck('training_only_alpha_'+key,abs(alpha-r.alpha)<1e-10)
 ck('inner_quadratic_minimizer_'+key,alpha*reward-alpha*alpha*cost>=max(0.,reward-cost)-1e-10)
proto=json.loads((O/'CALIBRATION_PROTOCOL.json').read_text());ck('decision_hash_unchanged',hashlib.sha256((O/'CALIBRATION_DECISIONS.csv').read_bytes()).hexdigest()==proto['decision_sha256'])
arc=np.load(R/'inputs/CALIBRATION_OUTCOMES.npz');sc=pd.read_csv(O/'CALIBRATION_SCORES.csv');comp=pd.read_csv(O/'PARTNER_COMPONENT_SCORES.csv');outerr=[]
for r in dec.itertuples():
 key=r.feature_id;y=arc[key+'__truth'];base=arc[key+'__M0'];orig=arc[key+'__original'];e=arc[key+'__edges'];n=int(arc[key+'__nodes'])
 for policy,alpha in [('zero',0.),('original',1.),('inner_gate',r.gate_alpha),('inner_shrink',r.alpha)]:
  p=base+alpha*(orig-base);loss=np.sum(P(e,n,y-p)**2);full=np.sum((y-p)**2);ref=sc[(sc.feature_id==key)&(sc.policy==policy)].iloc[0]
  err=float(max(abs(loss-ref.partner_sse),abs(full-ref.full_sse)));outerr.append(err);ck('outer_SVD_'+key+'_'+policy,err<1e-8,err)
 delta=orig-base;dp=P(e,n,delta);p=orig-(1-r.alpha)*dp;ref=comp[comp.feature_id==key].iloc[0]
 ck('partner_only_SVD_'+key,max(abs(np.sum((y-p)**2)-ref.full_sse),abs(np.sum(P(e,n,y-p)**2)-ref.partner_sse))<1e-8)
 # Predictor invariance to labels follows from using only base, original, edges and alpha above.
 ck('endpoint_component_preserved_'+key,np.max(abs((p-base)-P(e,n,p-base)-(delta-dp)))<1e-9)
perf=pd.read_csv(O/'REGRESSION_PERFORMANCE.csv');pr=pd.read_csv(O/'HELD_PREDICTIONS.csv')
for (model,scheme),v in pr.groupby(['regression','scheme']):
 vals=v.sort_values('row').predicted.to_numpy()
 for scope,ix in [('all',np.ones(len(d),bool)),('active',~d.zero_to_reference.to_numpy())]:
  m=metric(d[ix],vals[ix]);r=perf[(perf.regression==model)&(perf.scheme==scheme)&(perf.scope==scope)].iloc[0];ck('regression_metric_'+model+'_'+scheme+'_'+scope,max(abs(m[c]-r[c]) for c in m)<1e-8)
po=pd.read_csv(O/'MODEL_SELECTION_RECORDS.csv')
for pool,z in d.groupby('pool_id'):
 z=z.sort_values(['model','candidate_key']);ck('matched_baseline_'+hashlib.sha256(pool.encode()).hexdigest()[:10],z.baseline.nunique()==1 and np.max(abs(z.base_partner_sse/z.base_partner_sse.iloc[0]-1))<1e-8)
 for policy,col in [('inner_partner','inner_partner_gain'),('inner_full','inner_full_gain'),('inner_both_nonnegative','inner_partner_gain')]:
  val=z[col].to_numpy().copy()
  if policy=='inner_both_nonnegative':val[z.inner_full_gain<0]=-np.inf
  ix=int(np.argmax(val));chosen='M0' if val[ix]<=1e-10 else z.iloc[ix].model;ref=po[(po.pool_id==pool)&(po.policy==policy)].iloc[0]
  ck('selection_'+policy+'_'+hashlib.sha256(pool.encode()).hexdigest()[:10],chosen==ref.chosen_model)
coverage=pd.read_csv(O/'COVERAGE.csv');ck('all_historical_records_retained_in_coverage',len(coverage)==1280)
out={'status':'PASS','checks':len(checks),'inner_svd_max_abs_difference_pp':max(maxerrs),'outer_svd_max_sse_difference':max(outerr),'scope':'Independent SVD rescore of all 55 functional inner/outer configurations, exact archived regression metric and selection verification; no RNA predictor retraining or independent biological confirmation.','items':checks}
(O/'VERIFICATION.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='items'},indent=2))
