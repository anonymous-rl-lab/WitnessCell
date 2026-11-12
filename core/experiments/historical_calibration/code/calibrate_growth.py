"""Training-only bounded amplitude calibration; score only after persisting all decisions."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
R=Path(__file__).resolve().parents[1];O=R/'results';OLD=R.parent/'partner_gain_regression_20261108';S=R.parent/'gain_regression_sources';H=R.parent/'gain_regression_history'
d=pd.read_csv(O/'ENRICHED_INPUT.csv');z=d[d.model.isin(['rank1','rank1_selected'])&d.modality.eq('growth')].copy();assert len(z)==55
# No external outcome column is consulted to compute any calibration decision.
dec=z[['feature_id','panel','study','cell','cohort','query_ids','inner_reward','inner_cost','inner_partner_gain','inner_full_gain','inner_dims']].copy()
dec['alpha']=np.where(dec.inner_cost>1e-14,np.clip(dec.inner_reward/(2*dec.inner_cost),0,1),0)
dec['gate_alpha']=(dec.inner_partner_gain>0).astype(float)
dec['inner_predicted_gain']=dec.alpha*dec.inner_reward-dec.alpha**2*dec.inner_cost
dec.to_csv(O/'CALIBRATION_DECISIONS.csv',index=False)
protocol={'scope':'Retrospective development on recovered historical splits; not new preregistration or new independent confirmation','rule':'alpha=clip(inner_reward/(2*inner_cost),0,1); zero for zero correction; no threshold search','objective':'Minimize pooled foldwise endpoint-projected inner SSE along the archived rank-one correction','reference':'Original full-response-selected rank-one model, not best outer predictor','data':'All 55 eligible functional training configurations','decision_sha256':hashlib.sha256((O/'CALIBRATION_DECISIONS.csv').read_bytes()).hexdigest(),'outer_outcomes_used_for_alpha':False,'fitting_biological_predictor':False}
(O/'CALIBRATION_PROTOCOL.json').write_text(json.dumps(protocol,indent=2))

def X(e,n):
 x=np.zeros((len(e),n));x[np.arange(len(e))[:,None],e]=1;return x
def P(e,n,y):
 x=X(e,n);return y-x@np.linalg.lstsq(x,y,rcond=None)[0]
def predict(mods,singles,e):
 scale=np.std(singles,0);scale[scale==0]=1;i,j=e.T;cols=[]
 for k,m in enumerate(mods):
  a,ep,db,dc=[np.asarray(m[s]) for s in ['a','e','d_beta','d_c']]
  val=.5*(a[i]+a[j]+(m['beta']+db[i])*a[j]+(m['beta']+db[j])*a[i]+(m['c']+dc[i])*a[j]**2+(m['c']+dc[j])*a[i]**2)+m['e0']+ep[i]+ep[j];cols.append(val*scale[k])
 return np.column_stack(cols)
meta={a['label']:a for a in json.loads((OLD/'inputs/metadata.json').read_text())};rows=[];checks=[];reps=[];arrays={}
for row in dec.itertuples():
 key=row.feature_id
 if key in meta:
  a=np.load(OLD/'inputs'/meta[key]['outcome_file']);tr=np.load(OLD/'inputs'/meta[key]['train_file']);e=a['edges'];n=len(tr['genes']);y=a['y'];base=a['M0'];orig=a['rank1'];group=key.rsplit('_fold',1)[0]
 else:
  import re
  f,seed,arm=re.fullmatch(r'K562_layout_f(\d)_s(\d)_(.*)',key).groups();md=H/'bridge_restored/WitnessCell_NM_NCS_Bridge_20261031/growth_run/training_v2'/f'fold_{f}_seed_{seed}'/'results'/arm
  source=S/'horlbeck/WitnessCell_Horlbeck_TwoStage_Repro_20261018/stages/K562/cases'/f'outer{f}';tr=np.load(source/'inputs/training.npz');ev=np.load(source/'inputs/evaluation.npz');e=ev['edges'];y=ev['y'];n=len(tr['genes']);base=predict(json.loads((md/'M0_models.json').read_text()),tr['singles'],e);sel=json.loads((md/'selection.json').read_text())['rank1']
  if sel['selected']=='G0':orig=base.copy()
  else:
   m=json.loads((md/'rank1_model.json').read_text());u=np.asarray(m['u']);corr=m['lambda']*u[e[:,0]]*u[e[:,1]];orig=base+corr[:,None]
  group='K562_layout'
 r=P(e,n,y-base);delta=P(e,n,orig-base);den=np.sum(r*r);reward=2*np.sum(r*delta);cost=np.sum(delta*delta)
 old=z[z.feature_id==key].iloc[0];original_gain=100*(reward-cost)/den;checks.append(abs(original_gain-old.gain_pp));assert abs(original_gain-old.gain_pp)<1e-7,(key,original_gain,old.gain_pp)
 arrays[key+'__truth']=y;arrays[key+'__M0']=base;arrays[key+'__original']=orig;arrays[key+'__edges']=e;arrays[key+'__nodes']=np.array(n)
 for name,alpha in [('original',1.),('zero',0.),('inner_gate',row.gate_alpha),('inner_shrink',row.alpha)]:
  pred=base+alpha*(orig-base);p_sse=float(np.sum(P(e,n,y-pred)**2));f_sse=float(np.sum((y-pred)**2));f_base=float(np.sum((y-base)**2));gain=100*(1-p_sse/den)
  assert abs(gain-100*(alpha*reward-alpha**2*cost)/den)<1e-8
  rows.append(dict(feature_id=key,panel=row.panel,study=row.study,cell=row.cell,group=group,cohort=row.cohort,policy=name,alpha=alpha,partner_base_sse=den,partner_sse=p_sse,partner_gain_pp=gain,full_base_sse=f_base,full_sse=f_sse,full_gain_pp=100*(1-f_sse/f_base),outer_reward_pp=100*reward/den,outer_cost_pp=100*cost/den,query_pairs=len(e)))
  for j in range(y.shape[1]):
   b=np.sum(r[:,j]**2);loss=np.sum(P(e,n,y[:,j]-pred[:,j])**2);reps.append(dict(feature_id=key,group=group,policy=name,repeat=j+1,partner_base_sse=b,partner_sse=loss,gain_pp=100*(1-loss/b)))
np.savez_compressed(R/'inputs/CALIBRATION_OUTCOMES.npz',**arrays)
a=pd.DataFrame(rows);a.to_csv(O/'CALIBRATION_SCORES.csv',index=False);pd.DataFrame(reps).to_csv(O/'CALIBRATION_REPEATS.csv',index=False)
g=a.groupby(['group','policy']).agg(base=('partner_base_sse','sum'),loss=('partner_sse','sum'),full_base=('full_base_sse','sum'),full_loss=('full_sse','sum'),configs=('feature_id','nunique')).reset_index();g['partner_gain_pp']=100*(1-g.loss/g.base);g['full_gain_pp']=100*(1-g.full_loss/g.full_base);g.to_csv(O/'CALIBRATION_BY_CONDITION.csv',index=False)
# Layout queries are repeated across seed/arm; pooled numbers are descriptive configuration averages, not new sample size.
pan=a.groupby(['panel','policy']).partner_gain_pp.mean().reset_index();pan.to_csv(O/'CALIBRATION_BY_PANEL.csv',index=False)
(O/'CALIBRATION_REPLAY_CHECKS.json').write_text(json.dumps({'original_replay_max_gain_difference_pp':max(checks),'configurations':len(dec),'status':'PASS','decision_sha256_after_scoring':hashlib.sha256((O/'CALIBRATION_DECISIONS.csv').read_bytes()).hexdigest()},indent=2))
print(g.to_string(index=False))
