"""Compare identical eligible rows; whole-query purging and source-held model selection."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import numpy as np,pandas as pd,json
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'results';d=pd.read_csv(O/'ENRICHED_INPUT.csv');d['zero_to_reference']=d.zero_to_reference.fillna(False).astype(bool);d=d.reset_index(drop=True)

from common import weights, split, fit, metric, configs
perfs=[];rows=[];coefs={};splits=[]
for name,(cols,cats) in configs.items():
 for scheme,group in [('fit',None),('hold_query','cohort'),('hold_panel','panel'),('hold_study','study')]:
  p=np.empty(len(d))
  if group:
   for val in d[group].unique():
    tr,te=split(d,group,val);p[te]=fit(d[tr],d[te],cols,cats)[0]
    if name=='mean':splits.append(dict(scheme=scheme,hold=str(val),train=int(tr.sum()),test=int(te.sum()),purged=len(d)-int(tr.sum())-int(te.sum())))
  else:p,coefs[name]=fit(d,d,cols,cats)
  # A mathematically zero correction has known zero gain; identical rule for every comparator.
  p[d.zero_to_reference]=0.
  for scope,ix in [('all',np.ones(len(d),bool)),('active',~d.zero_to_reference.to_numpy())]:
   perfs.append(dict(regression=name,scheme=scheme,scope=scope,n=int(ix.sum()),**metric(d[ix],p[ix])))
  for i in range(len(d)):rows.append(dict(row=i,regression=name,scheme=scheme,predicted=p[i],observed=d.gain_pp.iloc[i]))
pd.DataFrame(perfs).to_csv(O/'REGRESSION_PERFORMANCE.csv',index=False);pr=pd.DataFrame(rows);pr.to_csv(O/'HELD_PREDICTIONS.csv',index=False);(O/'COEFFICIENTS.json').write_text(json.dumps(coefs,indent=2));pd.DataFrame(splits).to_csv(O/'SPLITS.csv',index=False)
d['analysis_weight']=weights(d);d.to_csv(O/'WEIGHTED_INPUT.csv',index=False)
# Baselines not fitted: the historical inner score as an uncalibrated forecast.
idperf=[]
for col in ['inner_partner_gain','inner_full_gain']:
 pp=d[col].to_numpy();pp[d.zero_to_reference]=0.;idperf.append(dict(forecast=col,**metric(d,pp)))
pd.DataFrame(idperf).to_csv(O/'UNCALIBRATED_FORECASTS.csv',index=False)
# Candidate selection inside matched pools, zero alternative included for every task.
policies=[];default_order=['rank1','rank1_selected','R_original','R','core_features_family_selected_P','input_probe_input_selected_P','network_v4_full_selected_delta','shared12_selected']
for pool,z0 in d.groupby('pool_id',sort=True):
 # Keep identical source aliases deterministic. Ties are broken by model name, never outcome.
 z=z0.sort_values(['model','candidate_key']).copy();assert z.baseline.nunique()==1,pool
 den=z.base_partner_sse.to_numpy();assert np.max(abs(den/den[0]-1))<1e-8,pool
 ref=next((z[z.model==m].iloc[0] for m in default_order if (z.model==m).any()),None)
 def record(policy,choice,forecast=np.nan):
  policies.append(dict(pool_id=pool,feature_id=z.feature_id.iloc[0],panel=z.panel.iloc[0],study=z.study.iloc[0],cell=z.cell.iloc[0],cohort=z.cohort.iloc[0],policy=policy,chosen_model='M0' if choice is None else choice['model'],candidate_id='' if choice is None else choice['candidate_key'],predicted_gain=float(forecast),observed_gain=0. if choice is None else float(choice['gain_pp']),candidate_count=len(z),oracle_gain=max(0.,float(z.gain_pp.max()))))
 def select(policy,values):
  k=int(np.argmax(values));record(policy,None if values[k]<=1e-10 else z.iloc[k],values[k])
 record('zero',None,0.)
 if ref is not None:record('historical_default',ref)
 select('inner_partner',z.inner_partner_gain.to_numpy());select('inner_full',z.inner_full_gain.to_numpy())
 # Fixed dual guard, no threshold tuning on outer outcomes.
 vals=np.where(z.inner_full_gain>=0,z.inner_partner_gain,-np.inf);select('inner_both_nonnegative',vals)
 for scheme in ['hold_query','hold_panel','hold_study']:
  for method in ['data_only','inner_partner','data_and_inner','validation_and_capacity','data_validation_capacity']:
   p=pr[(pr.regression==method)&(pr.scheme==scheme)].set_index('row');values=p.loc[z.index,'predicted'].to_numpy();select(scheme+'__'+method,values)
po=pd.DataFrame(policies);po['regret_pp']=po.oracle_gain-po.observed_gain;po.to_csv(O/'MODEL_SELECTION_RECORDS.csv',index=False)
# Average pools within feature configurations, then configurations within panel, then panels.
a=po.groupby(['policy','panel','feature_id']).agg(gain=('observed_gain','mean'),regret=('regret_pp','mean'),selected=('chosen_model',lambda s:float((s!='M0').mean())),negative=('observed_gain',lambda x:float((x< -1e-8).mean())),pools=('pool_id','size')).reset_index()
b=a.groupby(['policy','panel']).agg(gain=('gain','mean'),regret=('regret','mean'),selected=('selected','mean'),negative=('negative','mean'),configs=('feature_id','size'),pools=('pools','sum')).reset_index();b.to_csv(O/'SELECTION_BY_PANEL.csv',index=False)
c=b.groupby('policy').agg(mean_gain_pp=('gain','mean'),mean_regret_pp=('regret','mean'),selection_rate=('selected','mean'),negative_rate=('negative','mean'),panels=('panel','size')).reset_index();c.to_csv(O/'SELECTION_SUMMARY.csv',index=False)
# Compare within-configuration slopes to separate model matching from panel shifts.
ass=[]
for col in ['inner_partner_gain','inner_full_gain','amplitude','log_regularization','log_effective_df']:
 z=d[np.isfinite(d[col])].copy();x=z[col]-z.groupby('feature_id')[col].transform('mean');y=z.gain_pp-z.groupby('feature_id').gain_pp.transform('mean');w=weights(z);cov=w@(x*y);den=np.sqrt((w@(x*x))*(w@(y*y)))
 ass.append(dict(factor=col,rows=len(z),within_configuration_r=float(cov/den) if den>0 else np.nan))
pd.DataFrame(ass).to_csv(O/'WITHIN_CONFIGURATION_ASSOCIATIONS.csv',index=False)
print('n',len(d),'training configs',d.feature_id.nunique(),'panels',d.panel.nunique(),'studies',d.study.nunique(),'pools',po.pool_id.nunique());print(pd.DataFrame(perfs).query("scope=='all'").to_string(index=False));print(c.to_string(index=False))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True);colors=plt.get_cmap('tab10')
ax=axs[0,0]
for i,(panel,z) in enumerate(d.groupby('panel')):
 z=z.drop_duplicates('score_key');ax.scatter(z.inner_partner_gain,z.gain_pp,s=15,alpha=.55,color=colors(i),label=panel)
ax.axhline(0,color='.7',lw=.7);ax.axvline(0,color='.7',lw=.7);ax.set_xscale('symlog',linthresh=1);ax.set_yscale('symlog',linthresh=1);ax.set(title='A  Archived inner versus outer partner gain',xlabel='Inner gain (%) — symlog',ylabel='Outer gain (%) — symlog');ax.legend(fontsize=7,ncol=2,frameon=False)
ax=axs[0,1];pf=pd.DataFrame(perfs);models=['data_only','inner_partner','data_and_inner','validation_and_capacity'];labels=['Data','Inner','Data + inner','Inner + capacity']
for off,scheme,col in [(-.23,'fit','#cbd5e1'),(0,'hold_query','#0f766e'),(.23,'hold_study','#d97706')]:
 z=pf[(pf.scope=='all')&(pf.scheme==scheme)].set_index('regression');ax.bar(np.arange(4)+off,[z.loc[m,'rmse_pp'] for m in models],.22,label=scheme.replace('_',' '),color=col)
ax.set_xticks(range(4),labels,rotation=12);ax.set(title='B  Gain forecast error',ylabel='RMSE (percentage points)');ax.legend(fontsize=7,frameon=False)
ax=axs[1,0];pols=['historical_default','inner_partner','inner_both_nonnegative','hold_study__data_and_inner'];names=['Historical default','Inner partner','Inner dual guard','Study-held regression'];ss=c.set_index('policy');ax.barh(names,[ss.loc[k,'mean_gain_pp'] for k in pols],color=['#94a3b8','#14b8a6','#0f766e','#d97706']);ax.axvline(0,color='.6',lw=.7);ax.set(title='C  Candidate selection, panel-balanced',xlabel='Realized partner gain (%)')
ax=axs[1,1]
for off,pol,col in [(-.18,'historical_default','#94a3b8'),(.18,'inner_both_nonnegative','#0f766e')]:
 z=b[b.policy==pol].set_index('panel');panels=sorted(d.panel.unique());ax.bar(np.arange(len(panels))+off,[z.loc[k,'gain'] if k in z.index else np.nan for k in panels],.35,color=col,label=pol.replace('_',' '))
ax.axhline(0,color='.6',lw=.7);ax.set_xticks(range(len(panels)),panels,rotation=35);ax.set(title='D  Selection gains by data panel',ylabel='Realized partner gain (%)');ax.legend(fontsize=7,frameon=False)
fig.savefig(O/'Model_Matched_Gain.png',dpi=210);fig.savefig(O/'Model_Matched_Gain.pdf');plt.close(fig)
