import numpy as np,pandas as pd
def weights(z):
 # Alias equivalence includes the archived training features and score; retained in ledger.
 w=1/z.groupby('score_key').score_key.transform('size').to_numpy(float)
 w/=pd.Series(w,index=z.index).groupby(z.feature_id).transform('sum').to_numpy()
 w/=pd.Series(w,index=z.index).groupby(z.panel).transform('sum').to_numpy()
 return w/w.sum()
def split(z,group,value):
 te=(z[group]==value).to_numpy();tr=~te
 if group=='cohort':
  ids=set(';'.join(z.loc[te,'query_ids'].unique()).split(';'));tr &= ~z.query_ids.map(lambda v:bool(ids.intersection(v.split(';')))).to_numpy()
 return tr,te

def fit(tr,te,cols,cats):
 w=weights(tr);y=tr.gain_pp.to_numpy();ym=w@y;xx=[];tt=[];names=[];scalers={}
 for c in cols:
  x=tr[c].to_numpy(float);t=te[c].to_numpy(float);missing=~np.isfinite(x);mt=~np.isfinite(t);good=~missing
  mu=float(w[good]@x[good]/w[good].sum()) if good.any() else 0.;x=np.where(good,x,mu);t=np.where(mt,mu,t);sd=float(np.sqrt(w@((x-mu)**2)));sd=sd if sd>1e-8 else 1.
  xx.append((x-mu)/sd);tt.append((t-mu)/sd);names.append(c);scalers[c]={'mean':mu,'sd':sd}
  if missing.any():xx.append(missing.astype(float));tt.append(mt.astype(float));names.append(c+'_missing')
 for c in cats:
  for lv in sorted(tr[c].fillna('missing').unique()):
   xx.append((tr[c].fillna('missing')==lv).to_numpy(float));tt.append((te[c].fillna('missing')==lv).to_numpy(float));names.append(c+'='+str(lv))
 if not xx:return np.full(len(te),ym),{'intercept':float(ym)}
 x=np.column_stack(xx);t=np.column_stack(tt);xm=w@x;x-=xm;t-=xm;ww=w*tr.feature_id.nunique();beta=np.linalg.solve(x.T@(ww[:,None]*x)+np.eye(len(names)),x.T@(ww*(y-ym)))
 return ym+t@beta,dict(mean_y=float(ym),centers=dict(zip(names,xm.tolist())),coefficients=dict(zip(names,beta.tolist())),scalers=scalers,alpha=1.)
def metric(z,p):
 w=weights(z);y=z.gain_pp.to_numpy();err=p-y;return dict(rmse_pp=float(np.sqrt(w@(err*err))),mae_pp=float(w@abs(err)),R2=float(1-(w@(err*err))/max(w@((y-w@y)**2),1e-30)))
basecats=['model_family','modality'];configs={
 'mean':([],[]),
 'data_only':(['q_raw','c_raw','log_degree'],basecats),
 'inner_partner':(['inner_partner_gain'],[]),
 'inner_full':(['inner_full_gain'],[]),
 'data_and_inner':(['q_raw','c_raw','log_degree','inner_partner_gain','inner_full_gain'],basecats),
 'validation_and_capacity':(['inner_partner_gain','inner_full_gain','log_inner_dims','amplitude','log_regularization','log_effective_df','model_rank'],basecats),
 'data_validation_capacity':(['q_raw','c_raw','log_degree','inner_partner_gain','inner_full_gain','log_inner_dims','amplitude','log_regularization','log_effective_df','model_rank'],basecats),
 'validation_reward_cost':(['inner_reward','inner_cost','inner_full_gain','log_degree'],basecats)}
