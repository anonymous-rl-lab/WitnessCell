"""Retrospective fixed-direction oracles and separately labelled repeat correlations."""
from pathlib import Path
import sys,json, argparse
import numpy as np,pandas as pd
from han_attribution import Project
ap=argparse.ArgumentParser();ap.add_argument('--control',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();R=a.control;O=a.out;O.mkdir(exist_ok=True,parents=True)
meta=pd.read_csv(R/'inputs/CONFIGURATIONS.csv');inn=np.load(R/'inputs/INNER_GROWTH_REPLAY.npz');out=np.load(R/'inputs/CALIBRATION_OUTCOMES.npz');dec=pd.read_csv(R/'results/DECISIONS.csv').set_index('feature_id');rows=[];qrows=[]
def corr(v):return float(np.corrcoef(v[:,0],v[:,1])[0,1])
for row in meta.itertuples():
 k=row.feature_id;edges=out[k+'__edges'];n=int(out[k+'__nodes']);P=Project(edges,n);y=out[k+'__truth'];b=out[k+'__M0'];r=P.p(y-b);d=P.p(out[k+'__original']-b);A=float(np.sum(r*d));B=float(np.sum(d*d));E=float(np.sum(r*r));alpha=float(dec.loc[k,'alpha']);u=A/B if B>1e-14 else 0;an=max(0,u);ab=float(np.clip(u,0,1));gain=lambda t:100*(2*t*A-t*t*B)/E
 rows.append(dict(feature_id=k,correction_zero=B<=1e-14,actual_alpha=alpha,unconstrained_alpha=u,nonnegative_alpha=an,bounded_alpha=ab,actual_gain_pct=gain(alpha),nonnegative_oracle_gain_pct=gain(an),bounded_oracle_gain_pct=gain(ab),original_gain_pct=gain(1),partner_E=E,partner_A=A,partner_B=B))
 if k.startswith('K562_layout'):continue
 Pi=Project(inn[k+'__edges'],int(inn[k+'__nodes']));q=corr(Pi.p(inn[k+'__y']));qo=corr(P.p(y));res=corr(P.p(y-b))
 # Retain both metrics and the historical field without assuming its numerical convention.
 v=Pi.p(inn[k+'__y']);cos=float(np.sum(v[:,0]*v[:,1])/np.linalg.norm(v[:,0])/np.linalg.norm(v[:,1]))
 qrows.append(dict(feature_id=k,training_pairs=len(v),query_pairs=len(y),training_repeat_pearson=q,training_repeat_cosine=cos,archived_q_raw=row.q_raw,query_repeat_pearson=qo,query_residual_repeat_pearson=res,original_partner_gain_pct=gain(1),controlled_partner_gain_pct=gain(alpha)))
pd.DataFrame(rows).to_csv(O/'FUNCTIONAL_ORACLES.csv',index=False);pd.DataFrame(qrows).to_csv(O/'FUNCTIONAL_REPEATABILITY.csv',index=False)
print(pd.DataFrame(qrows).to_string(index=False));print(pd.DataFrame(rows).query("feature_id in ['HeLa_high','A549_high','293T_high','Jurkat_high']").to_string(index=False))
