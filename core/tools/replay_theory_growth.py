"""Recompute controlled-system and functional scores from original primitive records."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
from scipy.stats import t
R=Path(__file__).resolve().parents[1];ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out.resolve()
if out.exists() or out.is_relative_to(R):ap.error('Use a fresh output outside the release.')
out.mkdir(parents=True);scores=R/'paper/evidence/source_data/formal_bridge/scores'
# Exact finite-family population assertions, independent of fitted coefficients.
letters='ABCDEF';pilot=[()]+[(x,) for x in 'XABCDEF'];low=[('A','D'),('B','D'),('B','E'),('C','F')];high=[('A','B'),('B','C'),('D','E'),('D','F')];query=[('A','C'),('A','F'),('C','E'),('E','F')]
def state(intervention,x):
 x=0 if 'X' in intervention else x;z={k:0 if k in intervention else x for k in letters};A,B,C=[z[k] for k in 'ABC'];return x,z,[A or (B and C),B or (A and C),C or (A and B)]
for intervention in pilot+low:
 for x in [0,1]:
  u,z,g=state(intervention,x);assert g==[u]*3
assert sorted([sum(k in pair for pair in low) for k in letters])==sorted([sum(k in pair for pair in high) for k in letters])
weights=np.array([1,-1,-1,1]);targets=[]
for world in range(4):
 m=[]
 for pair in query:
  values=[]
  for x in [0,1]:
   u,z,g=state(pair,x);values.append(.7+np.array([.25,-.15,.4,-.2,.3,-.1])@np.array([z[k] for k in letters])+(u if world==0 else g[world-1]))
  m.append(np.mean(values))
 targets.append(float(weights@np.array(m)))
np.testing.assert_allclose(targets,[0,-.5,0,-.5],atol=1e-15);floor=float(np.var(targets));assert abs(floor-.0625)<1e-15
cases=pd.read_csv(scores/'case_scores.csv');assert np.allclose((cases.prediction-cases.truth)**2,cases.rectangle_SE,rtol=1e-12,atol=1e-14)
seed=cases.groupby(['seed','method'],as_index=False).rectangle_SE.mean();saved=pd.read_csv(scores/'seed_scores.csv');merged=seed.merge(saved,on=['seed','method']);np.testing.assert_allclose(merged.rectangle_SE,merged.four_world_rectangle_MSE,rtol=1e-12,atol=1e-14)
summary=[]
for method,ss in seed.groupby('method'):
 v=ss.rectangle_SE.to_numpy();se=v.std(ddof=1)/np.sqrt(len(v));d=t.ppf(.975,len(v)-1)*se;summary.append(dict(method=method,seeds=len(v),four_world_rectangle_MSE=v.mean(),MCSE=se,MC95_low=v.mean()-d,MC95_high=v.mean()+d))
ss=pd.DataFrame(summary);expected=pd.read_csv(scores/'summary.csv').set_index('method')
for x in summary:
 for key in ['four_world_rectangle_MSE','MCSE','MC95_low','MC95_high']:assert np.isclose(x[key],expected.loc[x['method'],key],rtol=1e-11,atol=1e-13),(x,key)
ss.to_csv(out/'THEORY_SUMMARY.csv',index=False)
# Recompute all repeat and pooled growth errors after endpoint projection.
df=pd.read_csv(R/'paper/evidence/source_data/publication/Horlbeck_predictions.csv');ref=pd.read_csv(R/'paper/evidence/source_data/publication/Horlbeck_figure_metrics.csv');rows=[]
models=['M0','M0_repeat_mean','rank1_selected','spectral_selected']
for (stage,budget,fold),group in df.groupby(['stage','budget','fold']):
 nodes=sorted(set(group.endpoint_a)|set(group.endpoint_b));ix={x:i for i,x in enumerate(nodes)}
 for rep in ['1','2','pooled']:
  cur=group if rep=='pooled' else group[group['repeat']==int(rep)];full={};projected={}
  for model in models:
   fs=ps=0.
   for _,gg in cur.groupby('repeat'):
    X=np.zeros((len(gg),len(nodes)));X[np.arange(len(gg)),[ix[x] for x in gg.endpoint_a]]=1;X[np.arange(len(gg)),[ix[x] for x in gg.endpoint_b]]=1;e=gg.observed.to_numpy()-gg[model].to_numpy();pe=e-X@np.linalg.lstsq(X,e,rcond=None)[0];fs+=float(e@e);ps+=float(pe@pe)
   full[model]=fs;projected[model]=ps
  for model in models:
   r=ref[(ref.stage==stage)&(ref.budget==budget)&(ref.fold.astype(str)==str(fold))&(ref['repeat'].astype(str)==rep)&(ref.model==model)]
   if r.empty:continue
   assert len(r)==1;rr=r.iloc[0];assert np.isclose(full[model],rr.sse,rtol=1e-11,atol=1e-12);assert np.isclose(projected[model],rr.projected_sse,rtol=1e-11,atol=1e-12)
   rows.append(dict(stage=stage,budget=budget,fold=int(fold),repeat=rep,model=model,sse=full[model],projected_sse=projected[model],gain=1-full[model]/full['M0'],projected_gain=1-projected[model]/projected['M0']))
assert len(rows)==84,('Missing functional score comparisons',len(rows))
for (stage,budget,rep,model),gg in pd.DataFrame(rows).groupby(['stage','budget','repeat','model']):
 rr=ref[(ref.stage==stage)&(ref.budget==budget)&(ref.fold=='pooled')&(ref['repeat'].astype(str)==rep)&(ref.model==model)]
 if rr.empty:continue
 assert len(rr)==1
 for field in ['sse','projected_sse']:assert np.isclose(gg[field].sum(),rr.iloc[0][field],rtol=1e-11,atol=1e-12)
 rows.append(dict(stage=stage,budget=budget,fold='pooled',repeat=rep,model=model,sse=gg.sse.sum(),projected_sse=gg.projected_sse.sum()))
assert len(rows)==len(ref),('Functional score coverage',len(rows),len(ref))
pd.DataFrame(rows).to_csv(out/'GROWTH_SCORES.csv',index=False)
report=dict(status='PASS',theory_targets=targets,theory_lower_bound=floor,theory_case_rows=len(cases),growth_prediction_rows=len(df),growth_score_rows=len(rows),base_model_fits=0,scope='Exact family identities and scoring saved synthetic predictions; independent least-squares endpoint projection of saved functional predictions. No new model fitting.')
(out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
