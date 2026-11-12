"""Independent direct-design checks for saturation baselines and gain accounting."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,required=True);ap.add_argument('--expected',type=Path);args=ap.parse_args()
checks=[]
for label in json.loads((R/'ANALYSIS_PLAN.json').read_text())['configurations']:
    tr=np.load(R/'inputs'/label/'training.npz');pred=np.load(R/'fitted'/label/'PREDICTIONS.npz')
    singles=tr['singles'];n=len(singles);scale=np.sqrt(np.mean(singles**2,axis=0));scale[scale==0]=1
    for name,mask,qe,expect in [('outer',np.ones(len(tr['edges']),dtype=bool),tr['query_edges'],pred['B'])]+[(f'inner{k}',tr['folds']!=k,tr['edges'][tr['folds']==k],pred['B_oof'][tr['folds']==k]) for k in range(3)]:
        state=json.loads((R/'fitted'/label/('baseline.json' if name=='outer' else name+'_baseline.json')).read_text());ed=tr['edges'][mask];y=tr['y'][mask];lam=state['lambda_value'];outputs=[]
        for j in range(y.shape[1]):
            z=singles[ed[:,0],j]+singles[ed[:,1],j];f=z/(1+lam*abs(z)/scale[j])-z
            x=np.zeros((len(ed),n));x[np.arange(len(ed))[:,None],ed]=1
            # A single augmented least-squares solve, independent of residualization.
            coef=np.linalg.lstsq(np.column_stack([x,f]),y[:,j]-z,rcond=None)[0]
            zq=singles[qe[:,0],j]+singles[qe[:,1],j];fq=zq/(1+lam*abs(zq)/scale[j])-zq
            outputs.append(zq+coef[qe[:,0]]+coef[qe[:,1]]+coef[-1]*fq)
        gap=float(np.max(abs(np.column_stack(outputs)-expect)));assert gap<1e-9,(label,name,gap)
        checks.append(dict(configuration=label,fit=name,direct_baseline_max_abs=gap))
ha=pd.read_csv(R/'results/HAN_ATTRIBUTION_POOLED.csv');old=pd.read_csv(args.base/'experiments/han_native/results/POOLED.csv');gaps=[]
for row in ha.itertuples():
    x=old[(old.panel==row.panel)&(old.baseline==row.baseline)&(old.method==row.family+'_original')].iloc[0]
    gaps.append(abs(row.total_gain_pct-x.partner_gain_vs_own_baseline_pct))
assert max(gaps)<1e-9
comparison=[]
if args.expected:
    for name in ['COMMON_REFERENCE.csv','SATURATION_SCORES.csv','SATURATION_POOLED.csv','SATURATION_BY_CONFIGURATION.csv','SATURATION_BY_REPEAT.csv','HAN_DIRECTION_ATTRIBUTION.csv','HAN_ATTRIBUTION_POOLED.csv','HAN_ATTRIBUTION_BY_PARTITION.csv','HAN_ATTRIBUTION_BY_REPEAT.csv']:
        a=pd.read_csv(R/'results'/name);b=pd.read_csv(args.expected/name);assert a.shape==b.shape and list(a.columns)==list(b.columns)
        mx=0.
        for c in a:
            if pd.api.types.is_numeric_dtype(a[c]):
                assert np.allclose(a[c],b[c],rtol=1e-10,atol=1e-10,equal_nan=True),(name,c)
                mx=max(mx,float(abs(a[c]-b[c]).max()))
            else:assert a[c].fillna('').equals(b[c].fillna('')),(name,c)
        comparison.append(dict(table=name,max_abs_difference=mx))
report=dict(status='PASS',independent_baselines=checks,han_gain_max_abs_difference=max(gaps),replay_tables=comparison)
(R/'results/INDEPENDENT_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
