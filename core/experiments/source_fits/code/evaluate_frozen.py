"""Independent scoring after both finite candidate batches were frozen."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import json
from pathlib import Path
import numpy as np
import pandas as pd
from fit_models import ROOT,load_input,dump,sha

cells=['HeLa','A549','293T','Jurkat'];out=ROOT/'evaluation';out.mkdir(exist_ok=True)
choices={}
for cell in cells:
    opts=[]
    for folder in ['results','results_phase2']:
        p=ROOT/folder/cell;fr=json.loads((p/'FREEZE.json').read_text())
        for name,h in fr['files'].items():assert sha(p/name)==h,(cell,name)
        s=json.loads((p/'SELECTION_BEFORE_OUTER.json').read_text())
        opts.append(dict(phase=folder,**s['selected']['selected']))
    choices[cell]=min(opts,key=lambda x:(x['partner_sse'],x['rna_sse'],x['key']))
dump(out/'COMBINED_SELECTION_BEFORE_OUTER.json',{'rule':'lowest pooled inner partner SSE among the two already frozen candidate batches','selected':choices,'new_outer_predictions_scored':False})
rows=[];checks=[];details=[]
formal=ROOT.parent/'shen_takeover/Shen2017_Formal_20261108'
hor=ROOT/'sources/horlbeck/WitnessCell_Horlbeck_TwoStage_Repro_20261018/stages/Jurkat'
for cell in cells:
    d=load_input(ROOT/'inputs'/f'{cell}.json');q=d['query'];n=len(d['genes']);m=len(q)
    if (ROOT/'references'/f'{cell}_original.npz').exists():
        original=np.load(ROOT/'references'/f'{cell}_original.npz');assert np.array_equal(original['genes'],d['genes'])
        assert np.array_equal(original['edges'],q);y=original['y']
    elif cell=='Jurkat':
        truth=np.load(hor/'sealed/evaluation_truth.npz');assert np.array_equal(truth['edges'],q);y=truth['y']
        original=np.load(hor/'results/high/predictions.npz')
    else:
        truth=np.load(formal/'prepared/sealed_evaluation'/f'{cell}.npz');assert np.array_equal(truth['genes'],d['genes'])
        assert np.array_equal(truth['pairs'],d['genes'][q]);y=truth['y']
        original=np.load(formal/'runs/formal'/cell/'results/high/predictions.npz')
    assert np.array_equal(original['edges'],q)
    z1=np.load(ROOT/'results'/cell/'PREDICTIONS.npz');z2=np.load(ROOT/'results_phase2'/cell/'PREDICTIONS.npz')
    base=z1['M0'];assert np.max(abs(base-z2['M0']))<1e-12
    assert np.max(abs(base-original['M0']))<1e-12
    x=np.zeros((m,n));x[np.arange(m)[:,None],q]=1
    def project(v):return v-x@np.linalg.lstsq(x,v,rcond=None)[0]
    e=y-base;pe=project(e);e2=np.sum(e*e);p2=np.sum(pe*pe)
    preds={'original_rank1':original['rank1_selected'],'original_spectral':original['spectral_selected'],
        'phase1_selected':z1['selected'],'next2_features':z1['next2_selected'],'projected_rank1':z1['projected_selected'],
        'native_decomposed':z1['native_decomp_selected'],'phase2_selected':z2['selected'],'sf12':z2['sf12_selected'],
        'network_v4':z2['network_selected'],'combined_selected':(z1 if choices[cell]['phase']=='results' else z2)['selected']}
    s2=json.loads((ROOT/'results_phase2'/cell/'SELECTION_BEFORE_OUTER.json').read_text())
    nr=s2['selected']['network_selected']
    if nr['key']!='M0':preds['network_without_calibration']=base+z2[nr['key']+'__uncalibrated_delta']
    for method,pred in preds.items():
        error=y-pred;partner=project(error);pdelta=project(pred-base)
        gain=100*(1-np.sum(partner**2)/p2)
        A=float(np.sum(pe*pdelta));B=float(np.sum(pdelta*pdelta));assert abs((p2-np.sum(partner**2))-(2*A-B))<1e-10
        row=dict(cell=cell,method=method,queries=m,full_gain_pct=float(100*(1-np.sum(error**2)/e2)),partner_gain_pct=float(gain),
            full_sse=float(np.sum(error**2)),partner_sse=float(np.sum(partner**2)),
            partner_gain_rep1=float(100*(1-np.sum(partner[:,0]**2)/np.sum(pe[:,0]**2))),
            partner_gain_rep2=float(100*(1-np.sum(partner[:,1]**2)/np.sum(pe[:,1]**2))),
            correction_cosine=A/np.sqrt(p2*B) if B>1e-24 else None,
            alignment_reward_2A=2*A,correction_cost_B=B,base_partner_sse=float(p2))
        rows.append(row)
    checks.append(dict(cell=cell,original_M0_max_abs=float(np.max(abs(base-original['M0']))),query_identity_equal=True,
        phase1_freeze_verified=True,phase2_freeze_verified=True,partner_projection_rank=int(np.linalg.matrix_rank(x))))
    np.savez_compressed(out/f'{cell}_scored_predictions.npz',genes=d['genes'],edges=q,y=y,M0=base,**preds)
pd.DataFrame(rows).to_csv(out/'SCORES.csv',index=False);dump(out/'CHECKS.json',checks)
print(pd.DataFrame(rows).query("method in ['original_rank1','next2_features','sf12','network_v4','combined_selected']")[['cell','method','full_gain_pct','partner_gain_pct','partner_gain_rep1','partner_gain_rep2']].to_string(index=False))
