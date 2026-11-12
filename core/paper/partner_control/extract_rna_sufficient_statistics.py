#!/usr/bin/env python3
"""One-time extraction from existing RNA histories; no fitting or selection.

Numerical replay uses RNA_SUFFICIENT_STATISTICS.csv and does not require these
historical source directories. Primary eleven-rectangle RNA results are untouched.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent;W=R.parents[1]
sys.path.insert(0,str(R/'software'))
from witnesscell_partner_control import project

def graph(names):
    pairs=[str(s).split('+') for s in names];nodes=sorted({s for p in pairs for s in p});ix={v:i for i,v in enumerate(nodes)}
    return np.array([[ix[a],ix[b]] for a,b in pairs]),len(nodes)
def stats(y,b,p,e,n):
    he,pe,dim=project(e,n,y-b);hd,pd,dim=project(e,n,p-b)
    return dict(endpoint_E=float(np.sum(he*he)),endpoint_A=float(np.sum(he*hd)),endpoint_B=float(np.sum(hd*hd)),partner_E=float(np.sum(pe*pe)),partner_A=float(np.sum(pe*pd)),partner_B=float(np.sum(pd*pd)),contrast_dim=dim)

def main():
    d=pd.read_csv(R/'inputs/ARCHIVED_ENRICHED_INPUT.csv')
    d=d[d.panel.eq('CaRPool') & d.source.str.endswith('/predictions.npz') & d.inner_origin.str.endswith('/selection.json')]
    U=W/'gain_regression_history/WitnessCell_R_Upstream_Selection_20261101/WitnessCell_R_Upstream_Selection_20261101/engine/source/adapters'
    a=np.load(U/'original_moments.npz');mask=a['score_gene_mask'].astype(bool);genes=a['genes'][mask].tolist();gi={g:i for i,g in enumerate(genes)};ci={s:i for i,s in enumerate(a['conditions'])};means=a['means'][:,mask]
    sources={str(U/'original_moments.npz'):hashlib.sha256((U/'original_moments.npz').read_bytes()).hexdigest()};rows=[]
    for src,group in d.groupby('source'):
        path=Path(src);outer=np.load(path);selection=json.loads((path.parent/'selection.json').read_text());ins=sorted(path.parent.glob('inner_selected_*.npz'))
        sources[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        by={}
        for ip in ins:
            t=np.load(ip);e,n=graph(t['conditions']);b=t['base'];y=t['y']
            if 'model_names' in t.files:
                pred={str(k):t['selected_prediction_values'][j] for k,j in zip(t['model_names'],t['model_prediction_index'])}
            else:pred={k[5:]:t[k] for k in t.files if k.startswith('pred_')}
            pred.update(R_original=t['R'],R=t['R'])
            for name,v in pred.items():by.setdefault(name,[]).append(stats(y,b,v,e,n))
            sources[str(ip)]=hashlib.sha256(ip.read_bytes()).hexdigest()
        names=outer['conditions'].tolist();cols=[gi[g] for g in outer['genes'].tolist()];y=means[[ci[s] for s in names]][:,cols]-outer['control'];b=outer['pred_effect_M0'];e,n=graph(names)
        for row in group.itertuples():
            if row.model not in by or len(by[row.model])!=3:raise ValueError('inner record unavailable: '+row.model)
            p=outer['pred_effect_'+row.model];ist={k:sum(s[k] for s in by[row.model]) for k in by[row.model][0]};ost=stats(y,b,p,e,n)
            aa=ist['partner_A'];bb=ist['partner_B'];alpha=float(np.clip(aa/bb,0,1)) if bb>1e-14 else 0.
            ig=100*(2*aa-bb)/ist['partner_E'];og=100*(2*ost['partner_A']-ost['partner_B'])/ost['partner_E']
            assert abs(ig-row.inner_partner_gain)<1e-7 and abs(og-row.gain_pp)<1e-7
            choice=selection.get('selected',{}).get(row.model,{})
            rows.append(dict(record_id=row.record_id,feature_id=row.feature_id,model=row.model,pool_id=row.pool_id,source=src,
                actual_family=choice.get('family',''),chosen_eta=choice.get('eta',np.nan),alpha=alpha,
                **{'inner_'+k:v for k,v in ist.items()},**{'outer_'+k:v for k,v in ost.items()},
                original_partner_gain_pp=og,inner_replay_error=abs(ig-row.inner_partner_gain),outer_replay_error=abs(og-row.gain_pp)))
    pd.DataFrame(rows).to_csv(R/'inputs/RNA_SUFFICIENT_STATISTICS.csv',index=False)
    (R/'inputs/RNA_SOURCE_HASHES.json').write_text(json.dumps(sources,indent=2)+'\n')
    print('RNA matched archive records',len(rows),'across source files',len(d.source.unique()))

if __name__=='__main__':main()
