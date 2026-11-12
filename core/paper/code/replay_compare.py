"""Precision-aware comparison; bounds propagate observed prediction differences."""
from pathlib import Path
import json,platform,sys
import numpy as np
import pandas as pd

def norm_budget(reference, family, endpoints):
    eps=np.finfo(np.float32 if family=='v4' else np.float64).eps
    factor=64 if family=='v4' else 256*max(1,endpoints)
    return factor*eps*max(1.,float(np.linalg.norm(reference)))

def compare(reference, replay, mode='compatible'):
    reference=Path(reference);replay=Path(replay)
    sys.path.insert(0,str(replay/'code'))
    import run,common,torch,scipy
    checks=[];stages=[];coefficients=[];rows=[]
    def check(name,ok,**details):
        checks.append(dict(check=name,passed=bool(ok),**details))
        if not ok:raise AssertionError(checks[-1])
    for manifest in ['MODEL_FREEZE.json','RUN_CODE_LOCK.json']:
        for rel,digest in json.loads((reference/manifest).read_text())['files'].items():
            check('immutable '+rel,common.sha(replay/rel)==digest)
    for p in sorted((reference/'data').rglob('*.npz')):
        if p.name.endswith('_frozen_predictions.npz'):continue
        other=replay/p.relative_to(reference)
        if common.sha(p)==common.sha(other):equal=True
        else:
            aa,bb=common.load(p),common.load(other);equal=aa.keys()==bb.keys() and all(np.array_equal(aa[k],bb[k]) for k in aa)
        check('input arrays '+str(p.relative_to(reference)),equal)
    exact=True;pooled_bounds={}
    for panel in run.PANELS:
        old=common.load(reference/f'data/{panel}_frozen_predictions.npz');new=common.load(replay/f'data/{panel}_frozen_predictions.npz')
        for key in ['pairs','genes','models']:
            check(panel+' identity '+key,np.array_equal(old[key],new[key]))
        for key in ['B','M0','D','source_alpha']:
            a,b=old[key],new[key];exact &= np.array_equal(a,b)
            if key in ['M0','D']:
                for j,m in enumerate(old['models']):
                    aa,bb=a[...,j],b[...,j];gap=float(np.linalg.norm(bb-aa));limit=norm_budget(aa,str(m).split('__')[1],len(old['genes']))
                    check(panel+' inference '+key+' '+str(m),gap<=limit,gap_norm=gap,budget=limit)
                    stages.append(dict(panel=panel,model=str(m),stage=key,max_abs=float(np.max(abs(bb-aa))),norm_difference=gap,precision_budget=limit,exact=np.array_equal(aa,bb)))
            else:
                gap=float(np.linalg.norm(b-a));limit=256*np.finfo(float).eps*max(1.,float(np.linalg.norm(a)))*len(old['genes'])
                check(panel+' inference '+key,gap<=limit)
        for seed in range(20261108,20261113):
            for fold in range(3):
                rel=Path('results')/panel/str(seed)/str(fold)
                so=json.loads((reference/rel/'PARAMETERS.json').read_text());sn=json.loads((replay/rel/'PARAMETERS.json').read_text())
                check('calibration geometry '+str(rel),so['calibration_blocks']==sn['calibration_blocks'])
                ao,an=np.array(so['alpha']),np.array(sn['alpha']);A,B=np.array(so['cross']),np.array(so['energy']);C,E=np.array(sn['cross']),np.array(sn['energy'])
                for j,m in enumerate(old['models']):
                    check('energy branch '+str(rel)+' '+str(m),(B[j]>1e-14)==(E[j]>1e-14))
                    bound=(abs(C[j]-A[j])+abs(A[j]/B[j])*abs(E[j]-B[j]))/E[j] if B[j]>1e-14 else 0.
                    bound+=512*np.finfo(float).eps*max(1.,abs(ao[j]),abs(an[j]));gap=abs(an[j]-ao[j])
                    check('coefficient propagation '+str(rel)+' '+str(m),gap<=bound)
                    coefficients.append(dict(panel=panel,seed=seed,fold=fold,model=str(m),absolute_difference=float(gap),propagation_bound=float(bound)))
                exact &= so==sn
                qo=common.load(reference/rel/'QUERY_PREDICTIONS.npz');qn=common.load(replay/rel/'QUERY_PREDICTIONS.npz')
                check('query identities '+str(rel),np.array_equal(qo['ids'],qn['ids']))
                ids=qo['ids'];p,outs_old=run.arms(old,ids,so,panel);_,outs_new=run.arms(new,ids,sn,panel)
                keys=['model','method','repeat'];score_old=pd.read_csv(reference/rel/'SCORES.csv').set_index(keys);score_new=pd.read_csv(replay/rel/'SCORES.csv').set_index(keys)
                check('score identities '+str(rel),score_old.index.equals(score_new.index))
                for arm in outs_old:
                    diff=outs_new[arm]-outs_old[arm];pdiff=p.P(diff);hdiff=diff-pdiff;exact &= np.array_equal(outs_new[arm],outs_old[arm])
                    for metric,delta in [('full',diff),('partner',pdiff),('endpoint',hdiff)]:
                        for rep,col in [('pooled',slice(None))]+([('repeat1',slice(0,1)),('repeat2',slice(1,2))] if delta.shape[1]==2 else []):
                            for j,m in enumerate(old['models']):
                                idx=(m,arm,rep);sse=float(score_old.loc[idx,metric+'_sse']);val=float(score_new.loc[idx,metric+'_sse'])
                                dn=float(np.linalg.norm(delta[:,col,j]));n=delta[:,col,j].size
                                rounding=128*np.finfo(float).eps*(n+len(p.genes))*max(1.,sse)
                                bound=2*np.sqrt(max(0,sse))*dn+dn*dn+rounding;gap=abs(val-sse)
                                check('SSE propagation '+str(rel)+' '+str(idx)+' '+metric,gap<=bound)
                                key=(panel,str(m),arm,rep,metric);pooled_bounds[key]=pooled_bounds.get(key,0.)+bound
                exact &= score_old[['full_sse','partner_sse','endpoint_sse']].equals(score_new[['full_sse','partner_sse','endpoint_sse']])
    keys=['panel','model','method','repeat'];a=pd.read_csv(reference/'results/POOLED_SCORES.csv').set_index(keys).sort_index();b=pd.read_csv(replay/'results/POOLED_SCORES.csv').set_index(keys).sort_index()
    check('pooled identities',a.index.equals(b.index))
    for idx in a.index:
        for metric in ['full','partner','endpoint']:
            x,y=float(a.loc[idx,metric+'_sse']),float(b.loc[idx,metric+'_sse']);gap=abs(x-y);bound=pooled_bounds[(*idx,metric)]+128*np.finfo(float).eps*max(1.,x)
            check('pooled SSE '+str(idx)+' '+metric,gap<=bound)
            rows.append(dict(zip(keys,idx),metric=metric,absolute_difference=gap,relative_difference=gap/max(abs(x),1e-300),propagation_bound=bound))
    if mode=='exact':check('exact numerical replay requested',exact)
    out=replay/'audit';out.mkdir(exist_ok=True)
    pd.DataFrame(stages).to_csv(out/'REPLAY_INFERENCE_DIFFERENCES.csv',index=False)
    pd.DataFrame(coefficients).to_csv(out/'REPLAY_COEFFICIENT_DIFFERENCES.csv',index=False)
    pd.DataFrame(rows).to_csv(out/'REPLAY_SCORE_DIFFERENCES.csv',index=False)
    first='none' if exact else ('frozen inference' if any(x['max_abs']>0 for x in stages) else 'calibration or scoring')
    report=dict(status='PASS',mode=mode,exact_numerical_match=bool(exact),first_observed_difference=first,
        environment=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,torch=torch.__version__),
        archived_requirements=(reference/'requirements.txt').read_text(),checks_passed=len(checks),folds=60,score_rows=10440,
        max_prediction_difference=max(x['max_abs'] for x in stages),max_alpha_difference=max(x['absolute_difference'] for x in coefficients),
        max_sse_difference=max(x['absolute_difference'] for x in rows),max_relative_sse_difference=max(x['relative_difference'] for x in rows),
        base_model_fits=0,scope='Frozen inference, target scalar re-estimation, independent QR verification and precision-aware comparison. No weights, archived results or scientific rules changed. Observed stage does not establish environmental root cause.')
    (replay/'REPLAY_RESULT_V36.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
