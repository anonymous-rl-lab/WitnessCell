"""Recompute the full neural comparison from frozen predictions; NumPy only. Requires NumPy only.

Initial export: python score.py --source PATH_TO_GSE220974_RUN
Portable replay: python score.py
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
def load(p):
    with np.load(p, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}
def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def export(source):
    m = load(source / 'data/moments.npz')
    conditions = np.char.replace(m['conditions'], '|', '+')
    ids = np.flatnonzero(np.char.find(conditions, '+') >= 0)
    pairs = conditions[ids]; mask = m['score_gene_mask']
    n, t, g = len(ids), len(m['target_names']), int(mask.sum())
    data = dict(conditions=pairs, targets=m['target_names'], genes=m['genes'][mask],
        objectives=m['objectives'][:, mask], truth=m['means'][ids][:, mask],
        half_U=np.array([m[k][ids] @ m['objectives'].T for k in ('first_means','second_means')]),
        fold=np.full(n,-1,int))
    records=[]
    for fold in range(5):
        for kind in ('classical', 'presage'):
            folder = source / f'results/detected49/{kind}/fold_{fold}'
            rec = json.loads((folder/'FREEZE.json').read_text())
            expected = rec.get('output_sha256', rec.get('predictions_sha256'))
            assert digest(folder/'predictions.npz') == expected
            z = load(folder/'predictions.npz'); q = z['query']
            assert np.array_equal(z['conditions'], pairs[q])
            assert np.array_equal(z['target_names'], data['targets'])
            if kind == 'classical':
                assert np.all(data['fold'][q] == -1); data['fold'][q] = fold
                assert not set(z['query']) & set(z['train'])
                if 'data/moments.npz' in rec['source_sha256']:
                    assert digest(source/'data/moments.npz') == rec['source_sha256']['data/moments.npz']
                us = {k[8:]:v for k,v in z.items() if k.startswith('utility_')}
                rs = {k[12:]:v for k,v in z.items() if k.startswith('pred_effect_')}
            else:
                assert np.all(data['fold'][q] == fold)
                assert not set(rec['query_ids']) & set(rec['training_ids'])
                us = {'PRESAGE_Combo_seed24':z['utility']}
                rs = {'PRESAGE_Combo_seed24':z['pred_effect'][:,mask]}
            for prefix, items, width in [('U__',us,t),('RNA__',rs,g)]:
                for name, values in items.items():
                    key = prefix+name
                    if key not in data: data[key] = np.full((n,width),np.nan)
                    assert np.isnan(data[key][q]).all(); data[key][q] = values
            records.append(dict(kind=kind, fold=fold, prediction_sha256=expected))
    for kind in ('gears','scouter'):
        prefix='GEARS' if kind=='gears' else 'Scouter'
        for seed in range(3):
            name=prefix+'_seed'+str(seed)
            data['U__'+name]=np.full((n,t),np.nan)
            data['RNA__'+name]=np.full((n,g),np.nan)
            for fold in range(5):
                folder=source/f'results/detected49/{kind}_seed{seed}/fold_{fold}'
                rec=json.loads((folder/'FREEZE.json').read_text());expected=rec['predictions_sha256']
                assert digest(folder/'predictions.npz')==expected
                z=load(folder/'predictions.npz');q=z['query']
                assert np.array_equal(z['conditions'],pairs[q])
                assert np.array_equal(z['target_names'],data['targets'])
                assert np.array_equal(z['genes'],m['genes'])
                assert np.all(data['fold'][q]==fold)
                assert not set(rec['query']) & set(rec['train'])
                assert np.isnan(data['U__'+name][q]).all()
                data['U__'+name][q]=z['utility'];data['RNA__'+name][q]=z['pred_effect'][:,mask]
                records.append(dict(kind=kind,seed=seed,fold=fold,prediction_sha256=expected))
        for key in ('U__','RNA__'):
            data[key+prefix+'_Ensemble']=np.mean([data[key+prefix+'_seed'+str(seed)] for seed in range(3)],axis=0)
    assert np.all(data['fold'] >= 0)
    for k,v in data.items():
        if k.startswith(('U__','RNA__')): assert np.isfinite(v).all()
    np.savez_compressed(ROOT/'scoring_inputs.npz', **data)
    (ROOT/'SOURCE_CHECKS.json').write_text(json.dumps(dict(status='PASS',
        prediction_files=records, source_package='WitnessCell_NextPro2_Continuation_20261102.zip',
        scope='All 30 neural jobs complete; 15 result rows including historical and individual-seed comparators',
        evaluation='Previously opened outcomes, recomputed 2026-11-05; no post-score tuning'),indent=2)+'\n')

def score():
    d=load(ROOT/'scoring_inputs.npz'); names=d['targets'].astype(str); pairs=d['conditions'].astype(str)
    x=np.array([[name in pair.split('+') for name in names] for pair in pairs],float)
    truth=d['truth']@d['objectives'].T
    geometry=[]; groups=[]
    for f in range(5):
        ids=np.flatnonzero(d['fold']==f); rank=np.linalg.matrix_rank(x[ids])
        left,_,_=np.linalg.svd(x[ids],full_matrices=True); basis=left[:,rank:]
        proj=np.eye(len(ids))-x[ids]@np.linalg.pinv(x[ids])
        np.testing.assert_allclose(basis@basis.T,proj,atol=1e-12)
        geometry.append((ids,basis,proj))
        for a,anchor in enumerate(names):
            pool=ids[x[ids,a]>0]
            if len(pool)>=2:
                pool=np.array(sorted(pool,key=lambda i:pairs[i]))
                groups.append((f,a,anchor,pool))
    models=sorted(k[3:] for k in d if k.startswith('U__')); assert len(models)==15
    records=[];summary=[];chosen_by_model={};checks=0
    degrees=sum(b.shape[1] for _,b,_ in geometry)
    for model in models:
        pred=d['U__'+model]; rows=[];chosen_by_model[model]=[]
        for f,a,anchor,pool in groups:
            utility=pred[pool,a]; chosen=pool[np.flatnonzero(utility>=utility.max()-1e-12)[0]]
            actual=truth[chosen,a]; best=truth[pool,a].max()
            row=dict(model=model,group=f'f{f}:{anchor}',fold=f,anchor=anchor,
                physical_anchor=anchor.rsplit('-',1)[0],candidate_count=len(pool),pair=pairs[chosen],
                actual_U=float(actual),best_U=float(best),regret=float(best-actual),
                hit=int(actual>=best-1e-12),half0_U=float(d['half_U'][0,chosen,a]),
                half1_U=float(d['half_U'][1,chosen,a]))
            rows.append(row);records.append(row);chosen_by_model[model].append(int(chosen))
        large=[r for r in rows if r['candidate_count']>=4]
        rec=dict(model=model,tasks=len(rows),mean_U=float(np.mean([r['actual_U'] for r in rows])),
            mean_regret=float(np.mean([r['regret'] for r in rows])),hits=sum(r['hit'] for r in rows),
            RNA_MSE='',partner_MSE='',partner_degrees=degrees,large_tasks=len(large),
            large_regret=float(np.mean([r['regret'] for r in large])),
            half0_U=float(np.mean([r['half0_U'] for r in rows])),half1_U=float(np.mean([r['half1_U'] for r in rows])))
        if 'RNA__'+model in d:
            response=d['RNA__'+model]
            np.testing.assert_allclose(response@d['objectives'].T,pred,rtol=1e-7,atol=1e-10)
            err=response-d['truth'];rec['RNA_MSE']=float(np.mean(err**2))
            svd_score=sum(np.sum((b.T@err[ids])**2) for ids,b,_ in geometry)/(degrees*err.shape[1])
            proj_score=sum(np.sum((p@err[ids])**2) for ids,_,p in geometry)/(degrees*err.shape[1])
            np.testing.assert_allclose(svd_score,proj_score,rtol=1e-12,atol=1e-14);checks+=1
            rec['partner_MSE']=float(svd_score)
        summary.append(rec)
    summary.sort(key=lambda r:r['mean_regret'])
    write_csv(ROOT/'SUMMARY.csv',summary);write_csv(ROOT/'CHOICES.csv',records)
    lookup={r['model']:r for r in summary};rg=lookup['R_G'];pr=lookup['PRESAGE_Combo_seed24'];tr=lookup['PairPick_TargetRank'];ad=lookup['SingleB_additive']
    identical=sum(a==b for a,b in zip(chosen_by_model['R_G'],chosen_by_model['PRESAGE_Combo_seed24']))
    facts=dict(models=len(models),tasks=len(groups),physical_anchors=len(set(a.rsplit('-',1)[0] for _,_,a,_ in groups)),
        RG_PRESAGE_identical_choices=identical,
        RG_regret_reduction_vs_additive_percent=100*(1-rg['mean_regret']/ad['mean_regret']),
        RG_RNA_reduction_vs_PRESAGE_percent=100*(1-rg['RNA_MSE']/pr['RNA_MSE']),
        TargetRank_regret_increase_vs_RG_percent=100*(tr['mean_regret']/rg['mean_regret']-1),
        RG_partner_reduction_vs_additive_percent=100*(1-rg['partner_MSE']/ad['partner_MSE']),
        projection_independent_checks=checks,partner_degrees=degrees)
    main_names=['R_G','GEARS_Ensemble','Scouter_Ensemble','StandardRidge_RNA','SingleB_additive']
    write_csv(ROOT/'MAIN_FIVE_MODELS.csv',[lookup[n] for n in main_names])
    for reference in ['GEARS_Ensemble','Scouter_Ensemble','StandardRidge_RNA','SingleB_additive']:
        rr=lookup[reference]
        facts['RG_vs_'+reference]={metric:100*(1-rg[metric]/rr[metric]) for metric in ['mean_regret','RNA_MSE','partner_MSE']}
    (ROOT/'FACTS.json').write_text(json.dumps(facts,indent=2)+'\n')
    print((ROOT/'SUMMARY.csv').read_text());print(json.dumps(facts,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path);a=ap.parse_args()
    if a.source:export(a.source)
    score()
