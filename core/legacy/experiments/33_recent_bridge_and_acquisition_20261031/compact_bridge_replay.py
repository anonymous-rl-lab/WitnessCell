"""Independent compact replay of the first bridge and its retrospective.

The frozen final predictions and scores are reproduced. Optional inner-CV G
matrix caches may be omitted, with their original SHA256 explicitly reported.
This does not claim that those caches have been refitted or regenerated.
No original scientific source or historical lock is modified.

python compact_bridge_replay.py --bridge /restored/WitnessCell_NM_NCS_Bridge_20261031 \
    --retrospective /restored/bridge_retrospective_20261031 --out /new/work
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse,csv,hashlib,json,shutil,subprocess,sys
from pathlib import Path
import numpy as np


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())


def verify_retained(bridge):
    record=read(bridge/'RELEASE_MANIFEST.json');omitted=[];checked=0
    for name,expected in record['files'].items():
        path=bridge/name
        if path.exists():
            assert sha(path)==expected,('changed retained file',name)
            checked+=1
        else:
            assert '/inner_models/' in name and name.endswith('.npz'),('unexpected missing file',name)
            omitted.append({'path':name,'original_sha256':expected,
                            'role':'reconstructible inner-CV matrix cache; not required for frozen prediction replay'})
    return {'retained_original_files_checked':checked,'omitted_inner_cache_files':len(omitted),
            'omissions':omitted,'all_retained_hashes_match':True,
            'all_original_files_present':not omitted}


def numeric_equal(a,b,label):
    if isinstance(a,(int,float)) and not isinstance(a,bool) and isinstance(b,(int,float)):
        assert np.isclose(a,b,rtol=1e-12,atol=1e-12,equal_nan=True),(label,a,b)
    elif isinstance(a,dict):
        assert set(a)==set(b),(label,'keys')
        for k in a:
            if k in ('source','opened_npz','npz_opened','elapsed_seconds'):continue
            numeric_equal(a[k],b[k],label+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(label,'length')
        for i,(x,y) in enumerate(zip(a,b)):numeric_equal(x,y,label+'/'+str(i))
    else: assert a==b,(label,a,b)


def compare_csv(a,b):
    with a.open(newline='') as f:x=list(csv.DictReader(f))
    with b.open(newline='') as f:y=list(csv.DictReader(f))
    assert len(x)==len(y),(str(a),'rows')
    for i,(u,v) in enumerate(zip(x,y)):
        assert set(u)==set(v)
        for k in u:
            if u[k]==v[k]:continue
            try: aa,bb=float(u[k]),float(v[k])
            except ValueError:raise AssertionError((str(a),i,k,u[k],v[k]))
            assert np.isclose(aa,bb,rtol=1e-12,atol=1e-12,equal_nan=True),(str(a),i,k,aa,bb)
    return len(x)


def growth_replay(bridge):
    run=bridge/'growth_run'
    assets=bridge/'growth_assets/extracted/WitnessCell_Horlbeck_TwoStage_Repro_20261018/stages/K562'
    sys.path.insert(0,str(assets/'core'))
    from frozen_models import FittedModel
    from run_smoke import corr_prediction
    locked=read(run/'PREDICTIONS_LOCKED_v2.json');assert locked['completed_pipelines']==45
    missing=[]
    for name,expected in locked['files'].items():
        p=run/name
        if p.exists():assert sha(p)==expected,name
        else:
            assert '/inner_models/' in name and name.endswith('.npz'),name
            missing.append(name)
    methods=('M0','M0_repeat_mean','rank1_selected','spectral_selected','rank1_regularized','spectral_regularized')
    oldrows=read(run/'scoring_v2/TASK_METRICS.json')
    maxdiff=0.;pipelines=0;score_rows=0
    for fold in range(3):
        with np.load(assets/f'cases/outer{fold}/inputs/training.npz',allow_pickle=False) as z:singles=z['singles'].copy()
        scales=np.std(singles,axis=0);scales[scales==0]=1
        with np.load(assets/f'cases/outer{fold}/inputs/evaluation.npz',allow_pickle=False) as z:
            edges=z['edges'].copy();ids=z['global_ids'].copy();truth=z['y'].copy()
        X=np.zeros((len(edges),len(singles)));X[np.arange(len(edges)),edges[:,0]]=1;X[np.arange(len(edges)),edges[:,1]]=1
        U,sv,_=np.linalg.svd(X,full_matrices=False)
        rank=int(np.sum(sv>max(X.shape)*np.finfo(float).eps*sv[0]));U=U[:,:rank]
        def project(x):return x-U@(U.T@x)
        for seed in range(5):
            task=f'fold_{fold}_seed_{seed}'
            for arm in ('high_local_support','low_local_support','random_degree_switch_chain'):
                p=run/f'training_v2/{task}/results/{arm}'
                models=[FittedModel.from_dict(x) for x in read(p/'M0_models.json')]
                selection=read(p/'selection.json')
                with np.load(run/f'layouts_v2/{task}/{arm}.npz',allow_pickle=False) as z:train=z['edges'].copy();query=z['query_edges'].copy()
                with np.load(p/'predictions.npz',allow_pickle=False) as z:saved={k:z[k].copy() for k in z.files}
                assert np.array_equal(edges,query) and np.array_equal(ids,saved['global_ids'])
                base=np.column_stack([models[t].predict(query)*scales[t] for t in range(2)])
                replay={'M0':base,'M0_repeat_mean':np.repeat(base.mean(axis=1)[:,None],2,axis=1)}
                for family in ('rank1','spectral'):
                    reg=base+corr_prediction(saved[family+'_G'],train,query,len(singles))[:,None]
                    replay[family+'_regularized']=reg
                    replay[family+'_selected']=base.copy() if selection[family]['selected']=='G0' else reg
                baseline_error=truth-base
                bfull=np.sum(baseline_error**2,axis=0);bpartner=np.sum(project(baseline_error)**2,axis=0)
                for method in methods:
                    value=replay[method];diff=float(np.max(np.abs(value-saved[method])));maxdiff=max(maxdiff,diff)
                    np.testing.assert_array_equal(value,saved[method])
                    error=truth-value;fs=np.sum(error**2,axis=0);ps=np.sum(project(error)**2,axis=0)
                    for rep in (0,1,'pooled'):
                        old=next(x for x in oldrows if x['fold']==fold and x['seed']==seed and x['arm']==arm and x['model']==method and x['repeat']==rep)
                        take=lambda v:float(v.sum() if rep=='pooled' else v[rep])
                        actual={'full_sse':take(fs),'partner_sse':take(ps),'M0_full_sse':take(bfull),'M0_partner_sse':take(bpartner),'net_partner_gain':take(bpartner)-take(ps)}
                        for k,v in actual.items():numeric_equal(v,old[k],f'{task}/{arm}/{method}/{rep}/{k}')
                        score_rows+=1
                pipelines+=1
    assert pipelines==45 and score_rows==len(oldrows)
    return {'final_pipelines':pipelines,'prediction_outputs':pipelines*len(methods),'score_rows':score_rows,
            'max_prediction_abs_difference':maxdiff,'inner_caches_not_replayed':len(missing),'no_fitting':True}


def main():
    p=argparse.ArgumentParser();p.add_argument('--bridge',type=Path,required=True)
    p.add_argument('--retrospective',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();src=args.bridge.resolve();retro=args.retrospective.resolve();out=args.out.resolve()
    assert not out.exists(),'Use a new output directory'
    initial=verify_retained(src)
    out.mkdir(parents=True)
    bridge=out/'bridge_experiment_20261031';ret=out/'bridge_retrospective_20261031'
    shutil.copytree(src,bridge);shutil.copytree(retro,ret)
    commands=[]
    def execute(root,script,*arguments):
        log=out/f'command_{len(commands):02d}.log'
        with log.open('w') as f:
            proc=subprocess.run([sys.executable,str(root/script),*arguments],cwd=root,stdout=f,stderr=subprocess.STDOUT)
        commands.append({'script':str((root/script).relative_to(out)),'returncode':proc.returncode,'log':log.name})
        assert proc.returncode==0,log.read_text()[-8000:]
    growth=growth_replay(bridge)
    causal=bridge/'causal_run';expected=causal/'scores_expected'
    (causal/'scores').rename(expected)
    execute(causal,'score_with_amendment.py')
    score_files=[];npz_arrays=0
    for old in sorted(expected.iterdir()):
        new=causal/'scores'/old.name
        if not new.exists():continue  # Later report-only additions are not scorer outputs.
        if old.suffix=='.csv':compare_csv(old,new);score_files.append(old.name)
        elif old.suffix=='.json':numeric_equal(read(old),read(new),old.name);score_files.append(old.name)
        elif old.suffix=='.npz':
            with np.load(old,allow_pickle=False) as x,np.load(new,allow_pickle=False) as y:
                assert x.files==y.files
                for k in x.files:assert np.array_equal(x[k],y[k]),(old.name,k);npz_arrays+=1
    assert (causal/'scores/SCORE_MANIFEST.json').exists()
    execute(ret,'causal_audit/diagnose.py')
    execute(ret,'causal_audit/crn_structure_audit.py')
    execute(ret,'growth_audit/audit_frozen_design.py')
    execute(ret,'growth_audit/exact_invariant.py')
    compared=[]
    for folder in ('causal_audit','growth_audit'):
        for old in sorted((retro/folder).glob('*')):
            new=ret/folder/old.name
            if not new.exists():continue
            if old.suffix=='.csv':compare_csv(old,new);compared.append(str(old.relative_to(retro)))
            elif old.suffix=='.json' and old.name not in ('MANIFEST.json','AUDIT_PROVENANCE.json'):
                numeric_equal(read(old),read(new),str(old.relative_to(retro)));compared.append(str(old.relative_to(retro)))
    result={'status':'PASS','verification_scope':'retained original hashes + frozen final weights to predictions + scores; no inner-cache refit claim',
            'compression_provenance':initial,'growth':growth,'causal':{'replayed_models':144,'matching_score_files':score_files,'truth_arrays_bitwise_equal':npz_arrays},
            'retrospective_artifacts_compared':compared,'commands':commands}
    (out/'COMPACT_BRIDGE_REPLAY_RESULT.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ('compression_provenance','causal','commands')},indent=2))


if __name__=='__main__':main()
