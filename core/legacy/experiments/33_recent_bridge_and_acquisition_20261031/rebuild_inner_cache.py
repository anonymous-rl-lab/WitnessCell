"""Rebuild omitted tuning and final caches using unchanged scientific models.

Rank-one matrices are exactly reconstructed from retained factor parameters.
Spectral matrices are refitted in their original warm-start order using retained
training inputs and baseline parameters. OOF/final predictions use the original
prediction function and saved selection. No evaluation labels are read.
Every NPZ must match the ORIGINAL prediction lock SHA256 before it is written.
Existing cache files are checked and are never replaced by different bytes.

Validated with Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and one BLAS thread.
Different numerical or ZIP-library environments may fail the exact hash check;
this is an explicit incompatibility, not permission to alter the original lock.

python rebuild_inner_cache.py --root /path/to/WitnessCell_NM_NCS_Bridge_20261031 --restore
python rebuild_inner_cache.py --root /path/to/bridge --out /new/cache/tree
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'
import argparse, hashlib, io, json, sys, time
from pathlib import Path
import numpy as np
import scipy

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def read(p):return json.loads(Path(p).read_text())

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',required=True,type=Path)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--restore',action='store_true')
    group.add_argument('--out',type=Path)
    parser.add_argument('--limit-arms',type=int)
    parser.add_argument('--predictions-only',action='store_true',
                        help='Rebuild OOF predictions from already restored G caches')
    parser.add_argument('--final-only',action='store_true',
                        help='Rebuild final predictions and final matrices only')
    parser.add_argument('--cache-root',type=Path,
                        help='Optional growth-run-shaped root containing restored inner G caches')
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    root=args.root.resolve();growth=root/'growth_run'
    assets=root/'growth_assets/extracted/WitnessCell_Horlbeck_TwoStage_Repro_20261018/stages/K562'
    plan=read(assets/'PLAN.json')
    for name,h in plan['core_sha256'].items():assert sha(assets/'core'/name)==h, name
    sys.path.insert(0,str(assets/'core'))
    from frozen_models import FittedModel
    from spectral_models import fit_spectral
    from run_smoke import corr_prediction
    lock=read(growth/'PREDICTIONS_LOCKED_v2.json')['files']
    expected={k:v for k,v in lock.items() if '/inner_models/' in k and k.endswith('.npz')}
    assert len(expected)==2160,len(expected)
    dest=growth if args.restore else args.out.resolve()
    if not args.restore:
        assert not dest.exists(),'Use a fresh destination'
        dest.mkdir(parents=True)
    begun=time.monotonic();results=[];count=0;original_array_comparisons=0;maxdiff=0.;pred_count=0;pred_arrays=0;final_count=0;final_arrays=0
    reports=sorted((growth/'training_v2').glob('fold_*_seed_*/results/*/inner_summary.json'))
    assert len(reports)==45,len(reports)
    if args.limit_arms:reports=reports[:args.limit_arms]
    for pp in reports:
        arm=pp.parent;task=arm.parent.parent.name;fold=int(task.split('_')[1])
        summary=read(pp);grid=summary['grid']
        with np.load(assets/f'cases/outer{fold}/inputs/training.npz',allow_pickle=False) as z:
            singles=z['singles']
        scales=np.std(singles,axis=0);scales[scales==0]=1.
        # This retained OOF archive stores exactly the authorized layout rows.
        data=load(arm/(grid[0]['key']+'_inner_predictions.npz'))
        edges,y,folds=data['edges'],data['y'],data['fold']
        n=len(singles);arm_started=time.monotonic();arm_count=0
        base_oof=np.zeros_like(y);predictions={row['key']:np.zeros_like(y) for row in grid}
        for f in ([] if args.final_only else range(3)):
            train=folds!=f;val=~train
            models=[FittedModel.from_dict(x) for x in read(arm/f'inner{f}_baseline.json')]
            base=np.column_stack([models[t].predict(edges[train])*scales[t] for t in range(2)])
            base_oof[val]=np.column_stack([models[t].predict(edges[val])*scales[t] for t in range(2)])
            residual=y[train]-base
            for nu in plan['nu']:
                warm=None
                for fraction in plan['tau_fractions']:
                    for family in ('rank1','spectral'):
                        key=f'{family}_t{fraction:g}_n{nu:g}'
                        row=next(x for x in grid if x['key']==key)
                        diag=next(x for x in row['diagnostics'] if x['inner_fold']==f)
                        rel=str((arm/'inner_models'/f'{key}_f{f}.npz').relative_to(growth))
                        if args.predictions_only:
                            cache=(args.cache_root or growth)/rel
                            assert sha(cache)==expected[rel],cache
                            G=load(cache)['G']
                        elif family=='rank1':G=diag['lambda']*np.outer(diag['u'],diag['u'])
                        else:
                            model=fit_spectral(n,edges[train],residual,diag['tau'],nu,
                                      max_iter=plan['optimizer_maxiter'],tol=plan['optimizer_tolerance'],
                                      deadline=None,initial=warm)
                            G=model['G'];warm=G
                            assert model['converged']==diag['converged'], (task,arm.name,key,f,'convergence')
                            assert model['iterations']==diag['iterations'], (task,arm.name,key,f,'iterations')
                        buff=io.BytesIO();np.savez_compressed(buff,G=G);raw=buff.getvalue()
                        actual=hashlib.sha256(raw).hexdigest()
                        original=growth/rel
                        if original.exists():
                            old=load(original)['G'];delta=float(np.max(np.abs(old-G)))
                            maxdiff=max(maxdiff,delta);original_array_comparisons+=1
                            assert np.array_equal(old,G),(rel,'array differs',delta)
                        assert actual==expected[rel], (rel,'NPZ SHA mismatch',actual,expected[rel])
                        target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True)
                        if target.exists():assert sha(target)==expected[rel],rel
                        else:target.write_bytes(raw)
                        correction=corr_prediction(G,edges[train],edges[val],n)
                        predictions[key][val]=base_oof[val]+correction[:,None]
                        count+=1;arm_count+=1
        for row in ([] if args.final_only else grid):
            key=row['key'];rel=str((arm/(key+'_inner_predictions.npz')).relative_to(growth))
            values={'edges':edges,'fold':folds,'y':y,'pred':predictions[key],'baseline':base_oof}
            buff=io.BytesIO();np.savez_compressed(buff,**values);raw=buff.getvalue()
            assert hashlib.sha256(raw).hexdigest()==lock[rel],(rel,'OOF NPZ SHA mismatch')
            original=growth/rel
            if original.exists():
                old=load(original)
                for k,a in values.items():
                    assert np.array_equal(old[k],a),(rel,k,'OOF array mismatch')
                    pred_arrays+=1
            target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True)
            if target.exists():assert sha(target)==lock[rel],rel
            else:target.write_bytes(raw)
            pred_count+=1
        layout=load(growth/'layouts_v2'/task/(arm.name+'.npz'))
        query_edges,query_ids=layout['query_edges'],layout['query_ids']
        models=[FittedModel.from_dict(x) for x in read(arm/'M0_models.json')]
        def baseline(ed):return np.column_stack([models[t].predict(ed)*scales[t] for t in range(2)])
        b=baseline(query_edges);residual=y-baseline(edges)
        selection=read(arm/'selection.json')
        payload={'edges':query_edges,'global_ids':query_ids,'M0':b,
                 'M0_repeat_mean':np.repeat(b.mean(axis=1)[:,None],2,axis=1)}
        for family in ('rank1','spectral'):
            diag=read(arm/(family+'_model.json'))
            if family=='rank1':G=diag['lambda']*np.outer(diag['u'],diag['u'])
            else:
                model=fit_spectral(n,edges,residual,diag['tau'],diag['nu'],
                           max_iter=plan['optimizer_maxiter'],tol=plan['optimizer_tolerance'],
                           deadline=None,initial=None)
                G=model['G']
                assert model['converged']==diag['converged'],(task,arm.name,'final convergence')
                assert model['iterations']==diag['iterations'],(task,arm.name,'final iterations')
            correction=corr_prediction(G,edges,query_edges,n)
            payload[family+'_regularized']=b+correction[:,None]
            payload[family+'_selected']=payload[family+'_regularized'] if selection[family]['selected']!='G0' else b.copy()
            payload[family+'_G']=G
        rel=str((arm/'predictions.npz').relative_to(growth))
        buff=io.BytesIO();np.savez_compressed(buff,**payload);raw=buff.getvalue()
        assert hashlib.sha256(raw).hexdigest()==lock[rel],(rel,'final NPZ SHA mismatch')
        original=growth/rel
        if original.exists():
            old=load(original)
            for k,a in payload.items():
                assert np.array_equal(old[k],a),(rel,k,'final array mismatch')
                final_arrays+=1
        target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():assert sha(target)==lock[rel],rel
        else:target.write_bytes(raw)
        final_count+=1
        result={'task':task,'arm':arm.name,'npz_sha256_exact':arm_count,
                'oof_npz_sha256_exact':0 if args.final_only else len(grid),
                'final_prediction_npz_sha256_exact':1,'seconds':time.monotonic()-arm_started}
        results.append(result);print(json.dumps(result),flush=True)
    report={'status':'PASS','mode':'exact factor reconstruction + original spectral warm-start refit',
       'root':str(root),'restored_cache_root':str(dest),'arms':len(results),
       'npz_files_sha256_exact':count,'original_array_comparisons':original_array_comparisons,
       'oof_npz_sha256_exact':pred_count,'oof_original_arrays_exact':pred_arrays,
       'final_prediction_npz_sha256_exact':final_count,'final_original_arrays_exact':final_arrays,
       'final_only':args.final_only,
       'predictions_only_from_restored_G':args.predictions_only,
       'max_array_difference':maxdiff,'evaluation_labels_read':False,
       'numpy':np.__version__,'scipy':scipy.__version__,'python':sys.version,
       'seconds':time.monotonic()-begun,'per_arm':results}
    output=args.report or dest/'INNER_CACHE_REBUILD_RESULT.json'
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k!='per_arm'},indent=2),flush=True)

if __name__=='__main__':main()
