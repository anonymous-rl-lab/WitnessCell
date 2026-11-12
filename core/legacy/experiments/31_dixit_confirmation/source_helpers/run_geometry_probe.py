#!/usr/bin/env python3
"""Graph-only selection audit plus synthetic kernel-transmission diagnosis."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse,hashlib,importlib.util,json,sys,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np

def label_seed(ids,pairs,outer,tag):
    s='geometry_v1|'+str(outer)+'|'+tag+'|'+'|'.join(sorted(pairs[i] for i in ids))
    return int.from_bytes(hashlib.sha256(s.encode()).digest()[:8],'little')

def graph_folds(ids,pairs,outer,x,mode='cycle3'):
    """Deterministic graph-only split. Public API; no response argument exists."""
    ids=np.asarray(ids,int)
    if mode=='balanced2':
        order=sorted(ids,key=lambda i:hashlib.sha256(('balanced2_v1|'+str(outer)+'|'+pairs[i]).encode()).hexdigest())
        return [np.array(sorted(order[j::2]),int) for j in range(2)]
    if mode!='cycle3':raise ValueError(mode)
    rng=np.random.default_rng(label_seed(ids,pairs,outer,mode));allset=set(ids)
    rankcache={}
    def rank(a):
        key=tuple(sorted(a))
        if key not in rankcache:rankcache[key]=int(np.linalg.matrix_rank(x[list(key)]))
        return rankcache[key]
    def objective(ff):
        train=[sorted(allset-set(z)) for z in ff]
        trainrank=[rank(t) for t in train]
        deg=[int(np.min(x[t].sum(axis=0))) for t in train]
        valrank=[rank(z) for z in ff]
        df=[len(z)-r for z,r in zip(ff,valrank)]
        return (min(trainrank),min(deg),min(df),sum(df),-sum(valrank))
    best=None;bestscore=None
    for restart in range(6):
        order=rng.permutation(ids);ff=[list(map(int,order[j::3])) for j in range(3)];score=objective(ff)
        for it in range(500):
            a,b=rng.choice(3,2,replace=False);i=int(rng.integers(len(ff[a])));j=int(rng.integers(len(ff[b])))
            ff[a][i],ff[b][j]=ff[b][j],ff[a][i];candidate=objective(ff)
            if candidate>=score:score=candidate
            else:ff[a][i],ff[b][j]=ff[b][j],ff[a][i]
        if bestscore is None or score>bestscore:bestscore=score;best=[z.copy() for z in ff]
    return [np.array(sorted(z),int) for z in best]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--asset-root',type=Path,required=True);ap.add_argument('--bridge-root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--tasks',default='Dixit,Wessels55');args=ap.parse_args()
    spec=importlib.util.spec_from_file_location('geometry_bridge',args.bridge_root/'run_bridge.py');b=importlib.util.module_from_spec(spec);sys.modules[spec.name]=b;spec.loader.exec_module(b)
    root=args.asset_root;out=args.out;out.mkdir(parents=True,exist_ok=True)
    if (out/'RUN_STARTED.json').exists():raise FileExistsError(out)
    gp=root/'experiment_assets/WitnessCell_Partner_Component_Run_20261018/model/global_pair_model.py'
    gm=b.module('global_pair_model',gp);old=b.module('legacy_rna_controls',root/'source_v19/WitnessCell_Paper_v19_Full/source_data/RNA_controls/run_controls.py')
    new=b.module('frozen_rna_confirmation',root/'work_B/run_frozen_rna_confirmation.py');rr=b.module('readout_control',root/'work_A/readout_control.py')
    science=SimpleNamespace(tune_baseline=gm.tune_baseline,baseline_predict=gm.baseline_predict,projected_blocks=gm.projected_blocks,sat=gm.sat,training_readout_control=rr.training_readout_control)
    b.dump(out/'RUN_STARTED.json',{'utc':b.stamp(),'plan_sha256':b.sha(Path(__file__).with_name('PLAN.txt')),'script_sha256':b.sha(__file__),
        'outer_scoring_after_all_predictions_frozen':True,'stage':'retrospective development diagnostics'})
    tick=time.monotonic();allprepared=[];geom=[];decisions=[];transmission=[]
    for taskname in args.tasks.split(','):
        task=b.load_task(taskname,root,gm,new,old);n,p=task.d.y.shape;lut={a:i for i,a in enumerate(task.pairs)}
        arrays={k:np.full((n,p),np.nan) for k in ['legacy_M0','new_zero',*[f'new_{mode}_{obj}' for mode in ('hash3','cycle3','balanced2') for obj in ('P','full')]]}
        foldids=np.full(n,-1,int)
        for fold in task.folds:
            f=fold['fold'];tr=np.array([lut[a] for a in fold['train_pairs']],int);q=np.array([lut[a] for a in fold['test_pairs']],int)
            vals={c:task.values[c].copy() for c in ['control',*task.nodes,*fold['train_pairs']]}
            halves={c:tuple(z.copy() for z in task.halves[c]) for c in ['control',*fold['train_pairs']]}
            model=new.MeasuredSystem(task.pairs,task.nodes,vals,halves,f,science)
            with np.load(args.bridge_root/'run1'/taskname/f'fold_{f}/PREDICTION_COMPONENTS.npz',allow_pickle=False) as z:comp={k:z[k] for k in z.files}
            if not np.array_equal(comp['query'],q):raise ValueError('query mismatch')
            arrays['legacy_M0'][q]=comp['legacy_M0'];arrays['new_zero'][q]=comp['new_M0'];foldids[q]=f
            for mode in ('hash3','cycle3','balanced2'):
                ff=new.folds(tr,task.pairs,f,'selection') if mode=='hash3' else graph_folds(tr,task.pairs,f,model.d.x,mode)
                score={obj:{lam:0. for lam in b.LAMBDAS} for obj in ('P','full','half1_P','half2_P')};zero={obj:0. for obj in score};totaldf=0;energies=np.zeros(3)
                for j,held in enumerate(ff):
                    fit=np.array(sorted(set(tr)-set(held)),int);base,g,w,pars,audit=model.component(fit,held,f'inner{j}')
                    pv=b.projection(model.d.x[held]);df=len(held)-int(np.linalg.matrix_rank(model.d.x[held]));totaldf+=df
                    truth=model.d.y[held];h1=np.stack([model.halves[task.pairs[i]][0]-model.halves['control'][0] for i in held]);h2=np.stack([model.halves[task.pairs[i]][1]-model.halves['control'][1] for i in held])
                    sat=gm.sat(model.d.additive[held],pars['alpha']);z1=pv@(h1-sat);z2=pv@(h2-sat)
                    energies+=np.array([np.sum(z1*z2),np.sum(z1*z1),np.sum(z2*z2)])
                    for obj,target in [('P',truth),('full',truth),('half1_P',h1),('half2_P',h2)]:
                        err=base-target;err=err if obj=='full' else pv@err;zero[obj]+=float(np.sum(err*err))
                        for lam in b.LAMBDAS:
                            er=base+g[lam]*w.weight-target;er=er if obj=='full' else pv@er;score[obj][lam]+=float(np.sum(er*er))
                    geom.append({'task':taskname,'outer':f,'mode':mode,'inner':j,'train_pairs':len(fit),'validation_pairs':len(held),'train_rank':int(np.linalg.matrix_rank(model.d.x[fit])),
                        'train_min_endpoint_degree':int(np.min(model.d.x[fit].sum(axis=0))),'validation_contrast_dimensions':df,
                        'validation_nodes':int(np.count_nonzero(model.d.x[held].sum(axis=0)))})
                choices={obj:min(b.LAMBDAS,key=lambda z:(score[obj][z],-z)) for obj in score}
                for obj in ('P','full'):arrays[f'new_{mode}_{obj}'][q]=comp['new_M0']+comp['G_OOF'][b.LAMBDAS.index(choices[obj])]*comp['weights']
                record={'task':taskname,'outer':f,'mode':mode,'choices':choices,'validation_contrast_dimensions':totaldf,'scores':score,'zero':zero,
                    'half_signal_cosine':float(energies[0]/np.sqrt(energies[1]*energies[2])) if energies[1]*energies[2]>0 else None,
                    'P_training_gain':1-score['P'][choices['P']]/zero['P'],'full_training_gain':1-score['full'][choices['full']]/zero['full'],
                    'half_choices_agree':choices['half1_P']==choices['half2_P'],'splits':[[task.pairs[i] for i in h] for h in ff]}
                decisions.append(record);b.dump(out/taskname/f'outer_{f}_{mode}_decision.json',record)
            if taskname=='Dixit':
                e=np.stack([task.values[a]-task.values['control'] for a in task.nodes]);gram=e@e.T/e.shape[1]
                edge=np.array([[task.nodes.index(a) for a in pair.split('+')] for pair in task.pairs]);pq=b.projection(model.d.x[q]);ptr=b.projection(model.d.x[tr])
                rng=np.random.default_rng(202610200 +f);syn={mode:[] for mode in ('identity','empirical','planted_aligned')}
                operator_norms={}
                for sim in range(100):
                    u=rng.normal(size=len(task.nodes));signal=u[edge[:,0]]*u[edge[:,1]];target=pq@signal[q];den=float(np.sum(target*target))
                    for mode,eg in [('identity',np.eye(len(task.nodes))),('empirical',gram),('planted_aligned',np.outer(u,u))]:
                        k=gm.symmetric_kernel(edge,eg);kt,kq,_,_=gm.projected_blocks(model.d.x[tr],model.d.x[q],k[np.ix_(tr,tr)],k[np.ix_(q,tr)],k[np.ix_(q,q)])
                        pred=kq@np.linalg.solve(kt+.1*np.eye(len(tr)),signal[tr]);pg=pq@pred
                        syn[mode].append(1-float(np.sum((pg-target)**2))/den)
                        if sim==0:operator_norms[mode]={'projected_cross_frobenius':float(np.linalg.norm(pq@kq@ptr)),
                            'projected_train_frobenius':float(np.linalg.norm(kt)),'kernel_diagonal_fraction':float(np.sum(np.diag(k)**2)/np.sum(k*k))}
                transmission.append({'outer':f,'lambda_fixed':.1,'draws':100,'operator':operator_norms,'synthetic_gain':{mode:{'mean':float(np.mean(a)),
                    'median':float(np.median(a)),'q05':float(np.quantile(a,.05)),'q95':float(np.quantile(a,.95)),'min':float(np.min(a)),'max':float(np.max(a))} for mode,a in syn.items()},
                    'interpretation':'Synthetic partner signal independent of empirical singles; aligned representation is known in this diagnostic and is not an empirical learned model.'})
            b.emit('geometry_fold_frozen',task=taskname,outer=f,seconds=time.monotonic()-tick)
        b.dump(out/taskname/'PREDICTIONS_FROZEN.json',{'utc':b.stamp(),'arrays':{k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in arrays.items()},'no_outer_scoring_yet':True})
        allprepared.append((task,arrays,foldids))
    b.dump(out/'ALL_PREDICTIONS_FROZEN.json',{'utc':b.stamp(),'all_tasks_frozen':True})
    results=[]
    for task,arrays,foldids in allprepared:
        score=b.evaluate(task.d,arrays,foldids,task.primary_blocks,gm);b.dump(out/task.name/'SCORES.json',score)
        results.extend({'task':task.name,'method':k,**v} for k,v in score.items())
    b.csv_write(out/'RESULTS.csv',results);b.csv_write(out/'GEOMETRY.csv',geom);b.dump(out/'DECISIONS.json',decisions);b.dump(out/'SYNTHETIC_TRANSMISSION.json',transmission)
    b.dump(out/'COMPLETE.json',{'utc':b.stamp(),'seconds':time.monotonic()-tick,'original_outer_components_unchanged':True,'outer_results_used_for_selection':False})
    b.emit('geometry_complete',seconds=time.monotonic()-tick)

if __name__=='__main__':main()
