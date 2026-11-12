#!/usr/bin/env python3
"""Read-only legacy/new RNA bridge and controlled residual diagnostics."""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import argparse,csv,datetime,hashlib,importlib.util,itertools,json,sys,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np

HERE=Path(__file__).resolve().parent
DEFAULT_ROOT=Path('/workspace/scratch/7517238e1d04/witnesscell_v20_20261019')
LAMBDAS=(.1,1.,10.)
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def emit(event,**v):print(json.dumps({'utc':stamp(),'event':event,**v}),flush=True)
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m
def canonical(p):return '+'.join(sorted(p.split('+')))
def projection(x):
    p=np.eye(len(x))-x@np.linalg.pinv(x);return (p+p.T)/2
def simple_parameters(p):return {k:p[k] for k in ('alpha','noise','gamma')}
def csv_write(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def rectangles(pairs,allowed):
    lookup={p:i for i,p in enumerate(pairs)};allowed=set(map(int,allowed))
    nodes=sorted({g for i in allowed for g in pairs[i].split('+')});found=set()
    for a,b,c,d in itertools.combinations(nodes,4):
        for excluded in (((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))):
            ex={canonical('+'.join(p)) for p in excluded}
            ps=[canonical('+'.join(p)) for p in itertools.combinations((a,b,c,d),2) if canonical('+'.join(p)) not in ex]
            if all(p in lookup and lookup[p] in allowed for p in ps):found.add(tuple(sorted(lookup[p] for p in ps)))
    return [list(b) for b in sorted(found)]

def load_task(name,root,gm,new,old):
    assets=root/'experiment_assets'
    if name=='Wessels55':
        cp=assets/'repo/experiments/22_metric_calibration_stress/assets/formal_contracts/Wessels.phase_m_contract.npz'
        dp=assets/'WitnessCell_Partner_Component_Run_20261018/model/frozen_inputs/B_FROZEN_DESIGN.json'
        design=json.loads(dp.read_text())
        with np.load(cp,allow_pickle=False) as z:s={k:z[k] for k in z.files}
        ci={str(c):i for i,c in enumerate(s['conditions'])}
        targets=sorted({g for c in ci if '+' in c for g in c.split('+')})
        keep=~np.isin(s['genes'].astype(str),targets);genes=s['genes'].astype(str)[keep]
        pairs=design['complete_RNA_single_pairs'];nodes=sorted(design['allowed_singles']);folds=design['outer_folds']
        variance=s['full_variances'];primary=design['model_panel_geometry']['picked']
        source_paths=[cp,dp];symbols=genes.copy()
    else:
        cp=root/'work_B/Dixit_external_ready_compact.npz';dp=root/'work_B/Dixit_labels_only_threefold_candidate.csv'
        with np.load(cp,allow_pickle=False) as z:s={k:z[k] for k in z.files}
        ci={str(c):i for i,c in enumerate(s['conditions'])};nodes=sorted(c for c in ci if c!='control' and '+' not in c)
        keep=~np.isin(s['gene_symbols'].astype(str),nodes);genes=s['genes'].astype(str)[keep];symbols=s['gene_symbols'].astype(str)[keep]
        with dp.open() as f:rows=list(csv.DictReader(f))
        pairs=[canonical(r['gene_a']+'+'+r['gene_b']) for r in rows];foldids=np.array([int(r['fold']) for r in rows])
        eval_blocks=[b for f in range(3) for b in rectangles(pairs,np.flatnonzero(foldids==f))]
        primary=[[pairs[i] for i in b] for b in eval_blocks]
        used=set();training_blocks=[]
        for b in sorted(eval_blocks,key=lambda ids:gm.digest([pairs[i] for i in ids])):
            if not used.intersection(b):training_blocks.append([pairs[i] for i in b]);used.update(b)
        folds=[]
        for f in range(3):
            tr=np.flatnonzero(foldids!=f);q=np.flatnonzero(foldids==f)
            inn=gm.preserved_folds(tr,pairs,training_blocks)
            folds.append({'fold':f,'train_pairs':[pairs[i] for i in tr],'test_pairs':[pairs[i] for i in q],
                'inner_test_pairs':[[pairs[i] for i in ids] for ids in inn],
                'test_rectangle_blocks':[[pairs[i] for i in b] for b in rectangles(pairs,q)]})
        design={'outer_folds':folds,'allowed_singles':nodes,'model_panel_geometry':{'picked':training_blocks},
            'legacy_inner_fold_adapter':'greedy edge-disjoint subset of original within-outer rectangles; preserved_folds'}
        variance=s['variances'];source_paths=[cp,dp]
    ctrl=s['full_means'][ci['control'],keep]
    values={c:s['full_means'][ci[c],keep].copy() for c in ['control',*nodes,*pairs]}
    half_values={c:(s['first_means'][ci[c],keep].copy(),s['second_means'][ci[c],keep].copy()) for c in ['control',*pairs]}
    pi=np.array([ci[p] for p in pairs]);single={g:values[g]-ctrl for g in nodes}
    d=gm.make_dataset(pairs,single,np.stack([values[p]-ctrl for p in pairs]),design['model_panel_geometry']['picked'])
    inp=old.Inputs(d,genes,ctrl,s['first_means'][pi][:,keep],s['second_means'][pi][:,keep],s['full_means'][pi][:,keep],variance[pi][:,keep],design)
    return SimpleNamespace(name=name,d=d,inp=inp,values=values,halves=half_values,genes=genes,symbols=symbols,
        pairs=pairs,nodes=nodes,folds=folds,primary_blocks=primary,source_paths=source_paths,design=design)

def paired_components(model,tr,q,stage,gm,new):
    base,goof,w,pars,audit=model.component(tr,q,stage)
    residual_in=model.d.y[tr]-gm.baseline_predict(model.d,tr,tr,pars)
    xt,xq=model.d.x[tr],model.d.x[q];k=model.d.pair_kernel
    kt,kq,_,_=gm.projected_blocks(xt,xq,k[np.ix_(tr,tr)],k[np.ix_(q,tr)],k[np.ix_(q,q)])
    eig,v=np.linalg.eigh(kt);left=kq@v
    gin={lam:((left/(eig+lam))@v.T)@residual_in for lam in LAMBDAS}
    return {'base':base,'G':{'in_sample':gin,'OOF':goof},'w':w,'pars':pars,'audit':audit}

def residual_diagnostic(model,tr,q,component,gm,new,task):
    d=model.d;pars=component['pars'];fitted=gm.baseline_predict(d,tr,tr,pars)
    satfull=gm.sat(d.additive[tr],pars['alpha']);endfull=fitted-satfull
    pos={i:j for j,i in enumerate(tr)};oof=np.empty_like(fitted);satoof=np.empty_like(fitted);records=[]
    for f,held in enumerate(new.folds(tr,model.pairs,model.outer,'final:residual')):
        fit=np.array([i for i in tr if i not in set(held)],int)
        value,pp=model.core(fit,held,'final:residual'+f':aux{f}')
        loc=[pos[i] for i in held];oof[loc]=value;satoof[loc]=gm.sat(d.additive[held],pp['alpha'])
        records.append({'fold':f,'train_pairs':[model.pairs[i] for i in fit],'query_pairs':[model.pairs[i] for i in held],
                        'parameters':simple_parameters(pp)})
    p=projection(d.x[tr]);rin=d.y[tr]-fitted;ro=d.y[tr]-oof
    terms={'projected_in_sample_residual':p@rin,'projected_OOF_residual':p@ro,
           'projected_residual_difference':p@(ro-rin),'projected_saturation_difference':p@(satfull-satoof),
           'projected_endpoint_difference':p@(endfull-(oof-satoof))}
    energy={k:float(np.sum(v*v)) for k,v in terms.items()};e=energy['projected_in_sample_residual']
    a=terms['projected_residual_difference'];b=terms['projected_in_sample_residual']
    den=float(np.linalg.norm(a)*np.linalg.norm(b));cross=float(np.sum(terms['projected_saturation_difference']*terms['projected_endpoint_difference']))
    prediction=[];pq=projection(d.x[q])
    for lam in LAMBDAS:
        gi=component['G']['in_sample'][lam];go=component['G']['OOF'][lam];delta=go-gi
        prediction.append({'lambda':lam,'G_in_energy':float(np.sum(gi*gi)),'G_OOF_energy':float(np.sum(go*go)),
            'G_difference_energy':float(np.sum(delta*delta)),'projected_query_G_difference_energy':float(np.sum((pq@delta)**2)),
            'difference_to_in_energy_ratio':float(np.sum(delta*delta)/np.sum(gi*gi)) if np.sum(gi*gi)>0 else None})
    return {'task':task,'outer_fold':model.outer,'training_pairs':[model.pairs[i] for i in tr],
        'test_pairs':[model.pairs[i] for i in q],'full_core_parameters':simple_parameters(pars),'OOF_core_blocks':records,
        'energies':energy,'difference_to_in_energy_ratio':energy['projected_residual_difference']/e if e>0 else None,
        'difference_to_in_norm_ratio':float(np.sqrt(energy['projected_residual_difference']/e)) if e>0 else None,
        'difference_vs_in_cosine':float(np.sum(a*b)/den) if den>0 else None,
        'saturation_endpoint_cross_product':cross,
        'decomposition_max_absolute_difference':float(np.max(np.abs(a-terms['projected_saturation_difference']-terms['projected_endpoint_difference']))),
        'full_fit_endpoint_projection_energy':float(np.sum((p@endfull)**2)),
        'G_differences':prediction,'interpretation':'Training-side diagnostic, not a test-result attribution or biological effect.'}

def isolated(model,tr,q,gm,new,task):
    inner=[]
    for f,held in enumerate(new.folds(tr,model.pairs,model.outer,'selection')):
        fit=np.array([i for i in tr if i not in set(held)],int)
        c=paired_components(model,fit,held,f'inner{f}',gm,new);inner.append((held,c))
    outer=paired_components(model,tr,q,'final',gm,new)
    choices=[]
    for residual in ('in_sample','OOF'):
        scores={metric:{lam:0. for lam in LAMBDAS} for metric in ('full','P','C')};counts={'full':0,'P':0,'C':0}
        for held,c in inner:
            truth=model.d.y[held];x=model.d.x[held];p=projection(x);df=len(held)-int(np.linalg.matrix_rank(x))
            blocks=rectangles(model.pairs,held);loc={j:i for i,j in enumerate(held)}
            counts['full']+=truth.size;counts['P']+=df*truth.shape[1];counts['C']+=len(blocks)*truth.shape[1]
            for lam in LAMBDAS:
                error=c['base']+c['G'][residual][lam]-truth
                scores['full'][lam]+=float(np.sum(error**2))
                if df>0:scores['P'][lam]+=float(np.sum((p@error)**2))
                for block in blocks:
                    _,sign=gm.block_contrast(model.d.y,block,model.d.x)
                    scores['C'][lam]+=float(np.sum((sign@error[[loc[i] for i in block]])**2))
        for objective in ('full','P','C'):
            used=objective if counts[objective]>0 else 'full'
            lam=min(LAMBDAS,key=lambda z:(scores[used][z],-z))
            choices.append({'residual':residual,'requested_objective':objective,'used_objective':used,
                'contrast_objective_available':counts[objective]>0,'lambda':lam,
                'scores':{str(k):v/counts[used] for k,v in scores[used].items()},'counts':counts})
    return outer,choices,residual_diagnostic(model,tr,q,outer,gm,new,task)

def evaluate(d,pred,foldids,blocks,gm):
    out={};blockids=[[d.pairs.index(p) for p in b] for b in blocks]
    for method,value in pred.items():
        err=value-d.y;full=float(np.mean(err**2));psse=0.;df=0
        for f in sorted(set(foldids)):
            q=np.flatnonzero(foldids==f);x=d.x[q];p=projection(x)
            psse+=float(np.sum((p@err[q])**2));df+=len(q)-int(np.linalg.matrix_rank(x))
        cs=[]
        for b in blockids:
            assert len(set(foldids[b]))==1
            _,sign=gm.block_contrast(d.y,b,d.x);cs.append(float(np.mean((sign@err[b])**2)))
        out[method]={'full_mse':full,'projected_sse':psse,'C_mse':float(np.mean(cs)) if cs else None,
                     'test_contrast_dimensions':df,'rectangles':len(cs)}
    for method,v in out.items():
        core='legacy_M0' if method.startswith('legacy_') else 'new_zero'
        for metric in ('full_mse','projected_sse','C_mse'):
            if v[metric] is not None:
                v[metric+'_gain_vs_own_core']=1-v[metric]/out[core][metric]
                v[metric+'_gain_vs_legacy_core']=1-v[metric]/out['legacy_M0'][metric]
    return out

def run_task(task,out,gm,new,old,science,root,tick):
    out.mkdir(parents=True,exist_ok=True);d=task.d;n,p=d.y.shape;lut={q:i for i,q in enumerate(d.pairs)}
    arrays={};foldids=np.full(n,-1,int);choices_all=[];diagnostics=[];metadata=[]
    def fill(name,q,v):
        if name not in arrays:arrays[name]=np.full((n,p),np.nan)
        arrays[name][q]=v
    for fold in task.folds:
        f=fold['fold'];start=time.monotonic();tr=np.array([lut[x] for x in fold['train_pairs']],int);q=np.array([lut[x] for x in fold['test_pairs']],int)
        fi=out/f'fold_{f}';fi.mkdir(exist_ok=True)
        # Legacy calls and unchanged old readout formula. No permutation is fitted.
        lr=old.fit_outer_controls(task.inp,fold,include_permutations=False)
        for method,value in lr['predictions'].items():fill('legacy_'+method,q,value)
        allowed=['control',*task.nodes,*fold['train_pairs']]
        vals={c:task.values[c].copy() for c in allowed}
        halves={c:tuple(x.copy() for x in task.halves[c]) for c in ['control',*fold['train_pairs']]}
        model=new.MeasuredSystem(task.pairs,task.nodes,vals,halves,f,science)
        nr,g,w,nmeta=model.fit_outer(tr,q)
        for method,value in nr.items():fill('new_'+method,q,value)
        emit('bridge_fold_predictions_ready',task=task.name,fold=f,seconds=time.monotonic()-start,
             total_seconds=time.monotonic()-tick,aggregate_scored=False)
        outer,choices,diag=isolated(model,tr,q,gm,new,task.name)
        for residual in ('in_sample','OOF'):
            for lam in LAMBDAS:
                gg=outer['G'][residual][lam]
                fill(f'fixed_{residual}_lambda{lam:g}_unit',q,outer['base']+gg)
                fill(f'fixed_{residual}_lambda{lam:g}_R',q,outer['base']+gg*outer['w'].weight)
        for c in choices:
            gg=outer['G'][c['residual']][c['lambda']]
            for control in ('unit','R'):
                pred=outer['base']+gg*(outer['w'].weight if control=='R' else 1.)
                fill(f"tuned_{c['residual']}_{c['requested_objective']}_{control}",q,pred)
        foldids[q]=f;choices_all.append({'fold':f,'choices':choices});diagnostics.append(diag)
        dump(fi/'RESIDUAL_DIAGNOSTICS.json',diag);dump(fi/'OBJECTIVE_SELECTION.json',choices)
        dump(fi/'LEGACY_FIT.json',lr['metadata']);dump(fi/'NEW_FIT.json',nmeta);dump(fi/'NEW_CORE_FITS.json',model.core_log)
        # Exact reusable ingredients, not all redundant completed prediction arrays.
        np.savez_compressed(fi/'PREDICTION_COMPONENTS.npz',query=q,legacy_M0=lr['predictions']['M0'],
            legacy_M1=lr['predictions']['M1'],legacy_R=lr['predictions']['R'],legacy_scalar=lr['predictions']['scalar'],
            new_M0=nr['zero'],new_amplitude=nr['amplitude_only'],new_additive=nr['additive'],new_saturation=nr['saturation_only'],
            G_in=np.stack([outer['G']['in_sample'][z] for z in LAMBDAS]),
            G_OOF=np.stack([outer['G']['OOF'][z] for z in LAMBDAS]),weights=w.weight)
        metadata.append({'fold':f,'legacy_lambda':lr['metadata']['core_fit']['lambda'],'new_lambda':nmeta['lambda'],
            'new_policy':nmeta['selected_policy'],'new_uniform_scalar':nmeta['uniform_scalar'],
            'legacy_core_parameters':simple_parameters(lr['metadata']['core_fit']['M0']),
            'new_core_parameters':simple_parameters(nmeta['M0_parameters'])})
        dump(out/'PROGRESS.json',{'utc':stamp(),'task':task.name,'completed_folds':f+1,'elapsed_seconds':time.monotonic()-tick,'aggregate_scored':False})
        emit('residual_probe_fold_complete',task=task.name,fold=f,seconds=time.monotonic()-start,
            residual_difference_to_in_energy_ratio=diag['difference_to_in_energy_ratio'])
        if time.monotonic()-tick>3500:raise TimeoutError('One-hour bound reached; partial unscored records retained')
    if any(not np.isfinite(a).all() for a in arrays.values()):raise ValueError('Incomplete prediction array')
    frozen={'utc':stamp(),'all_outer_predictions_completed_before_scores':True,'metadata':metadata,
        'arrays':{k:{'shape':list(v.shape),'dtype':str(v.dtype),'sha256':hashlib.sha256(v.tobytes()).hexdigest()} for k,v in arrays.items()},
        'components':{str(f.relative_to(out)):sha(f) for f in out.glob('fold_*/PREDICTION_COMPONENTS.npz')}}
    dump(out/'ALL_PREDICTIONS_FROZEN.json',frozen)
    scores=evaluate(d,arrays,foldids,task.primary_blocks,gm)
    replay={}
    if task.name=='Wessels55':
        expected=json.loads((root/'source_v19/WitnessCell_Paper_v19_Full/source_data/RNA_controls/CONTROL_RESULTS.json').read_text())
        for method in ('M0','M1','R','scalar'):
            replay[method]={'full_mse_abs_diff':abs(scores['legacy_'+method]['full_mse']-expected['response_MSE'][method]),
                'C_mse_abs_diff':abs(scores['legacy_'+method]['C_mse']-expected['C_MSE'][method])}
        maximum=max(v for one in replay.values() for v in one.values())
        if maximum>1e-12:raise ValueError('Legacy reproduction mismatch: '+repr(replay))
    else:
        with np.load(root/'work_B/run1/ALL_PREDICTIONS_FROZEN.npz',allow_pickle=False) as z:
            for method in ('zero','unit','uniform','amplitude_only','R_amplitude_noise','selected','additive','saturation_only'):
                replay[method]=float(np.max(np.abs(arrays['new_'+method]-z[method])))
        if max(replay.values())>1e-12:raise ValueError('v20 Dixit reproduction mismatch: '+repr(replay))
    dump(out/'REPLAY_CHECK.json',{'status':'PASS','differences':replay})
    dump(out/'SCORES.json',scores);dump(out/'RESIDUAL_DIAGNOSTICS_ALL.json',diagnostics)
    rows=[{'method':k,**v} for k,v in scores.items()];csv_write(out/'ALL_METHODS.csv',rows)
    bridge=[{'task':task.name,'procedure':proc,'reported_policy':method,
             **scores[method]} for proc,method in [('legacy','legacy_R'),('v20_fixed_R','new_R_amplitude_noise'),('v20_selected','new_selected')]]
    dump(out/'BRIDGE.json',bridge)
    dump(out/'COMPLETE.json',{'utc':stamp(),'status':'all fixed candidates scored; development diagnosis',
        'method_count':len(arrays),'postscore_candidate_selection':False,'elapsed_seconds':time.monotonic()-tick})
    emit('task_complete',task=task.name,replay='PASS',bridge=bridge,total_seconds=time.monotonic()-tick)
    return bridge

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--asset-root',type=Path,default=Path(os.environ.get('WITNESSCELL_ASSET_ROOT',str(DEFAULT_ROOT))))
    ap.add_argument('--out',type=Path,default=HERE/'run1');ap.add_argument('--task',choices=['both','Wessels55','Dixit'],default='both');args=ap.parse_args()
    root=args.asset_root.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'RUN_STARTED.json').exists():raise FileExistsError('Refusing to overwrite prior debug run')
    gp=root/'experiment_assets/WitnessCell_Partner_Component_Run_20261018/model/global_pair_model.py'
    op=root/'source_v19/WitnessCell_Paper_v19_Full/source_data/RNA_controls/run_controls.py'
    npth=root/'work_B/run_frozen_rna_confirmation.py';rp=root/'work_A/readout_control.py'
    gm=module('global_pair_model',gp);old=module('legacy_rna_controls',op);new=module('frozen_rna_confirmation',npth);readout=module('readout_control',rp)
    science=SimpleNamespace(tune_baseline=gm.tune_baseline,baseline_predict=gm.baseline_predict,
        projected_blocks=gm.projected_blocks,sat=gm.sat,training_readout_control=readout.training_readout_control)
    tasks=[load_task(n,root,gm,new,old) for n in (['Wessels55','Dixit'] if args.task=='both' else [args.task])]
    paths=[Path(__file__),HERE/'PROBE_PLAN_v1.txt',gp,op,npth,rp,*[p for task in tasks for p in task.source_paths]]
    dump(out/'RUN_STARTED.json',{'utc':stamp(),'status':'development diagnostic inputs/code frozen before new probe fits',
        'files':{str(p):sha(p) for p in dict.fromkeys(paths)},'known_prior_outcomes':True,'numpy':np.__version__,
        'max_runtime_seconds':3600,'no_postscore_model_choice':True})
    tick=time.monotonic();bridge=[];emit('bridge_started',tasks=[t.name for t in tasks])
    for task in tasks:
        dump(out/(task.name+'_GEOMETRY.json'),{'genes':len(task.genes),'pairs':task.pairs,'nodes':task.nodes,
            'outer_folds':task.folds,'fixed_evaluation_rectangles':task.primary_blocks,
            'legacy_protected_inner_blocks':task.d.blocks,'gene_axis_sha256':hashlib.sha256(task.genes.tobytes()).hexdigest()})
        bridge.extend(run_task(task,out/task.name,gm,new,old,science,root,tick))
    dump(out/'BRIDGE_ALL.json',bridge);csv_write(out/'BRIDGE_ALL.csv',bridge)
    dump(out/'COMPLETE.json',{'utc':stamp(),'seconds':time.monotonic()-tick,'status':'complete development bridge',
        'new_final_model_selected_on_test_results':False})
    emit('bridge_complete',seconds=time.monotonic()-tick)

if __name__=='__main__':main()
