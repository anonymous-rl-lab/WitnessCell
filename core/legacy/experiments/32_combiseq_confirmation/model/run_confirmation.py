#!/usr/bin/env python3
"""One frozen confirmation; prepare reads labels only, run reads outcomes."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse,csv,gzip,hashlib,importlib.util,json,sys,time
from datetime import datetime,timezone
from pathlib import Path
from types import SimpleNamespace
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
NOISES=(.01,.03,.1,.3,1.)
METHODS=('additive','scalar_saturation','endpoint_control','fixed_reference','system_fixed_reference')

def stamp():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def module(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def projection(x):
    p=np.eye(len(x))-x@np.linalg.pinv(x);return (p+p.T)/2
def csv_write(p,rows):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def cohort():
    with gzip.open(ROOT/'combi_metadata/GSE174695_large_scale_meta.csv.gz','rt') as f:rows=list(csv.DictReader(f))
    allowed=[];excluded=[]
    for row in rows:
        a,b=row['Drug1'],row['Drug2'];sample=row['']
        if '_ctrlReadsPerGene' in sample:excluded.append({'sample':sample,'reason':'extra control identity not matched'});continue
        if a==b and a!='DMSO':excluded.append({'sample':sample,'reason':'same compound two-channel dose'});continue
        cond='control' if a==b=='DMSO' else (b if a=='DMSO' else a if b=='DMSO' else '+'.join(sorted((a,b))))
        allowed.append({'sample':sample,'run':row['Run'],'condition':cond,'drug1':a,'drug2':b})
    nodes=sorted({r['condition'] for r in allowed if r['condition']!='control' and '+' not in r['condition']})
    pairs=sorted({r['condition'] for r in allowed if '+' in r['condition']})
    assert all(set(p.split('+'))<=set(nodes) for p in pairs)
    conditions=['control',*nodes,*pairs];runs=sorted({r['run'] for r in allowed})
    assert runs==['rep1','rep2','rep3']
    assert all(any(r['condition']==c and r['run']==t for r in allowed) for c in conditions for t in runs)
    x=np.zeros((len(pairs),len(nodes)))
    for i,p in enumerate(pairs):
        for g in p.split('+'):x[i,nodes.index(g)]=1
    return allowed,excluded,nodes,pairs,conditions,runs,x

def prepare():
    if (ROOT/'FREEZE.json').exists():raise FileExistsError('Already frozen')
    allowed,excluded,nodes,pairs,conditions,runs,x=cohort()
    geo=module('confirm_graph',ROOT/'model/frozen_graph_folds.py')
    order=sorted(range(len(pairs)),key=lambda i:hashlib.sha256(('CombiSeq_confirm_v1|'+pairs[i]).encode()).hexdigest())
    folds=[]
    for f in range(3):
        q=np.array(sorted(order[f::3]));tr=np.array(sorted(set(range(len(pairs)))-set(q)))
        inner=geo.graph_folds(tr,pairs,f,x,'cycle3')
        folds.append({'fold':f,'train':tr.tolist(),'query':q.tolist(),'inner':[z.tolist() for z in inner],
            'train_rank':int(np.linalg.matrix_rank(x[tr])),'query_rank':int(np.linalg.matrix_rank(x[q])),
            'query_contrast_dimensions':len(q)-int(np.linalg.matrix_rank(x[q]))})
    genes=(ROOT/'combi_metadata/matrix_gene_names_only.txt').read_text().splitlines()
    mapping=json.loads((ROOT/'model/FROZEN_PROGRAM_ID_MAPPING.json').read_text())
    ids=[r['ensembl'] for r in mapping['rows'] if r['ensembl'] in genes and genes.count(r['ensembl'])==1]
    assert len(ids)>=40
    members={'source':'GSE174695','nodes':nodes,'pairs':pairs,'conditions':conditions,'runs':runs,
        'allowed':allowed,'excluded':excluded,'program_ids':ids,'folds':folds,
        'n_pairs':len(pairs),'n_nodes':len(nodes),'n_program':len(ids),'full_graph_rank':int(np.linalg.matrix_rank(x)),
        'expression_outcomes_examined':False}
    dump(ROOT/'COHORT_AND_SPLITS.json',members)
    print(json.dumps({k:members[k] for k in ('n_pairs','n_nodes','n_program','full_graph_rank')},indent=2),flush=True)
    print('Outer contrast dimensions:',[f['query_contrast_dimensions'] for f in folds],flush=True)

def freeze():
    if (ROOT/'FREEZE.json').exists():raise FileExistsError('Already frozen')
    files=[ROOT/'CONFIRMATION_PROTOCOL.txt',ROOT/'COHORT_AND_SPLITS.json',*sorted((ROOT/'model').glob('*.py')),
        *sorted((ROOT/'model').glob('*.json')),ROOT/'combi_metadata/GSE174695_large_scale_meta.csv.gz',
        ROOT/'combi_metadata/GSE174695_large_scale_matrix.csv.gz',ROOT/'combi_metadata/matrix_gene_names_only.txt']
    for p in files:
        if not p.is_file():raise FileNotFoundError(p)
    dump(ROOT/'FREEZE.json',{'utc':stamp(),'status':'FROZEN BEFORE EXPRESSION PROCESSING',
        'expression_outcomes_examined':False,'files':{str(p.relative_to(ROOT)):sha(p) for p in files},
        'new_source':'Combi-seq / GSE174695','confirmation_scope':'cross-intervention-type RNA, refitted frozen procedure'})
    print('FREEZE',sha(ROOT/'FREEZE.json'),flush=True)

def load_task(members,out):
    """Source-specific normalization; no all-condition feature selection."""
    import pandas as pd
    p=ROOT/'combi_metadata/GSE174695_large_scale_matrix.csv.gz'
    frame=pd.read_csv(p,index_col=0)
    expected=(ROOT/'combi_metadata/matrix_gene_names_only.txt').read_text().splitlines()
    assert frame.index.astype(str).tolist()==expected
    counts=frame.to_numpy(dtype=np.float64,copy=True)
    if not np.all(np.isfinite(counts)) or counts.min()<0:raise ValueError('Invalid raw counts')
    cols={str(c):i for i,c in enumerate(frame.columns)}
    sums=counts.sum(axis=0)
    if np.any(sums<=0):raise ValueError('Zero-total sample: frozen integrity gate')
    logged=np.log1p(counts*(10000./sums)[None,:]);del counts
    conditions=members['conditions'];runs=members['runs'];records=members['allowed']
    byrun=np.empty((len(runs),len(conditions),len(expected)))
    for t,run in enumerate(runs):
        for i,c in enumerate(conditions):
            group=[z for z in records if z['condition']==c and z['run']==run]
            orientations=sorted({(z['drug1'],z['drug2']) for z in group})
            by_orientation=[]
            for a,b in orientations:
                samples=[cols[z['sample']] for z in group if z['drug1']==a and z['drug2']==b]
                by_orientation.append(logged[:,samples].mean(axis=1))
            byrun[t,i]=np.mean(by_orientation,axis=0)
    means=byrun.mean(axis=0);values={c:means[i].copy() for i,c in enumerate(conditions)}
    nodes=members['nodes'];pairs=members['pairs'];ci={c:i for i,c in enumerate(conditions)}
    x=np.zeros((len(pairs),len(nodes)))
    for i,pair in enumerate(pairs):
        for g in pair.split('+'):x[i,nodes.index(g)]=1
    y=np.stack([values[p]-values['control'] for p in pairs])
    pi=np.array([expected.index(g) for g in members['program_ids']],int)
    source_summary={'utc':stamp(),'n_gene_rows':len(expected),'n_raw_samples':len(frame.columns),
        'count_file_sha256':sha(p),'normalization':'log1p(10000*reads/sample_total)',
        'min_sample_reads':float(sums.min()),'median_sample_reads':float(np.median(sums)),
        'max_sample_reads':float(sums.max()),'gene_selection':'none; fixed program lookup only'}
    dump(out/'SOURCE_PROCESSING.json',source_summary)
    np.savez_compressed(out/'CONDITION_MEANS.npz',means=means,byrun=byrun,conditions=np.array(conditions),genes=np.array(expected),pi=pi)
    return SimpleNamespace(nodes=nodes,pairs=pairs,values=values,d=SimpleNamespace(y=y,x=x),
        byrun=byrun,ci=ci,pi=pi,genes=expected)

def endpoint_raw(state,qp,noise):
    xq,_=state.design(qp);xt=state.xt
    return (xq@xt.T/2)@np.linalg.solve(xt@xt.T/2+noise*np.eye(len(xt)),state.resid)

def fit(task,fold,core,gm):
    train=np.array(fold['train'],int);program=task.pi
    b=SimpleNamespace(projection=projection)
    lambdas=core.LAMBDAS;scores={lam:0. for lam in lambdas}
    cross={n:0. for n in NOISES};norm={n:0. for n in NOISES};zero=0.;details=[]
    for j,validation in enumerate(fold['inner']):
        va=np.array(validation,int);tr=np.array(sorted(set(train)-set(va)),int)
        state=core.FixedFit(task,tr,program,gm,b);qp=[task.pairs[i] for i in va]
        target=state.h(task.d.y[va]);base=state.base(qp);resid=target-base;zero+=float(np.sum(resid**2))
        for lam in lambdas:scores[lam]+=float(np.sum((state.predict(qp,lam,'global')-target)**2))
        for noise in NOISES:
            raw=endpoint_raw(state,qp,noise);cross[noise]+=float(np.sum(raw*resid));norm[noise]+=float(np.sum(raw**2))
        details.append({'inner':j,'fit':tr.tolist(),'validation':va.tolist(),'alpha':state.alpha,
            'variance_floor':state.floor,'allowed_normalization':list(state.allowed)})
    choice=min(lambdas,key=lambda z:(scores[z],-z))
    gamma={n:float(np.clip(cross[n]/norm[n],0.,1.)) if norm[n]>1e-12 else 0. for n in NOISES}
    end_scores={n:zero-2*gamma[n]*cross[n]+gamma[n]**2*norm[n] for n in NOISES}
    noise=min(NOISES,key=lambda n:(end_scores[n],n))
    state=core.FixedFit(task,train,program,gm,b)
    info={'lambda':choice,'noise':noise,'gamma':gamma[noise],'alpha':state.alpha,
        'partner_inner_scores':scores,'endpoint_inner_scores':end_scores,'endpoint_gammas':gamma,
        'scalar_saturation_inner_sse':zero,'inner':details}
    return state,info

def predictions(state,qp,info):
    endpoint=state.base(qp)+info['gamma']*endpoint_raw(state,qp,info['noise'])
    return {'additive':state.additive(qp).ravel(),'scalar_saturation':state.base(qp).ravel(),
        'endpoint_control':endpoint.ravel(),
        'fixed_reference':state.predict(qp,info['lambda'],'global').ravel(),
        'system_fixed_reference':(endpoint+state.correction(qp,info['lambda'],'global')).ravel()}

def run(out):
    fr=json.loads((ROOT/'FREEZE.json').read_text())
    for rel,h in fr['files'].items():
        if sha(ROOT/rel)!=h:raise ValueError('Frozen file changed: '+rel)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if (out/'EXPRESSION_ACCESS_START.json').exists():raise FileExistsError(out)
    dump(out/'EXPRESSION_ACCESS_START.json',{'utc':stamp(),'freeze_sha256':sha(ROOT/'FREEZE.json'),
        'first_expression_processing_in_this_protocol':True})
    tick=time.monotonic();members=json.loads((ROOT/'COHORT_AND_SPLITS.json').read_text())
    task=load_task(members,out)
    core=module('confirm_core',ROOT/'model/frozen_reference_core.py');gm=module('confirm_gm',ROOT/'model/frozen_legacy_baseline.py')
    n=len(task.pairs);pp={k:np.full(n,np.nan) for k in METHODS};states=[];foldids=np.full(n,-1,int)
    for fold in members['folds']:
        f=fold['fold'];q=np.array(fold['query'],int);qp=[task.pairs[i] for i in q]
        state,info=fit(task,fold,core,gm);pred=predictions(state,qp,info)
        for k in METHODS:pp[k][q]=pred[k]
        foldids[q]=f;state.save(out/f'FIT_{f}.npz');dump(out/f'FIT_{f}.json',info)
        states.append((fold,state,info));print(json.dumps({'stage':'outer_prediction_frozen','fold':f,'seconds':time.monotonic()-tick}),flush=True)
    np.savez_compressed(out/'PREDICTIONS.npz',pairs=np.array(task.pairs),foldids=foldids,**pp)
    dump(out/'ALL_PREDICTIONS_FROZEN.json',{'utc':stamp(),'prediction_sha256':sha(out/'PREDICTIONS.npz'),'heldout_losses_computed':False})
    audit=[]
    for fold,state,info in states:
        f=fold['fold'];q=np.array(fold['query'],int);qp=[task.pairs[i] for i in q];ref=pp['system_fixed_reference'][q]
        singleton=np.array([predictions(state,[p],info)['system_fixed_reference'].item() for p in qp])
        reversed_pred=predictions(state,['+'.join(p.split('+')[::-1]) for p in qp],info)['system_fixed_reference']
        chunks=np.concatenate([predictions(state,qp[j:j+7],info)['system_fixed_reference'] for j in range(0,len(qp),7)])
        shuffled=predictions(state,qp[::-1],info)['system_fixed_reference'][::-1]
        duplicated=predictions(state,qp+[qp[0]],info)['system_fixed_reference'][:-1]
        invariant=max(float(np.max(abs(ref-z))) for z in (singleton,reversed_pred,chunks,shuffled,duplicated))
        original_y=task.d.y[q].copy();original_v={task.pairs[i]:task.values[task.pairs[i]].copy() for i in q}
        rng=np.random.default_rng(91823+f);task.d.y[q]=rng.normal(0,1e8,original_y.shape)
        for i in q:task.values[task.pairs[i]]=rng.normal(0,1e8,len(task.genes))
        poisoned,info2=fit(task,fold,core,gm);pred2=predictions(poisoned,qp,info2)
        task.d.y[q]=original_y;task.values.update(original_v)
        diff=max(float(np.max(abs(pp[k][q]-pred2[k]))) for k in METHODS)
        audit.append({'fold':f,'query_invariance_max_abs':invariant,'poison_prediction_max_abs':diff,
            'poison_all_choices_identical':info==info2,'PASS':invariant<=1e-10 and diff==0 and info==info2})
    dump(out/'INVARIANCE_AND_LEAKAGE.json',audit)
    if not all(v['PASS'] for v in audit):
        dump(out/'COMPLETE.json',{'utc':stamp(),'decision':'INVALID_RUN','reason':'query invariance or leakage audit failed','heldout_losses_computed':False})
        raise RuntimeError('Frozen implementation audit failed; no scientific score issued')
    truth=np.full(n,np.nan);runtruth=np.full((3,n),np.nan);rows=[];foldrows=[];squares={k:np.zeros(n) for k in METHODS};psquares={k:np.zeros(n) for k in METHODS}
    for fold,state,info in states:
        f=fold['fold'];q=np.array(fold['query'],int);pq=projection(task.d.x[q]);truth[q]=state.h(task.d.y[q]).ravel()
        for t in range(3):
            dy=np.stack([task.byrun[t,task.ci[task.pairs[i]]]-task.byrun[t,task.ci['control']] for i in q])
            runtruth[t,q]=state.h(dy).ravel()
        for k in METHODS:
            err=pp[k][q]-truth[q];squares[k][q]=err**2;psquares[k][q]=(pq@err)**2
            foldrows.append({'fold':f,'method':k,'n_pairs':len(q),'full_sse':float(np.sum(err**2)),
                'P_sse':float(np.sum((pq@err)**2)),'P_dimensions':int(np.linalg.matrix_rank(pq,tol=1e-8))})
    for k in METHODS:
        rows.append({'method':k,'n_pairs':n,'full_sse':float(squares[k].sum()),'P_sse':float(psquares[k].sum()),
            'full_gain_vs_endpoint':float(1-squares[k].sum()/squares['endpoint_control'].sum()),
            'P_gain_vs_endpoint':float(1-psquares[k].sum()/psquares['endpoint_control'].sum()),
            'full_gain_vs_saturation':float(1-squares[k].sum()/squares['scalar_saturation'].sum()),
            'P_gain_vs_saturation':float(1-psquares[k].sum()/psquares['scalar_saturation'].sum())})
    runrows=[]
    for t in range(3):
        errors={k:float(np.sum((pp[k]-runtruth[t])**2)) for k in METHODS}
        for k in METHODS:runrows.append({'run':members['runs'][t],'method':k,'full_sse':errors[k],
            'gain_vs_endpoint':1-errors[k]/errors['endpoint_control'],'gain_vs_saturation':1-errors[k]/errors['scalar_saturation']})
    edges=np.array([[task.nodes.index(g) for g in p.split('+')] for p in task.pairs])
    rng=np.random.default_rng(20261020);wg=rng.exponential(1.,(4000,len(task.nodes)));ww=wg[:,edges[:,0]]*wg[:,edges[:,1]]
    interval={}
    for metric,sq in [('full',squares),('P',psquares)]:
        for control in ('endpoint_control','scalar_saturation'):
            gains=1-(ww@sq['system_fixed_reference'])/(ww@sq[control]);interval[metric+'_vs_'+control]=dict(zip(('q025','median','q975'),map(float,np.quantile(gains,[.025,.5,.975]))))
    deletions=[]
    for g in task.nodes:
        keep=np.array([g not in p.split('+') for p in task.pairs]);deletions.append({'omitted_compound':g,'n_pairs':int(keep.sum()),
            'full_gain_vs_endpoint':float(1-squares['system_fixed_reference'][keep].sum()/squares['endpoint_control'][keep].sum())})
    np.savez_compressed(out/'EVALUATION.npz',truth=truth,runtruth=runtruth,x=task.d.x,pairs=np.array(task.pairs),foldids=foldids,
        **{f'full_loss_{k}':v for k,v in squares.items()},**{f'P_loss_{k}':v for k,v in psquares.items()})
    csv_write(out/'RESULTS.csv',rows);csv_write(out/'FOLD_RESULTS.csv',foldrows);csv_write(out/'RUN_TARGET_RESULTS.csv',runrows);csv_write(out/'ENDPOINT_DELETION.csv',deletions)
    dump(out/'CONDITIONAL_STABILITY.json',{'interval_type':'conditional endpoint multiplier stability range; fixed OOF predictions; not unconditional biological CI','draws':4000,'seed':20261020,'results':interval})
    main=next(r for r in rows if r['method']=='system_fixed_reference');runpositive=sum(r['gain_vs_endpoint']>0 for r in runrows if r['method']=='system_fixed_reference')
    checks={'principal_gain_ge_2_percent':main['full_gain_vs_endpoint']>=.02,
        'conditional_lower_bound_positive':interval['full_vs_endpoint_control']['q025']>0,
        'full_gain_vs_saturation_positive':main['full_gain_vs_saturation']>0,
        'partner_gain_vs_endpoint_positive':main['P_gain_vs_endpoint']>0,
        'at_least_two_runs_positive':runpositive>=2,'invariance_and_leakage':all(v['PASS'] for v in audit)}
    decision='STRONG_SUPPORT' if all(checks.values()) else 'NOT_CONFIRMED' if main['full_gain_vs_endpoint']<=0 else 'PARTIAL_OR_INCONCLUSIVE'
    dump(out/'COMPLETE.json',{'utc':stamp(),'seconds':time.monotonic()-tick,'decision':decision,'checks':checks,'primary':main,
        'new_genetic_double_RNA_confirmation':False,'cross_intervention_RNA_confirmation':True,
        'post_score_scientific_tuning':False,'freeze_sha256':sha(ROOT/'FREEZE.json')})
    print(json.dumps({'decision':decision,'checks':checks,'primary':main,'stability':interval,'seconds':time.monotonic()-tick},indent=2),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','freeze','run']);a.add_argument('--out',type=Path);args=a.parse_args()
    if args.stage=='prepare':prepare()
    elif args.stage=='freeze':freeze()
    else:run(args.out or ROOT/'run1')
