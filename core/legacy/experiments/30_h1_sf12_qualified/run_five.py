"""Fixed five-rectangle extension of old12 and case_R. No new tuning.
Each task has its own sealed factory; cross-task training reuse is allowed.
This is historical-panel development evaluation, not independent confirmation.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,time,datetime,argparse,platform
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source'))
import numpy as np
from fold_factory import FoldFactory,projection,incidence,PARAMS
from legacy_model import old12,solve,RHO,OLD_C
from case_features import feature
OLD=ROOT.parent/'case_transfer_20261026'
ASSET=ROOT.parent/'h1_probe_20261026/training_inputs/WitnessCell_Neural_Partner_Pilot_20261021/joung_assets/ready/TRAINING_BATCH_ONLY.npz'
HALVES=OLD/'data/TRAINING_RAW_HALVES.npz'
KNOWLEDGE=OLD/'data/KNOWLEDGE_NODE.npz'
OUT=ROOT/'run'
NAMES=['M0','old12','case_R']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False))
def factory(query_ids=()):
    with np.load(KNOWLEDGE,allow_pickle=False) as z:k=z['knowledge_node'].copy()
    f=FoldFactory(asset_path=ASSET,knowledge_node=k)
    with np.load(HALVES,allow_pickle=False) as z:
        assert np.array_equal(z['pairs'],f.pairs) and np.array_equal(z['genes'],f.genes)
        f._raw_halves=z['halves'].copy()
    q=np.array(query_ids,dtype=int)
    f.asset['batch_values'][np.isin(f.asset['obs_condition'],len(f.nodes)+q)]=np.nan
    f._raw_halves[:,q]=np.nan
    return f
def task_specs():
    orig=json.loads((ROOT/'source/ORIGINAL_SPLIT.json').read_text())
    rect=orig['all_primary_candidates'];tasks=[]
    for row in rect:
        if not row['eligible']:continue
        q=row['query_ids'];tr=sorted(set(range(21))-set(q));eps=[]
        for ep in rect:
            va=ep['query_ids']
            if set(va)<=set(tr):eps.append(dict(train_ids=sorted(set(tr)-set(va)),val_ids=va,val_pairs=ep['query_pairs']))
        tasks.append(dict(name=f'Q{len(tasks)+1}',query_ids=q,query_pairs=row['query_pairs'],train_ids=tr,episodes=eps,
            minimum_query_Npi=row['minimum_query_precision_proxy'],minimum_qualified_degree=row['minimum_qualified_degree_after_holdout'],
            known_previous_result=q==orig['primary_query_ids']))
    assert len(tasks)==5 and tasks[3]['known_previous_result']
    assert [len(t['episodes']) for t in tasks]==[4,4,2,4,4]
    return tasks,orig['all_qualified_cases']
def freeze():
    tasks,qualified=task_specs();f=factory();c=np.array([.5,-.5,-.5,.5]);fullc=np.zeros((5,21))
    for i,t in enumerate(tasks):
        assert np.max(abs(c@incidence(f.edges(t['query_ids']),10)))<1e-12
        fullc[i,t['query_ids']]=c
    spec=dict(scope='Historical21 fixed-model extension; Q4 reproduction, other4 development evaluation.',tasks=tasks,
        qualified_cases=qualified,contrast=c.tolist(),methods=NAMES,M0_parameters=PARAMS,
        old12=dict(rho=RHO,c=OLD_C,refit='every outer/inner support graph'),
        case_R=dict(rho=RHO,parameters=1,postfit_amplitude=1.,endpoint_repair=False,
            loss='sum inner whole-rectangle projected SSE to Y−M0−old12; normalized ridge',
            episodes='all qualified rectangles inside support, including missing anchors',
            reliability='raw-composition technical halves, all support edges; NOT fixed7 calibrated halves'),
        selection='No hyperparameter/model/task selection from outer performance.',
        workflow='Freeze all5 predictions before outer scoring; seal each own outer labels. Cross-task training reuse allowed.',
        aggregation=dict(partner='1−sum task partnerSSE(model)/sum task partnerSSE(reference)',
            full_edge_balanced='Average each edge loss across tasks then average8edges; not prediction ensemble',
            diagnostics=['all5','exclude_Q4','per_task','median_positive_count','leave_one_task_out'],
            uncertainty='Descriptive. No independent-rectangle/RNA bootstrap; overlapping labels and fits.'),
        proceed_conditions=['pooled old12 partner gain>0','excludeQ4 partner gain>0','edge-balanced full gain>=0'],
        contrast_identity_rank=int(np.linalg.matrix_rank(fullc)),unique_query_edges=int(np.count_nonzero(np.any(fullc!=0,axis=0))),
        input_hashes={str(p):sha(p) for p in [ASSET,HALVES,KNOWLEDGE,ROOT/'source/ORIGINAL_SPLIT.json']},
        code_hashes={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'run_five.py',ROOT/'score_five.py',*(ROOT/'source').glob('*.py')]},
        compute=dict(device='CPU',threads=1,optimizer='exactsolve',epochs=None,seeds=[],smoke_limit_seconds=600,small_limit_seconds=3600),
        environment=dict(python=platform.python_version(),numpy=np.__version__))
    p=OUT/'PROTOCOL.json'
    if p.exists():assert json.loads(p.read_text())['spec']==spec,'Frozen specification changed'
    else:save(p,dict(frozen_utc=stamp(),spec=spec))
    return spec
def check_frozen():
    d=json.loads((OUT/'PROTOCOL.json').read_text())['spec']
    for n,h in d['code_hashes'].items():assert sha(ROOT/n)==h,n
    for p,h in d['input_hashes'].items():assert sha(p)==h,p
    return d
def design(f,fac,qualified):
    h,trace=feature(f.d['single_centered_exact'],f.d['train_edges'],f.R,f.d['eval_edges'],
        np.isin(f.d['train_pairs'],qualified),knowledge_node=f.d['knowledge_node'],singles_centered=True)
    w=fac.reliability(f.d['train_ids'])
    return h[:,:,:1]*w[None,:,None],dict(cases=trace,R_mean=float(w.mean()),R_nonzero=float(np.mean(w>0)))
def compute(t,spec,poison_counts=False):
    fac=factory(t['query_ids'])
    if poison_counts:fac.asset['counts'][len(fac.nodes)+np.array(t['query_ids'])]*=7919
    M=np.zeros((1,1));b=np.zeros(1);archive={};traces=[];nfits=0
    for i,ep in enumerate(t['episodes']):
        tr,va=ep['train_ids'],ep['val_ids'];assert not set(tr+va)&set(t['query_ids'])
        f=fac.get_fold(tr,va);nfits+=1
        d12,meta=old12(f);h,trace=design(f,fac,spec['qualified_cases'])
        target=fac.fixed_reference(va)-f.d['eval_base']-d12
        p=projection(f.d['eval_edges'],f.n);pr=p@target;ph=np.einsum('ab,bgp->agp',p,h)
        mi=np.einsum('egp,egq->pq',ph,ph);bi=np.einsum('egp,eg->p',ph,pr)
        M+=mi;b+=bi;archive[f'e{i}_M']=mi;archive[f'e{i}_b']=bi
        degrees=[len(a['case_rows']) for q in trace['cases']['query_support'] for a in q['anchors']]
        traces.append(dict(train_ids=tr,val_ids=va,old12=meta,case_R=trace,
            minimum_anchor_cases=min(degrees),target_partner_SSE=float(np.sum(pr*pr))))
    assert traces
    theta,info=solve(M,b);assert info['normal_error']<1e-10
    f=fac.get_fold(t['train_ids'],t['query_ids']);nfits+=1
    d12,meta=old12(f);h,trace=design(f,fac,spec['qualified_cases'])
    m0=f.d['eval_base'];old=m0+d12;pred=np.stack([m0,old,old+np.einsum('egp,p->eg',h,theta)])
    assert np.isfinite(pred).all();archive.update(M=M,b=b,theta=theta)
    return pred,h,archive,dict(task=t,case_R_fit=info,case_R_theta=theta.tolist(),episodes=traces,
        outer_old12=meta,outer_cases=trace,baseline_fits=nfits,own_outer_labels_masked=True,outer_scoring_truth_accessed=False)
def smoke():
    start=time.monotonic();spec=check_frozen();t=spec['tasks'][3];pred,h,stats,trace=compute(t,spec)
    with np.load(OLD/'run/FROZEN_PREDICTIONS.npz',allow_pickle=False) as z:
        expected=z['predictions'][[list(z['names']).index(n) for n in NAMES]]
    error=float(np.max(abs(pred-expected)))
    prev=json.loads((OLD/'run/FITS.json').read_text())['case_R']['theta'];te=float(np.max(abs(stats['theta']-prev)))
    poisoned,_,_,_=compute(t,spec,True);pe=float(np.max(abs(pred-poisoned)))
    assert error<1e-12 and te<1e-12 and pe==0 and time.monotonic()-start<600
    save(OUT/'SMOKE.json',dict(status='PASS',Q4_prediction_max_abs=error,Q4_theta_max_abs=te,
        withheld_counts_poison_max_abs=pe,query_expression_and_halves_masked=True,seconds=time.monotonic()-start,
        scope='KnownQ4 reproduction; no newQ1/Q2/Q3/Q5 outer scoring.'))
    print(json.dumps({'stage':'smoke','status':'PASS','seconds':time.monotonic()-start,'Q4_error':error}),flush=True)
def predict():
    start=time.monotonic();spec=check_frozen();assert json.loads((OUT/'SMOKE.json').read_text())['status']=='PASS'
    assert not (OUT/'ALL_PREDICTIONS_FROZEN.json').exists();files={};totalfits=0
    for t in spec['tasks']:
        t0=time.monotonic();pred,h,stats,trace=compute(t,spec);folder=OUT/t['name'];folder.mkdir(exist_ok=True)
        np.savez_compressed(folder/'PREDICTIONS.npz',predictions=pred,h_case_R=h,names=np.array(NAMES),
            query_ids=np.array(t['query_ids']),contrast=np.array(spec['contrast']))
        np.savez_compressed(folder/'TRAIN_STATS.npz',**stats)
        trace['seconds']=time.monotonic()-t0;save(folder/'TRACE.json',trace)
        for n in ['PREDICTIONS.npz','TRAIN_STATS.npz','TRACE.json']:files[str((folder/n).relative_to(OUT))]=sha(folder/n)
        totalfits+=trace['baseline_fits']
        print(json.dumps({'stage':'predict','task':t['name'],'episodes':len(t['episodes']),'seconds':trace['seconds'],'scored':False}),flush=True)
        assert time.monotonic()-start<3600
    save(OUT/'ALL_PREDICTIONS_FROZEN.json',dict(frozen_utc=stamp(),files=files,protocol_sha256=sha(OUT/'PROTOCOL.json'),
        seconds=time.monotonic()-start,baseline_fits=totalfits,old12_fits=totalfits,case_R_fits=5,
        own_outer_scoring_not_performed=True,cross_task_training_label_reuse=True))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['freeze','smoke','predict']);a=ap.parse_args();globals()[a.stage]()
