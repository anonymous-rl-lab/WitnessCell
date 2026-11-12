"""Frozen, quality-screened empirical-case development experiment.

No new hyperparameter search: normalized ridge rho=.1, no post-fit amplitude.
Every coefficient-training episode hides a WHOLE qualified rectangle and
rebuilds measurements, original M0, shared12 and empirical case features.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import json,hashlib,time,datetime,argparse,platform
import numpy as np
from fold_factory import FoldFactory,incidence,projection
from case_features import feature,smoke as feature_smoke

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'run'
RHO=.1
OLD_C=.6717359152620326
SEEDS=[20261026021,20261026022,20261026023]
ARMS=['case_scalar','case_R','case_transfer']+[f'wrong_{i}' for i in range(3)]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False))
def loadsplit():return json.loads((ROOT/'audit/PRIMARY_SPLIT.json').read_text())
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sealed_factory(split):
    f=FoldFactory()
    q=np.array(split['primary_query_ids'])
    rows=np.isin(f.asset['obs_condition'],len(f.nodes)+q)
    f.asset['batch_values'][rows]=np.nan
    # Prevent later auxiliary reliability from ever decoding a usable primary
    # target. Only support-side rows are used in all subsequent computations.
    with np.load(ROOT/'data/TRAINING_RAW_HALVES.npz',allow_pickle=False) as z:
        f._raw_halves=z['halves'].copy()
    f._raw_halves[:,q]=np.nan
    return f

def old12(f):
    s=f.d['single_centered_exact'];scale=np.sqrt(np.mean(s*s,axis=0))
    floor=max(1e-10,float(np.median(scale[scale>1e-12]))*1e-3)
    scale=np.maximum(scale,floor);t=np.tanh(s/scale)
    x,_=f.node_features(2,'knowledge');tr=f.d['train_edges'];qe=f.d['eval_edges']
    prod=x[tr[:,0]]*x[tr[:,1]];center=prod.mean(0)
    ps=np.maximum(np.sqrt(np.mean((prod-center)**2,axis=0)),1e-8)
    def make(edges):
        a,b=t[edges[:,0]],t[edges[:,1]]
        ap,an=np.maximum(a,0),np.maximum(-a,0)
        bp,bn=np.maximum(b,0),np.maximum(-b,0)
        op=ap*bn+an*bp
        loc=np.stack([ap*bp,an*bn,op,op*np.tanh(a+b)],axis=-1)*scale[None,:,None]
        mod=np.column_stack([np.ones(len(edges)),(x[edges[:,0]]*x[edges[:,1]]-center)/ps])
        return (loc[:,:,:,None]*mod[:,None,None,:]).reshape(len(edges),f.G,12)
    ft,fq=make(tr),make(qe)
    d=np.einsum('ab,bgp->agp',f.P,ft,optimize=True)
    e=fq-np.einsum('ab,bgp->agp',incidence(qe,f.n)@f.Pinv,ft,optimize=True)
    m=np.einsum('egp,egq->pq',d,d,optimize=True)
    b=np.einsum('egp,eg->p',d,f.P@f.R,optimize=True)
    th,meta=solve(m,b)
    return np.einsum('egp,p->eg',e,th*OLD_C),dict(theta=th.tolist(),**meta)

def solve(m,b):
    m=(m+m.T)/2
    sc=np.sqrt(np.maximum(np.diag(m),0));sc=np.maximum(sc,max(float(sc.max(initial=0))*1e-10,1e-15))
    n=m/sc[:,None]/sc[None,:];rhs=b/sc
    w=np.linalg.solve(n+RHO*np.eye(len(b)),rhs);theta=w/sc
    normal=float(np.linalg.norm((n+RHO*np.eye(len(b)))@w-rhs)/max(np.linalg.norm(rhs),1e-30))
    ev=np.linalg.eigvalsh(n)
    return theta,dict(rho=RHO,normal_error=normal,rank=int(np.sum(ev>max(float(ev.max(initial=0)),1.)*1e-10)),
                     normalized_eigenvalues=ev.tolist(),effective_df=float(np.trace(np.linalg.solve(n+RHO*np.eye(len(b)),n))))

def designs(f,fac,split):
    eligible=np.isin(f.d['train_pairs'],split['all_qualified_cases'])
    args=(f.d['single_centered_exact'],f.d['train_edges'],f.R,f.d['eval_edges'],eligible)
    kw=dict(knowledge_node=f.d['knowledge_node'],singles_centered=True)
    h,meta=feature(*args,**kw)
    rw=fac.reliability(f.d['train_ids'])
    out={'case_scalar':h[:,:,:1],'case_R':h[:,:,:1]*rw[None,:,None],'case_transfer':h}
    traces={'correct':meta,'R_weights_mean':float(rw.mean()),'R_weights_nonzero_fraction':float(np.mean(rw>0))}
    for i,seed in enumerate(SEEDS):
        out[f'wrong_{i}'],traces[f'wrong_{i}']=feature(*args,perm_seed=seed,**kw)
    return out,traces

def protocol():
    OUT.mkdir(exist_ok=True);sp=loadsplit();f=FoldFactory()
    spec=dict(scope='Quality-screened development holdout within historically used21; NOT independent confirmation.',
        primary_query_pairs=sp['primary_query_pairs'],primary_train_pairs=sp['primary_train_pairs'],
        data_qualification=sp['qualification_rule'],rho=RHO,postfit_amplitude=1.,
        old12=dict(rho=.1,c=OLD_C,coefficients='refit per support graph',status='historical development hyperparameters fixed'),
        arms=ARMS,permutation_seeds=SEEDS,main_arm='case_transfer',
        case_definition='raw empirical Y−M0; includes baseline error and noise, not pure biological interaction',
        weights='exp(.5*cos(panelcentered single B,C)+.5*cos(fixedGenePT64 B,C)), temperature1',
        multipliers=['1','d=T_C−T_B','T_A*d','(T_C+T_B)*d','T_A*(T_C+T_B)*d'],
        training_episodes=[x['val_pairs'] for x in sp['all_nested_candidates']],
        episode_rule='all4 qualified rectangles, even missing-anchor cases; support=T minus Q, full pipeline rebuilt',
        training_loss='sum_Q ||P_Q(Y_Q−M0_S−old12_S−H_Q theta)||² + diagonal normalized ridge .1',
        selection='NONE. Only1fullysupported validation rectangle; nesting leaves no qualified rectangle training. Fixedcapacitytest, no productionupgrade.',
        missing_anchor='that side0, other sideweight.5; no rescale',
        R_control='raw-composition technical halves, not calibrated7library repeats',
        input_hashes={str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else p.name:sha(p) for p in [ROOT/'audit/PRIMARY_SPLIT.json',f.asset_path,ROOT/'data/TRAINING_RAW_HALVES.npz']},
        code_hashes={n:sha(ROOT/n) for n in ['run_case.py','fold_factory.py','case_features.py']},
        compute=dict(device='CPU',threads=1,epochs='exactsolve',fits=6,smoke_limit_seconds=600,small_limit_seconds=3600),
        environment=dict(python=platform.python_version(),numpy=np.__version__))
    path=OUT/'PROTOCOL.json'
    if path.exists():assert json.loads(path.read_text())['spec']==spec,'Frozen specification changed'
    else:save(path,dict(frozen_utc=stamp(),spec=spec))
    return sp

def smoke():
    start=time.monotonic();sp=protocol();fac=sealed_factory(sp);tr=sp['primary_train_ids'];qe=sp['primary_query_ids']
    f=fac.get_fold(tr,qe);base,oldmeta=old12(f);ds,traces=designs(f,fac,sp)
    checks=feature_smoke()
    # Coefficients can be recovered in a known feature-space problem. This is
    # an interface/solver check, not a claim about true biological signals.
    h=ds['case_transfer'];p=projection(f.d['eval_edges'],f.n);d=np.einsum('ab,bgp->agp',p,h)
    known=np.array([.12,-.08,.04,.02,-.01]);target=np.einsum('egp,p->eg',d,known)
    m=np.einsum('egp,egq->pq',d,d);b=np.einsum('egp,eg->p',d,target)
    recovered=np.linalg.lstsq(d.reshape(-1,5),target.ravel(),rcond=1e-12)[0]
    err=float(np.max(abs(d.reshape(-1,5)@recovered-target.ravel())))
    zero,_=solve(m,np.zeros(5));assert err<1e-10 and np.max(abs(zero))==0
    # M0 and old12 reproduction uses archived full21 features without loading
    # archived query truth. It is separate from the primary sealed fitting.
    full=FoldFactory().get_fold(np.arange(21),['CDX1+HNF4A','CDX1+MSGN1','HNF4A+NHLH1','MSGN1+NHLH1'])
    od,_=old12(full)
    archive=ROOT.parent/'h1_probe_20261026/probe_run/FROZEN_PREDICTIONS.npz'
    with np.load(archive,allow_pickle=False) as z:
        i=list(z['method_names']).index('old12');old_error=float(np.max(abs(od-z['deltas'][i])))
    assert old_error<1e-12
    save(OUT/'SMOKE.json',dict(status='PASS',feature_checks=checks,known_feature_reconstruction_error=err,
        old12_archive_max_abs=old_error,primary_query_labels_unavailable_globally=False,
        sealed_candidate_branch_primary_labels_unavailable=True,
        historical_reproduction_access='Separate full21 old12 reproduction includes primary pairs as historical training rows; no resulting predictions or coefficients flow into candidate fit, selection or primary prediction.',
        primary_design_rank=int(np.linalg.matrix_rank(m)),
        primary_support=traces['correct'],seconds=time.monotonic()-start))
    print(json.dumps({'stage':'smoke','status':'PASS','seconds':time.monotonic()-start,'old12_error':old_error}),flush=True)

def train():
    start=time.monotonic();sp=protocol();assert (OUT/'SMOKE.json').exists();assert not (OUT/'FIT_FREEZE.json').exists()
    fac=sealed_factory(sp);stats={a:dict(M=np.zeros((1 if a in ARMS[:2] else 5,)*2),b=np.zeros(1 if a in ARMS[:2] else 5)) for a in ARMS}
    episodes=[];archive={}
    for i,ep in enumerate(sp['all_nested_candidates']):
        tr,va=ep['train_ids'],ep['val_ids'];assert not set(tr+va)&set(sp['primary_query_ids'])
        f=fac.get_fold(tr,va);d12,info12=old12(f);y=fac.fixed_reference(va)
        r=y-f.d['eval_base']-d12;p=projection(f.d['eval_edges'],f.n);r=p@r
        dd,trace=designs(f,fac,sp);row=dict(index=i,train_ids=tr,val_ids=va,trace=trace,old12=info12,methods={})
        for arm,h in dd.items():
            h=np.einsum('ab,bgp->agp',p,h,optimize=True)
            m=np.einsum('egp,egq->pq',h,h,optimize=True);b=np.einsum('egp,eg->p',h,r,optimize=True)
            stats[arm]['M']+=m;stats[arm]['b']+=b
            archive[f'e{i}_{arm}_M']=m;archive[f'e{i}_{arm}_b']=b
            row['methods'][arm]=dict(target_SSE=float(np.sum(r*r)),feature_energy=np.diag(m).tolist())
        episodes.append(row)
        print(json.dumps({'stage':'training_episodes','episode':i+1,'of':len(sp['all_nested_candidates']),'seconds':time.monotonic()-start}),flush=True)
        assert time.monotonic()-start<3600
    fits={}
    for arm,s in stats.items():
        theta,info=solve(s['M'],s['b']);assert info['normal_error']<1e-9
        fits[arm]=dict(theta=theta.tolist(),**info)
        archive[arm+'_M']=s['M'];archive[arm+'_b']=s['b']
    np.savez_compressed(OUT/'TRAIN_STATS.npz',**archive)
    save(OUT/'EPISODES.json',episodes);save(OUT/'FITS.json',fits)
    save(OUT/'FIT_FREEZE.json',dict(time_utc=stamp(),primary_truth_read=False,seconds=time.monotonic()-start,
        files={n:sha(OUT/n) for n in ['PROTOCOL.json','SMOKE.json','TRAIN_STATS.npz','FITS.json','EPISODES.json']}))
    print(json.dumps({'stage':'train','seconds':time.monotonic()-start,'fits':fits}),flush=True)

def predict():
    start=time.monotonic();sp=protocol();freeze=json.loads((OUT/'FIT_FREEZE.json').read_text())
    for n,digest in freeze['files'].items():assert sha(OUT/n)==digest
    assert not (OUT/'PREDICTION_FREEZE.json').exists()
    fac=sealed_factory(sp);f=fac.get_fold(sp['primary_train_ids'],sp['primary_query_ids'])
    d12,meta12=old12(f);dd,trace=designs(f,fac,sp);fits=json.loads((OUT/'FITS.json').read_text())
    preds={'M0':f.d['eval_base'],'old12':f.d['eval_base']+d12};hfiles={}
    for arm,h in dd.items():
        delta=np.einsum('egp,p->eg',h,np.array(fits[arm]['theta']),optimize=True)
        preds[arm]=preds['old12']+delta;hfiles[arm]=h
    np.savez_compressed(OUT/'FROZEN_PREDICTIONS.npz',names=np.array(list(preds)),predictions=np.stack(list(preds.values())),
        pairs=np.array(sp['primary_query_pairs']),genes=fac.genes,contrast=np.array(sp['primary_contrast']))
    np.savez_compressed(OUT/'QUERY_FEATURES.npz',**hfiles)
    save(OUT/'PREDICT_TRACE.json',dict(support=trace,old12=meta12))
    save(OUT/'PREDICTION_FREEZE.json',dict(time_utc=stamp(),primary_truth_read=False,seconds=time.monotonic()-start,
        files={n:sha(OUT/n) for n in ['FITS.json','PROTOCOL.json','FROZEN_PREDICTIONS.npz','QUERY_FEATURES.npz','PREDICT_TRACE.json']}))
    print(json.dumps({'stage':'predict','seconds':time.monotonic()-start,'primary_truth_read':False}),flush=True)

def score():
    sp=loadsplit();freeze=json.loads((OUT/'PREDICTION_FREEZE.json').read_text())
    for n,digest in freeze['files'].items():assert sha(OUT/n)==digest
    fac=FoldFactory();y=fac.fixed_reference(sp['primary_query_ids'])
    with np.load(OUT/'FROZEN_PREDICTIONS.npz',allow_pickle=False) as z:
        names=z['names'];pred=z['predictions'];v=z['contrast'];genes=z['genes']
    p=np.outer(v,v);m0=pred[0];old=pred[1];r0=y-m0;ro=y-old
    denomF=float(np.sum(r0*r0));denomC=float(np.sum((p@r0)**2));oldF=float(np.sum(ro*ro));oldC=float(np.sum((p@ro)**2))
    scores={};qstats={}
    with np.load(OUT/'QUERY_FEATURES.npz',allow_pickle=False) as z:
        for arm in ARMS:
            h=z[arm];ph=np.einsum('ab,bgp->agp',p,h)
            qstats[arm+'_MF']=np.einsum('egp,egq->pq',h,h)
            qstats[arm+'_bF']=np.einsum('egp,eg->p',h,ro)
            qstats[arm+'_MC']=np.einsum('egp,egq->pq',ph,ph)
            qstats[arm+'_bC']=np.einsum('egp,eg->p',ph,p@ro)
    for name,yp in zip(names,pred):
        e=y-yp;delta=yp-old;pd=p@delta;pr=p@ro
        sf=float(np.sum(e*e));sc=float(np.sum((p@e)**2));de=float(np.sum(pd*pd));dot=float(np.sum(pr*pd))
        gainF=oldF-sf;gainC=oldC-sc
        scores[str(name)]=dict(full_SSE=sf,partner_SSE=sc,full_gain_vs_M0_pct=100*(1-sf/denomF),
            partner_gain_vs_M0_pct=100*(1-sc/denomC),full_gain_vs_old12_pct=100*(1-sf/oldF),
            partner_gain_vs_old12_pct=100*(1-sc/oldC),case_direction_cosine=dot/np.sqrt(de*oldC) if de>1e-30 else None,
            case_delta_partner_energy=de,case_dot=dot,full_SSE_reduction_vs_old12=gainF,
            partner_SSE_reduction_vs_old12=gainC,endpoint_SSE_reduction_vs_old12=gainF-gainC,
            per_query_full_gain_vs_old12_pct=(100*(1-np.sum(e*e,axis=1)/np.sum(ro*ro,axis=1))).tolist())
    qstats.update(old_full_SSE=np.array(oldF),old_partner_SSE=np.array(oldC),m0_full_SSE=np.array(denomF),m0_partner_SSE=np.array(denomC))
    np.savez_compressed(OUT/'EVALUATION_STATS.npz',**qstats)
    np.savez_compressed(OUT/'REFERENCE.npz',y=y,genes=genes,pairs=np.array(sp['primary_query_pairs']),contrast=v)
    save(OUT/'SCORES.json',dict(scope='fixed-rule quality-screened development test; one partner contrast',methods=scores))
    print(json.dumps({'stage':'score','methods':scores},ensure_ascii=False),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','smoke','train','predict','score'])
    args=parser.parse_args()
    {'freeze':protocol,'smoke':smoke,'train':train,'predict':predict,'score':score}[args.stage]()
