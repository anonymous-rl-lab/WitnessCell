"""R core upgrades selected only within the permitted training conditions."""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'source'
sys.path.insert(0,str(SOURCE/'adapters'))
import rna_interface as ri
from kernel_features import build_kernels,mix_kernels

METHODS=('R_original','R_lambda_matched','R_calibrated','R_12','R_neural','R_plus_guarded','R_plus_partner',
         'R_phase1_primary','R_phase1_calibrated')
CONFIG=dict(version='R-core-gated-20261101-phase2',primary='R_plus_guarded',outer_folds=list(range(5)),inner_folds=3,
    score_genes=4995,lambdas=[.1,1.,10.],eta=[.25,1.],
    families=['native','sf12','neural','joint'],powers=[1,2],candidate_count=42,
    kernel_mixture='KR + eta*Kaux, divided by 1+eta; joint Kaux=(K12+KN)/2; each component independently train-projected/trace-normalized',
    correction='A=(w**power)*(operator @ training_M0_residual), power in {1,2}',
    beta='one nonnegative global scalar; no upper cap; jointly selected with lambda and geometry on pooled inner rectangles',
    guard='pooled inner full RNA SSE no higher than separately nested original R; quadratic feasible interval; numerical relative tolerance only',
    fallback='independent original R with inner-fit-specific nested lambda, not a disguised shared lambda beta=1 candidate',
    R_lambda_matched='native kernel; power=1; beta=1; select lambda by actual final weighted partner SSE; identical to phase1',
    R_calibrated='native geometry with joint lambda/beta and RNA guard; original R fallback',
    R_12='native or sf12 geometry with joint lambda/beta and RNA guard; original R fallback',
    R_neural='native or neural geometry with joint lambda/beta and RNA guard; original R fallback',
    R_plus_partner='all geometries; beta>=0; no RNA guard; predefined secondary tradeoff, never replaces primary',
    tie='original R; power one; fewer extra geometries; smaller eta; larger lambda; beta nearest one; deterministic family order',
    numerical_relative_tolerance=float(64*np.finfo(float).eps),
    outer_query_outcomes_used=False,development_panel=True,
    phase='Post-phase1 development extension motivated by weak-readout deterioration under global amplification; not independent confirmation',
    phase1_reference_columns='R_phase1_primary/R_phase1_calibrated recomputed by selecting only power=1 candidates',
    neural_scope='fixed shared nonlinear interaction representation; not a fitted full-feedback mechanism')

def write_json(path,value):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n');tmp.replace(p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def array_sha(a):
    a=np.ascontiguousarray(a);h=hashlib.sha256();h.update(str(a.dtype).encode());h.update(str(a.shape).encode());h.update(a.tobytes());return h.hexdigest()
def lock():
    names=['train_upgrade.py','kernel_features.py','source/adapters/rna_interface.py',
           'source/adapters/sf12_original/original_code/global_pair_model.py',
           'source/adapters/sf12_original/original_code/legacy_shared_core.py',
           'source/adapters/sf12_original/original_code/twelve_baseline.py']
    record=dict(config=CONFIG,code_sha256={n:sha(ROOT/n) for n in names})
    path=ROOT/'LOCKED_PROTOCOL.json'
    if path.exists():assert json.loads(path.read_text())==record,'Training code/protocol changed after lock'
    else:write_json(path,record)
    return sha(path)
def candidates():
    geom=[('native',0.)]+[(family,eta) for family in ('sf12','neural','joint') for eta in CONFIG['eta']]
    return [dict(name=f'{family}_e{eta:g}_l{lam:g}_p{power}',family=family,eta=eta,lam=lam,power=power)
            for family,eta in geom for lam in CONFIG['lambdas'] for power in CONFIG['powers']]
def moments_load():
    with np.load(SOURCE/'adapters/original_moments.npz',allow_pickle=False) as z:
        return {k:z[k].copy() for k in ('genes','conditions','first_means','second_means')}
def reliability(common,score,train,moments):
    assert np.array_equal(common['genes'],moments['genes'])
    lookup={str(c):i for i,c in enumerate(moments['conditions'])};ids=[lookup[score.pairs[i]] for i in train]
    x=score.x[train];P=np.eye(len(train))-x@np.linalg.pinv(x);mask=common['score_gene_mask']
    a=P@moments['first_means'][ids][:,mask];b=P@moments['second_means'][ids][:,mask]
    num=np.sum(a*b,axis=0);den=.5*np.sum(a*a+b*b,axis=0)
    return np.divide(num,den,out=np.zeros_like(num),where=den>0).clip(0,1)
def clone_fit(score,ids):
    y=np.full_like(score.y,np.nan);y[ids]=score.y[ids]
    return ri.gm.Dataset(score.pairs,score.nodes,y,score.additive,score.x,score.pair_kernel,score.blocks)
def fit_R(common,score,train,query,tune_folds,moments):
    d=clone_fit(score,train);assert np.isnan(d.y[query]).all()
    pred,meta=ri.gm.fit_outer(d,train,query,tune_folds)
    w=reliability(common,score,train,moments)
    return dict(base=pred['M0'],R=pred['M0']+w*(pred['M1']-pred['M0']),w=w,metadata=meta)
def C_matrix(score,ids):
    rows=[];blocks=[];lookup={j:i for i,j in enumerate(ids)}
    for block in ri.gm.available_blocks(score,ids):
        _,sign=ri.gm.block_contrast(np.zeros((len(score.pairs),1)),block,score.x)
        c=np.zeros(len(ids));c[[lookup[i] for i in block]]=sign
        assert np.max(abs(c@score.x[ids]))<1e-12
        rows.append(c);blocks.append([score.pairs[i] for i in block])
    return np.array(rows).reshape(-1,len(ids)),blocks
def operator(bundle,family,eta,lam):
    tt,qt=mix_kernels(bundle,family,eta)
    eig,V=np.linalg.eigh(tt)
    assert eig.min()>-1e-8 and np.isfinite(eig).all()
    return ((qt@V)/(eig+lam))@V.T
def stats_of(base,y,corr,C):
    e=base-y;ce=C@e;cc=C@corr
    return np.array([(ce*ce).sum(),(ce*cc).sum(),(cc*cc).sum(),(e*e).sum(),(e*corr).sum(),(corr*corr).sum()],float)
def solve_beta(stats,reference_RNA=None):
    """Exact scalar quadratic minimization on beta>=0, with an optional RNA interval."""
    p0,p1,p2,r0,r1,r2=map(float,stats)
    if not np.isfinite(stats).all() or p2<0 or r2<0:raise ValueError('Invalid squared-error statistics')
    lo=0.;hi=np.inf;reason='unconstrained';tol=CONFIG['numerical_relative_tolerance']
    if reference_RNA is not None:
        ref=float(reference_RNA);D=r0-ref
        if r2==0:
            if r1==0:
                if D>tol*max(abs(r0),abs(ref),np.finfo(float).tiny):
                    return dict(eligible=False,reason='constant_RNA_exceeds_reference',beta=None,partner_sse=None,rna_sse=None,interval_low=None,interval_high=None)
            elif r1>0:hi=-D/(2*r1)
            else:lo=max(0.,-D/(2*r1))
        else:
            disc=r1*r1-r2*D
            dtol=tol*max(r1*r1,abs(r2*D),np.finfo(float).tiny)
            if disc < -dtol:
                return dict(eligible=False,reason='empty_RNA_interval',beta=None,partner_sse=None,rna_sse=None,interval_low=None,interval_high=None)
            # Treat either sign of a round-off-sized double-root discriminant
            # identically, preserving the solution under a change of units.
            if abs(disc)<=dtol:disc=0.
            sq=float(np.sqrt(max(0.,disc)))
            q=-r1-np.copysign(sq,r1)
            if q==0:roots=(-r1/r2,-r1/r2)
            else:roots=(q/r2,D/q)
            low,high=sorted(roots);lo=max(0.,low);hi=high
        if hi<lo:
            return dict(eligible=False,reason='RNA_interval_outside_nonnegative',beta=None,partner_sse=None,rna_sse=None,interval_low=None,interval_high=None)
        reason='RNA_feasible_interval'
    opt=1. if p2==0 else -p1/p2
    beta=float(np.clip(opt,lo,hi))
    # Round-off at a computed feasible boundary is corrected toward its interior.
    if reference_RNA is not None:
        for _ in range(8):
            delta=r0+2*beta*r1+beta*beta*r2-reference_RNA
            scale=max(abs(r0),abs(reference_RNA),abs(2*beta*r1),abs(beta*beta*r2),np.finfo(float).tiny)
            if delta<=tol*scale:break
            middle=(lo+hi)/2 if np.isfinite(hi) else max(lo,beta+1.)
            beta=float(np.nextafter(beta,middle))
        if delta>tol*scale:raise ArithmeticError('Quadratic interval returned infeasible beta')
    partner=p0+2*beta*p1+beta*beta*p2;rna=r0+2*beta*r1+beta*beta*r2
    return dict(eligible=True,reason=reason,beta=beta,partner_sse=partner,rna_sse=rna,
                interval_low=float(lo),interval_high=None if not np.isfinite(hi) else float(hi),
                guard_delta=None if reference_RNA is None else float(rna-reference_RNA),
                stats=dict(zip(('C_e2','C_ea','C_a2','RNA_e2','RNA_ea','RNA_a2'),map(float,stats))))
def score_at_one(stats):
    p0,p1,p2,r0,r1,r2=map(float,stats)
    return dict(eligible=True,beta=1.,partner_sse=p0+2*p1+p2,rna_sse=r0+2*r1+r2,reason='fixed_beta_one')
def pick(rows):
    good=[r for r in rows if r['eligible']];assert good,'No eligible predictor, including original-R fallback'
    best=min(r['partner_sse'] for r in good)
    tied=[r for r in good if r['partner_sse']-best<=CONFIG['numerical_relative_tolerance']*max(abs(best),abs(r['partner_sse']),np.finfo(float).tiny)]
    def priority(r):
        if r['kind']=='original_R':return (0,0,0,0,0,0,0)
        family=r['family'];complexity={'native':0,'sf12':1,'neural':1,'joint':2}[family]
        return (1,r['power']!=1,complexity,r['eta'],-r['lam'],abs(r['beta']-1),CONFIG['families'].index(family))
    return min(tied,key=priority)
def candidate_corrections(bundle,case,mask,w):
    residual=case['train_residual'][:,mask]
    return {c['name']:(w**c['power'])*(operator(bundle,c['family'],c['eta'],c['lam'])@residual) for c in candidates()}

def smoke():
    protocol=lock();started=time.monotonic();common,data,full,score,tr,q,inner=ri.outer_inputs(0)
    prep=ri.prepare_fold(0);case=prep['inner'][0];fit=case['train_global_indices'];val=case['query_global_indices'];mom=moments_load()
    rf=fit_R(common,score,fit,val,case['baseline_tune_folds'],mom)
    mask=common['score_gene_mask'];assert np.max(abs(rf['base']-case['val_base'][:,mask]))<1e-12
    bundle=build_kernels(common,data['single_effects'],score,fit,val)
    corrs=candidate_corrections(bundle,case,mask,rf['w']);assert len(corrs)==42
    for candidate in candidates():
        if candidate['power']==1:
            phase1=rf['w']*(operator(bundle,candidate['family'],candidate['eta'],candidate['lam'])@case['train_residual'][:,mask])
            assert np.array_equal(corrs[candidate['name']],phase1)
    C,blocks=C_matrix(score,val);assert len(C)
    er=rf['R']-case['val_y'][:,mask];ref=float((er*er).sum())
    for a in corrs.values():
        st=stats_of(rf['base'],case['val_y'][:,mask],a,C)
        for guard in (None,ref):
            s=solve_beta(st,guard)
            if s['eligible']:
                beta=s['beta'];pred=rf['base']+beta*a
                assert np.isclose(((C@(pred-case['val_y'][:,mask]))**2).sum(),s['partner_sse'],rtol=1e-12,atol=1e-12)
                if guard is not None:assert ((pred-case['val_y'][:,mask])**2).sum()<=guard+1e-10*max(1,guard)
    # A valid one-dimensional example requires beta=3, including RNA guard.
    st=np.array([9.,-3.,1.,9.,-3.,1.]);assert solve_beta(st,None)['beta']==3.
    assert solve_beta(st,1.)['beta']==3.
    assert not solve_beta(np.array([1.,0.,1.,4.,0.,1.]),1.)['eligible']
    for scale in (1e-9,1e-3,1e3,1e9):assert abs(solve_beta(st*scale*scale,scale*scale)['beta']-3)<1e-12
    poisoned=ri.gm.Dataset(score.pairs,score.nodes,score.y.copy(),score.additive,score.x,score.pair_kernel,score.blocks)
    poisoned.y[val]=1e9;mp={k:v.copy() for k,v in mom.items()};lookup={str(c):i for i,c in enumerate(mp['conditions'])}
    nonfit=[lookup[n] for n in score.pairs if n not in {score.pairs[i] for i in fit}]
    mp['first_means'][nonfit]=1e12;mp['second_means'][nonfit]=-1e12
    rp=fit_R(common,poisoned,fit,val,case['baseline_tune_folds'],mp)
    bp=build_kernels(common,data['single_effects'],poisoned,fit,val)
    assert np.array_equal(rf['R'],rp['R']) and np.array_equal(rf['w'],rp['w'])
    for name in ('R','SF12','NEURAL'):
        for key in ('ktt','kqt'):assert np.array_equal(bundle[name][key],bp[name][key])
    final=fit_R(common,score,tr,q,inner,mom)
    with np.load(SOURCE/'adapters/REFERENCE_M0_R.npz',allow_pickle=False) as z:
        names=z['conditions'].astype(str);lookup={n:i for i,n in enumerate(names)}
        ii=[lookup[str(n)] for n in data['query_names']]
        maxdiff=float(np.max(abs(final['R']-z['pred_effect_R'][ii])))
    assert maxdiff<1e-12
    native=operator(bundle,'native',0.,rf['metadata']['lambda'])
    native_r=rf['base']+rf['w']*(native@case['train_residual'][:,mask])
    assert np.max(abs(native_r-rf['R']))<1e-12
    report=dict(status='PASS',protocol_sha256=protocol,seconds=time.monotonic()-started,candidates=42,
                all_power_one_corrections_equal_phase1=True,
                beta_above_one_and_empty_interval_checks=True,unit_scaling_checks=4,
                validation_and_nonfit_moment_poison_difference=0.,kernel_poison_difference=0.,
                original_R_max_abs_difference=maxdiff,native_operator_original_R_difference=float(np.max(abs(native_r-rf['R']))),
                kernel_audit=bundle['audit'],outer_query_outcomes_used=False)
    write_json(ROOT/'SMOKE_AUDIT.json',report);print(json.dumps(report),flush=True)

def run_fold(fold):
    started=time.monotonic();protocol=lock();out=ROOT/'results'/f'fold_{fold}';out.mkdir(parents=True,exist_ok=True)
    if (out/'FREEZE.json').exists():
        old=json.loads((out/'FREEZE.json').read_text());assert old['protocol_sha256']==protocol
        for n,h in old['files'].items():assert sha(out/n)==h
        print(json.dumps(dict(stage='already_frozen',fold=fold)),flush=True);return
    common,data,full,score,tr,q,inner=ri.outer_inputs(fold);prep=ri.prepare_fold(fold);mask=common['score_gene_mask'];mom=moments_load()
    accum={c['name']:np.zeros(6) for c in candidates()};cache=[];inner_records=[];R_C=0.;R_RNA=0.;nrect=0
    for j,case in enumerate(prep['inner']):
        fit=case['train_global_indices'];val=case['query_global_indices']
        rf=fit_R(common,score,fit,val,case['baseline_tune_folds'],mom);base=case['val_base'][:,mask];y=case['val_y'][:,mask]
        assert np.max(abs(base-rf['base']))<1e-12
        bundle=build_kernels(common,data['single_effects'],score,fit,val)
        corrs=candidate_corrections(bundle,case,mask,rf['w']);C,blocks=C_matrix(score,val);nrect+=len(C)
        er=rf['R']-y;R_C+=float(((C@er)**2).sum());R_RNA+=float((er*er).sum())
        for name,a in corrs.items():accum[name]+=stats_of(base,y,a,C)
        cache.append(dict(base=base,y=y,R=rf['R'],C=C,corrs=corrs,conditions=case['query_names']))
        inner_records.append(dict(inner=j,train_names=case['train_names'].tolist(),validation_names=case['query_names'].tolist(),
                rectangles=blocks,R=rf['metadata'],kernel_audit=bundle['audit'],reliability_weight_sha256=array_sha(rf['w'])))
    assert nrect>0
    fallback=dict(kind='original_R',name='original_R',eligible=True,beta=1.,partner_sse=R_C,rna_sse=R_RNA)
    guard_rows=[];partner_rows=[];fixed_rows=[]
    for c in candidates():
        ident=dict(kind='candidate',**c)
        guard_rows.append(dict(**ident,**solve_beta(accum[c['name']],R_RNA)))
        partner_rows.append(dict(**ident,**solve_beta(accum[c['name']],None)))
        if c['family']=='native' and c['power']==1:fixed_rows.append(dict(**ident,**score_at_one(accum[c['name']])))
    selected=dict(R_original=fallback,R_lambda_matched=pick(fixed_rows),
                  R_calibrated=pick([fallback]+[r for r in guard_rows if r['family']=='native']),
                  R_12=pick([fallback]+[r for r in guard_rows if r['family'] in ('native','sf12')]),
                  R_neural=pick([fallback]+[r for r in guard_rows if r['family'] in ('native','neural')]),
                  R_plus_guarded=pick([fallback]+guard_rows),R_plus_partner=pick([fallback]+partner_rows),
                  R_phase1_primary=pick([fallback]+[r for r in guard_rows if r['power']==1]),
                  R_phase1_calibrated=pick([fallback]+[r for r in guard_rows if r['family']=='native' and r['power']==1]))
    final=prep['final'];rf=fit_R(common,score,tr,q,inner,mom);base=final['query_base'][:,mask]
    assert np.max(abs(base-rf['base']))<1e-12
    bundle=build_kernels(common,data['single_effects'],score,tr,q);residual=final['train_residual'][:,mask]
    state=dict(train_residual=residual,query_base=base,query_R=rf['R'],reliability_weights=rf['w'],
               train_names=final['train_names'],query_names=final['query_names'],genes=common['genes'][mask],target_names=common['target_names'])
    used={r['name']:r for r in selected.values() if r['kind']=='candidate'};ops={}
    for name,r in used.items():
        op=operator(bundle,r['family'],r['eta'],r['lam']);ops[name]=op;state['operator__'+name]=op
    pred=dict(conditions=final['query_names'],genes=common['genes'][mask],control=common['control'][mask],
              fold=np.full(len(q),fold,dtype=np.int64),pred_effect_M0=base)
    for name,r in selected.items():
        pred['pred_effect_'+name]=rf['R'].copy() if r['kind']=='original_R' else base+r['beta']*((rf['w']**r['power'])*(ops[r['name']]@residual))
    np.savez_compressed(out/'state.npz',**state);np.savez_compressed(out/'predictions.npz',**pred)
    for j,cache_j in enumerate(cache):
        arrays={k:v for k,v in cache_j.items() if k!='corrs'};unique=[];indices=[]
        for name in METHODS:
            r=selected[name];p=cache_j['R'] if r['kind']=='original_R' else cache_j['base']+r['beta']*cache_j['corrs'][r['name']]
            idx=next((i for i,v in enumerate(unique) if np.array_equal(v,p)),None)
            if idx is None:idx=len(unique);unique.append(p)
            indices.append(idx)
        arrays.update(model_names=np.array(METHODS),model_prediction_index=np.array(indices),selected_prediction_values=np.stack(unique))
        np.savez_compressed(out/f'inner_selected_{j}.npz',**arrays)
    selection=dict(fold=fold,protocol_sha256=protocol,seconds=time.monotonic()-started,selected=selected,
        guarded_candidates=guard_rows,partner_candidates=partner_rows,lambda_only_candidates=fixed_rows,
        candidate_statistics={k:v.tolist() for k,v in accum.items()},inner=inner_records,
        inner_original_R_partner_sse=R_C,inner_original_R_RNA_sse=R_RNA,inner_rectangle_count=nrect,
        final_R=rf['metadata'],final_kernel_audit=bundle['audit'],outer_query_outcomes_used=False)
    write_json(out/'selection.json',selection)
    files=['state.npz','predictions.npz','selection.json']+[f'inner_selected_{j}.npz' for j in range(3)]
    write_json(out/'FREEZE.json',dict(status='FROZEN',fold=fold,protocol_sha256=protocol,
              outer_query_outcomes_used=False,files={n:sha(out/n) for n in files}))
    print(json.dumps(dict(stage='fold_frozen',fold=fold,seconds=selection['seconds'],
          selected={k:dict(name=v['name'],beta=v['beta']) for k,v in selected.items()})),flush=True)

def aggregate():
    protocol=lock();parts=[];freezes={}
    for fold in range(5):
        p=ROOT/'results'/f'fold_{fold}';r=json.loads((p/'FREEZE.json').read_text())
        assert r['protocol_sha256']==protocol and not r['outer_query_outcomes_used']
        for n,h in r['files'].items():assert sha(p/n)==h
        freezes[str(fold)]=sha(p/'FREEZE.json')
        with np.load(p/'predictions.npz',allow_pickle=False) as z:parts.append({k:z[k] for k in z.files})
    joined={k:parts[0][k] for k in ('genes','control')}
    for p in parts:assert np.array_equal(p['genes'],joined['genes']) and np.array_equal(p['control'],joined['control'])
    keys=['conditions','fold']+['pred_effect_M0']+['pred_effect_'+m for m in METHODS]
    for k in keys:joined[k]=np.concatenate([p[k] for p in parts])
    assert len(set(joined['conditions']))==55
    path=ROOT/'upgrade_predictions_frozen.npz';np.savez_compressed(path,**joined)
    write_json(ROOT/'ALL_PREDICTIONS_FROZEN.json',dict(status='FROZEN',protocol_sha256=protocol,fold_freezes=freezes,
         predictions_file=path.name,predictions_sha256=sha(path),conditions_sha256=array_sha(joined['conditions']),
         genes_sha256=array_sha(joined['genes']),methods=['M0']+list(METHODS),rows=55,genes=4995,
         units='effect relative to common control',outer_query_outcomes_used=False))
    print('ALL_FIVE_FOLDS_FROZEN',flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--smoke',action='store_true');ap.add_argument('--fold',type=int,choices=range(5));ap.add_argument('--aggregate',action='store_true');ap.add_argument('--lock',action='store_true')
    args=ap.parse_args()
    if args.lock:print(lock())
    elif args.smoke:smoke()
    elif args.fold is not None:run_fold(args.fold)
    elif args.aggregate:aggregate()
    else:ap.error('choose --smoke, --fold, --aggregate, or --lock')
