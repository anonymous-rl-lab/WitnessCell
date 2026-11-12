"""Bounded, nested full-M0 comparison. Full means train; independent held-out
conditions are scored against each uncorrected technical-half reference.
The two references do NOT represent two independent biological experiments.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[k]='1'
import csv,hashlib,json,time,datetime,itertools
from pathlib import Path
import numpy as np
import origin_free_model_frozen as m

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BP=ROOT/'batch_preprocessor_20261020'
ARMS=('raw','B1_direct7','B2_cal9','B2_cal7')
ALPHAS=(*m.ALPHAS,None)
def dump(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for c in iter(lambda:f.read(1<<20),b''):h.update(c)
    return h.hexdigest()
def csvwrite(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def gamma_coefficients(counts,qi,bi,train,pi):
    B=counts.shape[1];O=len(qi);L=np.zeros((B,B));F=np.zeros((B,O))
    for q in train:
        p=counts[q]/counts[q].sum();L+=np.diag(p)-np.outer(p,p)
        for o in np.flatnonzero(qi==q):
            F[:,o]=-p[bi[o]]*p;F[bi[o],o]+=p[bi[o]]
    ev=np.linalg.eigvalsh(L)
    assert np.linalg.matrix_rank(L,tol=1e-10)==B-1
    A=np.linalg.pinv(L,rcond=1e-12)@F
    A-=np.ones((B,1))@(pi@A)[None,:]
    assert np.max(abs(A.sum(1)))<1e-10
    return A,{'rank':B-1,'min_positive_eigenvalue':float(ev[1]),'gauge_max':float(np.max(abs(pi@A)))}

def main():
    started=time.monotonic()
    z=np.load(BP/'inputs/compact_inputs.npz');hz=np.load(BP/'stage2_inputs/compact_inputs.npz')
    names=z['condition_ids'].astype(str);genes=z['genes'];qi=z['observation_condition'];bi=z['observation_batch']
    counts=z['counts'];pi=z['pi'];folds=z['fold_ids'];singles=np.flatnonzero(folds<0);pairs=np.flatnonzero(folds>=0)
    target_batches=np.flatnonzero(pi);counts7=counts[:,target_batches];pi7=pi[target_batches]
    obs7=np.flatnonzero(np.isin(bi,target_batches));batchmap={int(b):i for i,b in enumerate(target_batches)}
    eligible=set(json.loads((BP/'stage2_metadata_plan.json').read_text())['stage2_eligible_pairs'])
    eval_ids=np.array([q for q in pairs if names[q] in eligible]);nodes=sorted({g for q in pairs for g in names[q].split('+')})
    pairnames=names[pairs].tolist();X=m.incidence(pairnames,nodes)
    single_order=np.array([next(q for q in singles if g in names[q].split('+')) for g in nodes])
    inner_by_outer={f:m.label_folds(np.flatnonzero(folds[pairs]!=f),pairnames,k=3) for f in range(5)}
    protocol={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'purpose':'development comparison of complete tuned M0 across four measurement estimators',
      'training_measurements':'full allowed condition-by-library means; never held-out query measurements',
      'evaluation':'same held-out combinations scored independently against each raw half direct7pi reference, then condition-equal mean',
      'clarification':'Full-depth training; not half-to-opposite-half training. Two technical references, not independent biological replication.',
      'arms':list(ARMS),'pi':pi.tolist(),'source_conditions':names.tolist(),'eligible_test_conditions':names[eval_ids].tolist(),
      'outer_folds':folds.tolist(),'inner_validation_local_pair_ids':{str(f):[v.tolist() for v in fs] for f,fs in inner_by_outer.items()},
      'baseline':{'source_sha256':sha(HERE/'origin_free_model_frozen.py'),'alphas':ALPHAS,'noises':m.NOISES,'gamma_range':[0,1]},
      'inner_target':'uncorrected full fixed7pi conditional means, same all arms; validation labels used only candidate risk',
      'nested_refit':'Each innerfit recalibrates only innertrainingpairs+10singles, then rebuilds single origin and additive backbone.',
      'candidate_selection':'Analytic pooled OOF gamma; minimum risk, then larger noise, then alpha order; original source tie rules.',
      'full_M0_identity':'Same formal family and candidate grid; common-reference objective and nested preprocessing are new.',
      'run_bound_seconds':300,'pilot':'outer fold0 all arms then remaining fixed folds if runtime/invariants permit, regardless of gains',
      'never_access':'original15query expression, real test labels, other data files',
      'code_sha256':sha(__file__)}
    dump(HERE/'PROTOCOL_BEFORE_FIT.json',protocol)
    # Expression arrays are accessed only after above protocol is written.
    Y=np.load(BP/'inputs/allowed_batch_means.npy',mmap_mode='r')
    HY=np.load(BP/'stage2_inputs/allowed_half_means.npy',mmap_mode='r');G=Y.shape[1]
    p=counts/counts.sum(1)[:,None];p7=counts7/counts7.sum(1)[:,None]
    pool=np.zeros((len(names),len(qi)));std=np.zeros_like(pool);pool7=np.zeros_like(pool)
    for o,(q,b) in enumerate(zip(qi,bi)):
        pool[q,o]=p[q,b];std[q,o]=pi[b]
        if b in batchmap:pool7[q,o]=p7[q,batchmap[int(b)]]
    assert np.max(abs(std.sum(1)-1))<1e-12
    fulltarget=std@Y
    reference={}
    for q in eval_ids:
      for h in (0,1):
        C=np.zeros(len(hz['half']));rows=np.flatnonzero((hz['half']==h)&(hz['condition']==q))
        C[rows]=pi[hz['batch'][rows]];assert abs(C.sum()-1)<1e-12
        reference[(int(q),h)]=C@HY
    prep_cache={};prep_audit=[];prep_states={}
    def preprocess(arm,train):
        key=(arm,tuple(map(int,train)))
        if key in prep_cache:return prep_cache[key]
        permitted=np.r_[singles,pairs[train]]
        excluded=np.setdiff1d(pairs,pairs[train])
        if arm=='raw':L=pool.copy();design={}
        elif arm=='B1_direct7':L=std.copy();design={}
        elif arm=='B2_cal9':
            A,design=gamma_coefficients(counts,qi,bi,permitted,pi)
            assert np.max(abs(A[:,np.isin(qi,excluded)]))==0
            L=pool-(p-pi)@A
        elif arm=='B2_cal7':
            A7,design=gamma_coefficients(counts7,qi[obs7],np.array([batchmap[int(b)] for b in bi[obs7]]),permitted,pi7)
            A=np.zeros((7,len(qi)));A[:,obs7]=A7
            assert np.max(abs(A[:,np.isin(qi,excluded)]))==0
            L=pool7-(p7-pi7)@A
        use=np.r_[single_order,pairs[train]]
        # Do not materialize calibrated query means. Only training/single rows.
        assert np.max(abs(L[use][:,np.isin(qi,excluded)]))==0
        v=L[use]@Y;single=v[:10];origin=single.mean(0);e=single-origin
        fit_y=v[10:]-origin;additive=.5*X@e
        result=(origin,fit_y,additive)
        prep_cache[key]=result
        statekey=f'{arm}_'+hashlib.sha256(np.asarray(train,dtype=np.int64).tobytes()).hexdigest()[:10]
        prep_states[statekey]=L[use]
        prep_audit.append({'key':statekey,'arm':arm,'train_global':pairs[train].tolist(),'train_conditions':names[pairs[train]].tolist(),
                           'excluded_pairs':names[excluded].tolist(),'source_nonzero_excluded':False,**design})
        return result
    def back(additive,q,a):return np.zeros((len(q),G)) if a is None else m.sat(additive[q],a)
    selection=[];stats={};candidates=[];score_rows=[];pred_rows=[];pred_meta=[];examples={};max_direct=0.;max_fit_direct=0.
    for f in range(5):
      train=np.flatnonzero(folds[pairs]!=f);query=np.array([i for i,q in enumerate(pairs) if folds[q]==f and names[q] in eligible])
      for arm in ARMS:
        if time.monotonic()-started>300:raise TimeoutError('frozen 300s compute bound')
        dot=np.zeros((7,5));raw2=np.zeros_like(dot);truth2=np.zeros(7)
        for innernum,va in enumerate(inner_by_outer[f]):
          fit=np.array(sorted(set(train)-set(va)),int);origin,yy,add=preprocess(arm,fit)
          xt=X[fit];xm=xt.mean(0);xc=xt-xm;xv=X[va]-xm;xx=xv.T@xv
          invs=[np.linalg.inv(xc.T@xc+2*eta*np.eye(10)) for eta in m.NOISES]
          for ai,alpha in enumerate(ALPHAS):
            r=yy-back(add,fit,alpha);u=r.mean(0);r-=u
            target=fulltarget[pairs[va]]-origin-back(add,va,alpha)-u
            q=xc.T@r;t=xv.T@target;cross=t@q.T;qq=q@q.T
            truth2[ai]+=float(np.sum(target*target))
            for ni,inv in enumerate(invs):
              ds=float(np.sum(inv*cross));rs=float(np.sum((inv@qq@inv.T)*xx));dot[ai,ni]+=ds;raw2[ai,ni]+=rs
              if ai in (0,6) and ni in (0,4):
                direct=xv@inv@q
                direct_stats=np.array([np.sum(target*direct),np.sum(direct*direct)])
                scale=max(1.,abs(ds),abs(rs));max_direct=max(max_direct,float(np.max(abs(direct_stats-[ds,rs]))/scale))
        gamma=np.divide(dot,raw2,out=np.zeros_like(dot),where=raw2>1e-12).clip(0,1)
        loss=truth2[:,None]-2*gamma*dot+gamma**2*raw2
        best=min(np.ndindex(loss.shape),key=lambda ij:(loss[ij],-m.NOISES[ij[1]],ij[0]))
        ai,ni=best;pars={'alpha':ALPHAS[ai],'noise':m.NOISES[ni],'gamma':float(gamma[ai,ni])}
        key=f'fold{f}__{arm}';stats[key+'__truth2']=truth2;stats[key+'__dot']=dot;stats[key+'__raw2']=raw2
        for aidx,nidx in np.ndindex(loss.shape):
          candidates.append({'fold':f,'arm':arm,'alpha':ALPHAS[aidx],'noise':m.NOISES[nidx],'gamma':float(gamma[aidx,nidx]),
                             'oof_sse':float(loss[aidx,nidx]),'oof_mse':float(loss[aidx,nidx]/(len(train)*G)),
                             'selected':bool((aidx,nidx)==best)})
        selection.append({'fold':f,'arm':arm,'selected':pars,'index':[int(ai),int(ni)],'oof_mse':float(loss[best]/(len(train)*G)),
                          'outertrain':names[pairs[train]].tolist(),'outertest':names[pairs[query]].tolist()})
        origin,yy,add=preprocess(arm,train);xt=X[train];xm=xt.mean(0);xc=xt-xm
        r=yy-back(add,train,pars['alpha']);u=r.mean(0)
        beta=np.linalg.solve(xc.T@xc+2*pars['noise']*np.eye(10),xc.T@(r-u))
        pred=origin+back(add,query,pars['alpha'])+u+pars['gamma']*(X[query]-xm)@beta
        # Independent comparison to the unmodified original formal baseline.
        singles_abs=add[:0] # no-op, reconstruction below uses saved single rows.
        _,_,_=origin,yy,add
        # Data-like namespace avoids ever filling query truth with measurements.
        from types import SimpleNamespace
        all_y=np.zeros((25,G));all_y[train]=yy
        d=SimpleNamespace(x=X,y=all_y,additive=add)
        native=origin+m.baseline(d,train,query,pars)
        max_fit_direct=max(max_fit_direct,float(np.max(abs(pred-native))))
        for j,loc in enumerate(query):
          q=int(pairs[loc]);pred_rows.append(pred[j].copy());pred_meta.append({'type':'prediction','fold':f,'arm':arm,'condition':names[q],'global_id':q})
          for h in (0,1):
            target=reference[q,h];v=float(np.mean((pred[j]-target)**2))
            score_rows.append({'fold':f,'arm':arm,'condition':names[q],'half':h,'mse':v})
          examples[key+'__'+str(q)]=np.stack([pred[j,z['example_idx']],reference[q,0][z['example_idx']],reference[q,1][z['example_idx']]])
        print(json.dumps({'fold':f,'arm':arm,'selected':pars,'n_scored':len(query),'seconds':round(time.monotonic()-started,3)}),flush=True)
      if f==0:
        dump(HERE/'PILOT.json',{'seconds':time.monotonic()-started,'direct_quad_relative_error':max_direct,'native_baseline_max_abs':max_fit_direct,
                              'decision':'continue fixed fivefold schedule; no result-dependent changes'})
    assert max_direct<1e-10 and max_fit_direct<1e-10
    # Compact gene-space Gram is sufficient for exact all-gene score replay.
    for (q,h),v in reference.items():pred_rows.append(v);pred_meta.append({'type':'reference','condition':names[q],'global_id':q,'half':h})
    matrix=np.asarray(pred_rows);matrix-=matrix.mean(0)
    gram=matrix@matrix.T/G
    np.savez_compressed(HERE/'SCORE_GRAM.npz',gram=gram)
    dump(HERE/'SCORE_GRAM_ROWS.json',pred_meta)
    np.savez_compressed(HERE/'SELECTION_QUADRATICS.npz',**stats)
    np.savez_compressed(HERE/'PREPROCESS_COEFFICIENTS.npz',**prep_states)
    np.savez_compressed(HERE/'FIVE_GENE_EXAMPLES.npz',**examples)
    dump(HERE/'SELECTION.json',selection);dump(HERE/'PREPROCESS_FITS.json',prep_audit)
    csvwrite(HERE/'ALL_CANDIDATES.csv',candidates);csvwrite(HERE/'PER_CONDITION_HALF_ERRORS.csv',score_rows)
    summary={arm:float(np.mean([r['mse'] for r in score_rows if r['arm']==arm])) for arm in ARMS}
    halfsummary={arm:{str(h):float(np.mean([r['mse'] for r in score_rows if r['arm']==arm and r['half']==h])) for h in(0,1)} for arm in ARMS}
    result={'status':'complete','all_gene_count':G,'n_eligible_conditions':len(eval_ids),'arms_mse':summary,'half_mse':halfsummary,
            'gains_vs_raw_pct':{a:100*(1-summary[a]/summary['raw']) for a in ARMS},
            'cal9_gain_vs_direct7_pct':100*(1-summary['B2_cal9']/summary['B1_direct7']),
            'cal7_gain_vs_direct7_pct':100*(1-summary['B2_cal7']/summary['B1_direct7']),
            'direct_quadratic_relative_error':max_direct,'original_source_prediction_max_abs':max_fit_direct,
            'compute_seconds':time.monotonic()-started,'old15query_read':False,'fit_count':len(prep_audit)}
    dump(HERE/'SUMMARY.json',result)
    sources=[BP/'inputs/compact_inputs.npz',BP/'inputs/allowed_batch_means.npy',BP/'stage2_inputs/compact_inputs.npz',BP/'stage2_inputs/allowed_half_means.npy']
    dump(HERE/'SOURCE_INPUT_HASHES.json',{str(p.relative_to(ROOT)):sha(p) for p in sources})
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
