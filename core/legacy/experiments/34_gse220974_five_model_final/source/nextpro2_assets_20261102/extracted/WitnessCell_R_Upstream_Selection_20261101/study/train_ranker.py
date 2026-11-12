#!/usr/bin/env python3
"""Train a small graph-feature residual ranker without outer query answers.

The objective is fixed by reproducible single-intervention response directions.
The first-stage R-G predictions are the unchanged frozen outer predictions.
The second stage learns 14 utility residuals on genuinely out-of-fit inner R-G
predictions. Its leave-inner-block calibration is internal calibration, not a
fully nested performance estimate: nuisance fits for other blocks overlap.
"""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parent
DEFAULT_SOURCE=ROOT.parent/'r_edge_stability_fine_mix_20261101'
RIDGES=(.1,1.,10.)
FRACTIONS=(0.,.25,.5,1.)
GRAPH_RIDGE=1.
EDGE_THRESHOLD=.5
CONFIG=dict(objective='L1-normalized, half-agreement-weighted A single effect on the common 4994 non-target RNA columns',
    graph_ridge=GRAPH_RIDGE,edge_threshold=EDGE_THRESHOLD,graph_hops='coupled multi-hop port solver',
    feature_projection='remove train-fit endpoint-linear component; mean training squared row norm one',
    ranker='14 graph features -> 14 fixed-program utility residuals; no intercept; ridge',
    ranker_ridges=list(RIDGES),correction_fractions=list(FRACTIONS),
    calibration='Leave one inner combination block out on OOF features and OOF R-G residuals; internal calibration, not a fully nested estimate',
    calibration_target='mean of within-A centered utility-residual MSEs across inner-block/A groups with at least two partners; only combinations containing A enter that program group',
    tie='smaller correction fraction, then stronger ridge',
    final_fit='fit selected ridge on OOF graph features and OOF residuals, apply to final outer-training graph features',
    outer_query_outcomes_used=False)

def configure(source_root=DEFAULT_SOURCE):
    global SOURCE,b,ri,ef,es,rm,br
    SOURCE=Path(source_root).resolve()
    sys.path.insert(0,str(SOURCE))
    import baseline_core as _b
    import edge_features as _ef
    import edge_stability as _es
    import relation_model as _rm
    import bridge_features as _br
    b=_b;ri=b.ri;ef=_ef;es=_es;rm=_rm;br=_br
    return SOURCE

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')

def objective_vectors(common,moments):
    """No double-condition row enters an objective direction."""
    names=common['target_names'].astype(str);genes=common['genes'].astype(str)
    lookup={str(x):i for i,x in enumerate(moments['conditions'])}
    ids=[lookup[n] for n in names];ci=lookup['control']
    full=moments['means'][ids]-moments['means'][ci]
    h1=moments['first_means'][ids]-moments['first_means'][ci]
    h2=moments['second_means'][ids]-moments['second_means'][ci]
    den=h1*h1+h2*h2
    agreement=np.divide(2*h1*h2,den,out=np.zeros_like(den),where=den>0).clip(0,1)
    agreement[h1*h2<=0]=0
    mask=common['score_gene_mask'].astype(bool)&~np.isin(genes,names)
    assert mask.sum()==4994
    weighted=full*agreement;weighted[:,~mask]=0.
    norms=np.abs(weighted).sum(axis=1)
    assert np.all(norms>0)
    V=weighted/norms[:,None]
    assert np.max(abs(np.abs(V).sum(axis=1)-1))<1e-12
    return V,dict(agreement=agreement,full_single_effect=full,common_mask=mask,l1_norms=norms)

def training_moments(fit_names):
    """Erase all non-training double values before passing moments downstream."""
    m=br.load_moments();names=m['conditions'].astype(str)
    forbidden=np.array(['+' in n and n not in set(fit_names) for n in names])
    for key,value in m.items():
        if key not in ('conditions','genes') and value.ndim and len(value)==len(names):
            if value.dtype.kind in 'fc':value[forbidden]=np.nan
            elif value.dtype.kind in 'iu':value[forbidden]=0
    return m

def native_blocks(score,fit,query):
    k=score.pair_kernel;x=score.x
    tt,qt,_,_=ri.gm.projected_blocks(x[fit],x[query],k[np.ix_(fit,fit)],
        k[np.ix_(query,fit)],k[np.ix_(query,query)])
    return tt,qt

def native_corrections(common,full,score,fit,query,tune_folds,moments):
    """Six native R-G candidates, fitted only using fit responses."""
    rf=b.fit_R(common,score,fit,query,tune_folds,moments)
    d=b.clone_fit(score,fit)
    pars=rf['metadata']['baseline'] if 'baseline' in rf['metadata'] else ri.gm.tune_baseline(d,fit,tune_folds)
    base_fit=ri.gm.baseline_predict(d,fit,fit,pars)
    residual=score.y[fit]-base_fit
    tt,qt=native_blocks(score,fit,query);eig,Q=np.linalg.eigh(tt)
    candidates={}
    for lam in (.1,1.,10.):
        raw=(((qt@Q)/(eig+lam))@Q.T)@residual
        for power in (1,2):
            candidates[(lam,power)]=(rf['w']**power)*raw
    return rf,candidates

def fit_native_rg(common,data,full,score,train,query,tune_folds,moments):
    """Rebuild native R-G selection strictly inside the supplied train set."""
    train=np.asarray(train,int);query=np.asarray(query,int)
    assert not set(train)&set(query)
    assert all(set(v)<=set(train) for v in tune_folds)
    stats={(lam,power):np.zeros(6) for lam in (.1,1.,10.) for power in (1,2)}
    old_C=0.;old_RNA=0.;nrect=0;inner_meta=[]
    for val in tune_folds:
        val=np.asarray(val,int);fit=np.array(sorted(set(train)-set(val)),int)
        if not len(val):continue
        aux=ri.gm.preserved_folds(fit,score.pairs,score.blocks)
        rf,corrections=native_corrections(common,full,score,fit,val,aux,moments)
        C,blocks=b.C_matrix(score,val);nrect+=len(C)
        y=score.y[val];assert np.isfinite(y).all()
        er=rf['R']-y;old_C+=float(((C@er)**2).sum());old_RNA+=float((er*er).sum())
        for key,corr in corrections.items():stats[key]+=b.stats_of(rf['base'],y,corr,C)
        inner_meta.append(dict(fit=fit.tolist(),validation=val.tolist(),rectangles=len(C)))
    fallback=dict(kind='original_R',name='original_R',eligible=True,beta=1.,partner_sse=old_C,rna_sse=old_RNA)
    rows=[fallback]
    for (lam,power),st in stats.items():
        rows.append(dict(kind='candidate',family='native',eta=0.,lam=lam,power=power,
            name=f'native_l{lam:g}_p{power}',**b.solve_beta(st,old_RNA)))
    selected=b.pick(rows)
    rf,corrections=native_corrections(common,full,score,train,query,tune_folds,moments)
    pred=rf['R'] if selected['kind']=='original_R' else rf['base']+selected['beta']*corrections[(selected['lam'],selected['power'])]
    assert np.isfinite(pred).all()
    return pred,dict(selected=selected,candidates=rows,inner=inner_meta,inner_rectangles=nrect,
                     fitting_ids=train.tolist(),query_ids=query.tolist(),outer_query_outcomes_used=False)

def graph_features(common,data,full,fit,moments):
    """Canonical 14-port coordinates, with fitted endpoint component removed."""
    fit=np.sort(np.asarray(fit,int))
    inp,evidence=ef.fitting_evidence(common,data,full,fit,moments)
    stable=es.fit_edge_stability(**evidence,ridge=GRAPH_RIDGE,seed=ef.BOOTSTRAP_SEED,
        n_bootstrap=ef.BOOTSTRAP_COUNT,mode='masked')
    masked=es.threshold_relation(stable['W'],stable['edge_stability'],EDGE_THRESHOLD)['W']
    E=inp['single_ports']/inp['scales']
    pairs=np.array([np.flatnonzero(u) for u in full.x],int)
    zero=rm.predict_pair_features(np.zeros_like(masked),E,pairs,solve_ridge=br.SOLVE_RIDGE)['prediction']
    pred=rm.predict_pair_features(masked,E,pairs,solve_ridge=br.SOLVE_RIDGE)['prediction']
    raw=pred-zero
    endpoint_coef=np.linalg.pinv(full.x[fit])@raw[fit]
    F=raw-full.x@endpoint_coef
    norm=float(np.sqrt(np.mean(np.sum(F[fit]*F[fit],axis=1))))
    raw_norm=float(np.sqrt(np.mean(np.sum(raw[fit]*raw[fit],axis=1))))
    if norm*norm<=1e-12*max(raw_norm*raw_norm,1e-300):F=np.zeros_like(F);norm=0.
    else:F=F/norm
    assert np.isfinite(F).all()
    assert np.max(abs(full.x[fit].T@F[fit]))<1e-8
    state={k:v for k,v in stable.items() if isinstance(v,np.ndarray)}
    state.update(masked_W=masked,scales=inp['scales'],single_ports=inp['single_ports'],
        raw_features=raw,projected_features=F,endpoint_coefficients=endpoint_coef,
        feature_normalization=np.array(norm),fit_ids=fit)
    return F,state

def ridge_fit(X,Y,lam):
    assert X.ndim==2 and Y.ndim==2 and len(X)==len(Y)
    return np.linalg.solve(X.T@X+float(lam)*np.eye(X.shape[1]),X.T@Y)

def ranking_calibration_loss(error,groups,incidence):
    """Equal weight for each eligible inner-block/A choice task."""
    values=[];contrasts=0
    for group in np.unique(groups):
        for a in range(error.shape[1]):
            chosen=(groups==group)&(incidence[:,a]>0)
            if chosen.sum()<2:continue
            e=error[chosen,a];values.append(float(np.mean((e-e.mean())**2)))
            contrasts+=int(chosen.sum()-1)
    assert values,'No training-side fixed-A groups for ranker calibration'
    return float(np.mean(values)),len(values),contrasts

def calibrate(X,Y,groups,incidence):
    rows=[];preds={}
    for lam in RIDGES:
        cv=np.zeros_like(Y)
        for group in np.unique(groups):
            va=groups==group;tr=~va
            cv[va]=X[va]@ridge_fit(X[tr],Y[tr],lam)
        preds[lam]=cv
        for fraction in FRACTIONS:
            error=fraction*cv-Y
            loss,ngroups,ncontrasts=ranking_calibration_loss(error,groups,incidence)
            rows.append(dict(ridge=lam,fraction=fraction,sse=loss,
                all_program_sse=float((error*error).sum()),choice_groups=ngroups,centered_contrasts=ncontrasts))
    best=min(x['sse'] for x in rows)
    tol=64*np.finfo(float).eps*max(best,np.finfo(float).tiny)
    selected=min((r for r in rows if r['sse']<=best+tol),key=lambda r:(r['fraction'],-r['ridge']))
    return selected,rows,preds[selected['ridge']]

def run_fold(fold,out_root,reference_root=None):
    started=time.monotonic();out=Path(out_root)/f'fold_{fold}';out.mkdir(parents=True,exist_ok=True)
    common,data,full,score,tr,q,inner=ri.outer_inputs(fold)
    assert np.isnan(full.y[q]).all() and np.isnan(score.y[q]).all()
    moments=training_moments([full.pairs[i] for i in tr])
    V,vaudit=objective_vectors(common,moments);mask=common['score_gene_mask'].astype(bool);Vm=V[:,mask]
    assert np.all(V[:,~mask]==0)
    X=[];Y=[];RG=[];truth=[];ids=[];groups=[];inner_meta=[];inner_states={}
    for j,val in enumerate(inner):
        fit=np.array(sorted(set(tr)-set(val)),int)
        aux=ri.gm.preserved_folds(fit,score.pairs,score.blocks)
        pred,meta=fit_native_rg(common,data,full,score,fit,val,aux,moments)
        F,state=graph_features(common,data,full,fit,moments)
        actual=score.y[val]@Vm.T;native=pred@Vm.T
        X.append(F[val]);Y.append(actual-native);RG.append(native);truth.append(actual)
        ids.extend(val.tolist());groups.extend([j]*len(val))
        inner_meta.append(dict(inner=j,train_ids=fit.tolist(),validation_ids=val.tolist(),native_rg=meta,
            retained_edges=int(np.count_nonzero(state['masked_W']))))
        inner_states.update({f'inner{j}__{k}':v for k,v in state.items()})
    X=np.concatenate(X);Y=np.concatenate(Y);RG=np.concatenate(RG);truth=np.concatenate(truth)
    ids=np.array(ids,int);groups=np.array(groups,int)
    assert set(ids)==set(tr) and len(ids)==len(set(ids))
    chosen,grid,cv=calibrate(X,Y,groups,full.x[ids])
    coef=ridge_fit(X,Y,chosen['ridge'])
    F,state=graph_features(common,data,full,tr,moments)
    reference_root=Path(reference_root) if reference_root else SOURCE
    reference_path=reference_root/'results'/f'fold_{fold}'/'predictions.npz'
    with np.load(reference_path,allow_pickle=False) as z:
        assert z['conditions'].astype(str).tolist()==[full.pairs[i] for i in q]
        assert np.array_equal(z['genes'],common['genes'][mask])
        outer_rg=z['pred_effect_R_G'].copy()
    utility_rg=outer_rg@Vm.T
    delta=chosen['fraction']*(F[q]@coef)
    utility_updated=utility_rg+delta
    np.savez_compressed(out/'predictions.npz',query_names=np.array([full.pairs[i] for i in q]),
        target_names=common['target_names'],query_ids=q,fold=np.array(fold),
        utility_R_G=utility_rg,utility_updated=utility_updated,utility_correction=delta,
        utility_single=vaudit['full_single_effect']@V.T,
        W=state['W'],stability=state['edge_stability'],masked_W=state['masked_W'])
    np.savez_compressed(out/'state.npz',objective_vectors=V,objective_score_vectors=Vm,
        objective_full_single_effect=vaudit['full_single_effect'],objective_agreement=vaudit['agreement'],
        objective_l1_norms=vaudit['l1_norms'],common_mask=vaudit['common_mask'],genes=common['genes'],
        target_names=common['target_names'],all_pair_names=common['all_pair_names'],train_ids=tr,query_ids=q,
        oof_ids=ids,oof_groups=groups,oof_incidence=full.x[ids],oof_features=X,oof_residual=Y,oof_R_G=RG,oof_truth=truth,
        oof_calibration_prediction=cv,coefficient=coef,chosen_ridge=np.array(chosen['ridge']),
        chosen_fraction=np.array(chosen['fraction']),final_features=F,
        frozen_outer_R_G_prediction=outer_rg,**{f'final__{k}':v for k,v in state.items()},**inner_states)
    meta=dict(fold=fold,config=CONFIG,selected=chosen,calibration_candidates=grid,inner=inner_meta,
        seconds=time.monotonic()-started,outer_query_outcomes_used=False,
        first_stage_reference=str(reference_path),first_stage_reference_sha256=sha(reference_path),
        train_names=[full.pairs[i] for i in tr],query_names=[full.pairs[i] for i in q],
        endpoint_projection_max_abs=float(np.max(abs(full.x[tr].T@F[tr]))),
        calibration_is_not_fully_nested=True,development_panel=True)
    write_json(out/'selection.json',meta)
    write_json(out/'FREEZE.json',dict(status='FROZEN_BEFORE_OUTER_UTILITY_SCORING',fold=fold,
        code_sha256=sha(Path(__file__)),outer_query_outcomes_used=False,
        files={n:sha(out/n) for n in ('predictions.npz','state.npz','selection.json')}))
    print(json.dumps(dict(fold=fold,selected=chosen,seconds=meta['seconds'],frozen=True)),flush=True)
    return meta

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-root',type=Path,default=DEFAULT_SOURCE)
    ap.add_argument('--output',type=Path,default=ROOT/'ranker')
    ap.add_argument('--reference-root',type=Path)
    ap.add_argument('--fold',type=int,action='append')
    args=ap.parse_args();configure(args.source_root)
    started=time.monotonic();records=[]
    for fold in (args.fold if args.fold is not None else range(5)):
        if fold not in range(5):raise ValueError('fold must be 0..4')
        records.append(run_fold(fold,args.output,args.reference_root))
    write_json(args.output/'TRAINING_RUN.json',dict(config=CONFIG,seconds=time.monotonic()-started,
        folds=[r['fold'] for r in records],outer_query_outcomes_used=False,code_sha256=sha(Path(__file__))))

if __name__=='__main__':main()
