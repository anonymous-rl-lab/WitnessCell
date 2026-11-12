#!/usr/bin/env python3
"""Standard endpoint ridge, two predeclared inner-selection objectives.

This is a transparent task-specific linear baseline, not an implementation of
the NM2025 GWP model. Outputs are effects relative to the common control.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_name]='1'
from pathlib import Path
import argparse,hashlib,importlib.util,json,time
import numpy as np

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'r_edge_stability_fine_mix_20261101'/'source'/'adapters'
OBJECTIVES=ROOT.parent/'pairpick_20261101'/'ranker'
RIDGES=(0.,.001,.01,.1,1.,10.,100.)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+'\n')
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def fit_ridge(u,y,ridge):
    """Minimize sum_i ||y_i-b-u_i W||² + ridge*||W||².

    One row per condition, equal weights. Centering leaves the intercept free.
    """
    u=np.asarray(u,float);y=np.asarray(y,float)
    if u.ndim!=2 or y.ndim!=2 or len(u)!=len(y) or not len(u):raise ValueError('nonempty 2D paired inputs required')
    if not np.isfinite(u).all() or not np.isfinite(y).all():raise ValueError('nonfinite input')
    if isinstance(ridge,(bool,np.bool_)) or not np.isfinite(ridge) or ridge<0:raise ValueError('nonnegative real ridge required')
    ux=u.mean(axis=0);ym=y.mean(axis=0);x=u-ux;v=y-ym
    if ridge==0:w=np.linalg.lstsq(x,v,rcond=None)[0]
    else:w=np.linalg.solve(x.T@x+ridge*np.eye(u.shape[1]),x.T@v)
    return np.vstack([ym-ux@w,w])

def predict(u,coef):
    u=np.asarray(u,float);coef=np.asarray(coef,float)
    if u.ndim!=2 or coef.ndim!=2 or u.shape[1]+1!=coef.shape[0] or not np.isfinite(u).all() or not np.isfinite(coef).all():raise ValueError('prediction shapes or values')
    return coef[0]+u@coef[1:]

def fit_subset(full,singles,fit,ridge):
    fit=np.asarray(fit,int);n=full.x.shape[1]
    if len(fit)!=len(set(fit.tolist())):raise ValueError('duplicate training condition')
    u=np.vstack([np.zeros((1,n)),np.eye(n),full.x[fit]])
    y=np.vstack([np.zeros((1,singles.shape[1])),singles,full.y[fit]])
    return fit_ridge(u,y,ridge)

def choose_max(scores):
    best=max(scores.values());return min(k for k,v in scores.items() if v>=best-1e-12)

def decision_rows(names,pair_names,predicted,truth,V,block):
    p=predicted@V.T;t=truth@V.T;rows=[]
    for ai,A in enumerate(names):
        ids=[i for i,pair in enumerate(pair_names) if A in str(pair).split('+')]
        if len(ids)<2:continue
        partner={next(x for x in str(pair_names[i]).split('+') if x!=A):i for i in ids}
        chosen=choose_max({b:float(p[i,ai]) for b,i in partner.items()})
        rows.append(dict(block=block,A=A,candidates=sorted(partner),selected=chosen,
                         predicted_utility=float(p[partner[chosen],ai]),actual_utility=float(t[partner[chosen],ai])))
    return rows

def select_ridge(records,key,maximize=False):
    optimum=(max if maximize else min)(r[key] for r in records)
    tol=64*np.finfo(float).eps*max(1.,max(abs(r[key]) for r in records))
    return max((r for r in records if abs(r[key]-optimum)<=tol),key=lambda r:r['ridge'])

def smoke(interface):
    checks=[];rng=np.random.default_rng(42)
    u=np.vstack([np.zeros((1,3)),np.eye(3),[[1,1,0],[0,1,1]]]);w=rng.normal(size=(3,5));b=rng.normal(size=5)
    coef=fit_ridge(u,b+u@w,0.);np.testing.assert_allclose(predict(np.array([[1,0,1]]),coef),b+np.array([[1,0,1]])@w,atol=1e-12)
    checks.append('heldout_additive_synthetic_recovered')
    constant=np.tile(b,(len(u),1));coef=fit_ridge(u,constant,100.)
    np.testing.assert_allclose(coef[0],b,atol=1e-14);np.testing.assert_allclose(coef[1:],0,atol=1e-14)
    checks.append('intercept_unpenalized_even_under_strong_ridge')
    common,data,full,score,tr,q,inner=interface.outer_inputs(0);fit=np.setdiff1d(tr,inner[0])
    original=fit_subset(full,data['single_effects'],fit,1.)
    backup=full.y.copy();forbidden=np.setdiff1d(np.arange(len(full.y)),fit)
    full.y[forbidden]=rng.normal(size=full.y[forbidden].shape)*1e8
    changed=fit_subset(full,data['single_effects'],fit,1.)
    np.testing.assert_array_equal(original,changed);full.y[:]=backup
    assert np.isnan(full.y[q]).all() and int(common['score_gene_mask'].sum())==4995
    checks.append('all_forbidden_double_labels_poison_invariant_and_outer_masked')
    errors=0
    for x,y,lam in [(u,constant,True),(u,constant,-1),(u,np.full_like(constant,np.nan),1),(u,constant[:1],1)]:
        try:fit_ridge(x,y,lam)
        except ValueError:errors+=1
    assert errors==4
    checks.append('invalid_type_value_and_shape_rejected')
    return dict(status='PASS',operations=4,checks=checks,outer_query_labels_read=False,
                poison_scope='all double outcomes outside inner fit including inner validation and outer query; fit invariance only')

def run_fold(fold,interface,source,objectives,output):
    folder=output/f'fold_{fold}';folder.mkdir(parents=True,exist_ok=True)
    if (folder/'FREEZE.json').exists():raise FileExistsError('already frozen: '+str(folder))
    common,data,full,score,tr,q,inner=interface.outer_inputs(fold)
    mask=common['score_gene_mask'].astype(bool);singles=data['single_effects'];names=common['target_names'].astype(str).tolist()
    objective_file=objectives/f'fold_{fold}'/'state.npz'
    with np.load(objective_file,allow_pickle=False) as z:
        V=z['objective_vectors'].copy();genes=z['genes'].copy();target_names=z['target_names'].copy();common_mask=z['common_mask'].copy()
    np.testing.assert_array_equal(genes,common['genes']);np.testing.assert_array_equal(target_names,common['target_names'])
    assert V.shape==(14,5020) and np.isfinite(V).all() and mask.sum()==4995 and common_mask.sum()==4994
    assert np.all(V[:,~common_mask]==0) and np.isnan(full.y[q]).all()
    records=[];all_inner=[]
    for lam in RIDGES:
        sse=0.;elements=0;decisions=[];folds=[]
        for block,val in enumerate(inner):
            fit=np.setdiff1d(tr,val);coef=fit_subset(full,singles,fit,lam);prediction=predict(full.x[val],coef)
            err=(prediction[:,mask]-full.y[val][:,mask])**2;sse+=float(err.sum());elements+=err.size
            rows=decision_rows(names,[full.pairs[i] for i in val],prediction,full.y[val],V,block);decisions.extend(rows)
            folds.append(dict(block=block,fit_ids=fit.tolist(),validation_ids=val.tolist(),rna_mse=float(err.mean()),choice_groups=len(rows)))
        assert decisions and elements>0
        record=dict(ridge=lam,rna_mse=sse/elements,choice_utility=float(np.mean([r['actual_utility'] for r in decisions])),
                    inner_groups=len(decisions),validation_elements=elements)
        records.append(record);all_inner.append(dict(ridge=lam,folds=folds,decisions=decisions))
    selected={'linear_rna':select_ridge(records,'rna_mse'),'linear_choice':select_ridge(records,'choice_utility',True)}
    pred={};coefs={}
    for model,selection in selected.items():
        coef=fit_subset(full,singles,tr,selection['ridge']);pred['pred_effect_'+model]=predict(full.x[q],coef);coefs['coef_'+model]=coef
    np.savez_compressed(folder/'predictions.npz',conditions=np.asarray([full.pairs[i] for i in q]),genes=common['genes'],fold=np.repeat(fold,len(q)),
                        control=common['control'],query_ids=q,target_names=common['target_names'],**pred,**coefs)
    meta=dict(fold=fold,algorithm='standard endpoint ridge; task-specific implementation, not NM2025 GWP',
        objective='sum condition squared errors + lambda*sum endpoint coefficient squared; intercept unpenalized',
        design='control zero + 14 singles + permitted training doubles; one equal-weight row per condition',
        output='effect relative to common control; full 5020 genes',ridges=list(RIDGES),selected=selected,candidates=records,inner=all_inner,
        train_ids=tr.tolist(),query_ids=q.tolist(),train_names=[full.pairs[i] for i in tr],query_names=[full.pairs[i] for i in q],
        training_condition_count=len(tr)+15,selection_masks=dict(rna=4995,choice=4994),
        ties='within 64*float64-eps*max(1,abs(values)); stronger ridge; predicted partner numeric ties 1e-12 alphabetical',
        objective_source_used_keys=['objective_vectors','genes','target_names','common_mask'],objective_source_sha256=sha(objective_file),
        input_sha256={name:sha(source/'sf12_original'/'inputs'/name) for name in ['common.npz',f'fold_{fold}.npz','frozen_original_design.json']},
        outer_query_outcomes_used=False)
    write_json(folder/'selection.json',meta)
    write_json(folder/'FREEZE.json',dict(status='FROZEN_BEFORE_EXTERNAL_SCORING',code_sha256=sha(__file__),outer_query_outcomes_used=False,
        files={p.name:sha(p) for p in [folder/'predictions.npz',folder/'selection.json']}))
    print(json.dumps(dict(fold=fold,selected=selected,query_predictions=len(q),outer_query_outcomes_used=False)),flush=True)
    return common,q,pred,meta

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,default=SOURCE)
    parser.add_argument('--objectives',type=Path,default=OBJECTIVES);parser.add_argument('--output',type=Path,default=ROOT/'linear')
    parser.add_argument('--smoke-only',action='store_true');a=parser.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    interface=load('linear_rna_interface',a.source/'rna_interface.py');checks=smoke(interface);write_json(a.output/'SMOKE.json',checks)
    if a.smoke_only:print(json.dumps(checks));return
    start=time.monotonic();results=[run_fold(f,interface,a.source,a.objectives,a.output) for f in range(5)]
    common=results[0][0];combined={k:np.full((55,len(common['genes'])),np.nan) for k in ('pred_effect_linear_rna','pred_effect_linear_choice')}
    assignment=np.full(55,-1,int)
    for fold,(c,q,p,meta) in enumerate(results):
        assert np.all(assignment[q]==-1);assignment[q]=fold
        for k in combined:combined[k][q]=p[k]
    assert np.all(assignment>=0) and all(np.isfinite(v).all() for v in combined.values())
    np.savez_compressed(a.output/'predictions.npz',conditions=common['all_pair_names'],genes=common['genes'],fold=assignment,control=common['control'],**combined)
    write_json(a.output/'TRAINING_RUN.json',dict(status='FROZEN_BEFORE_EXTERNAL_SCORING',code_sha256=sha(__file__),seconds=time.monotonic()-start,
        outer_query_outcomes_used=False,operations=checks['operations'],fold_count=5,variant_count=2,
        merged_prediction_sha256=sha(a.output/'predictions.npz'),fold_freeze_sha256={str(f):sha(a.output/f'fold_{f}'/'FREEZE.json') for f in range(5)}))
    print(json.dumps(dict(status='FROZEN',output=str(a.output),seconds=time.monotonic()-start)))

if __name__=='__main__':main()
