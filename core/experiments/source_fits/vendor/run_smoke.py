"""Bounded, one-fold two-budget development; evaluation labels read only by score.py."""
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_k]='1'
from pathlib import Path
import json,time,hashlib,datetime
import numpy as np
from frozen_models import fit_model
from spectral_models import fit_spectral,fit_rank1,tau_zero

ROOT=Path(__file__).resolve().parent
def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    return x
def dump(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(clean(x),indent=2,allow_nan=False))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def baseline(singles,edges,y,deadline):
    scales=np.std(singles,axis=0);scales[scales==0]=1
    models=[fit_model(singles[:,t]/scales[t],edges,y[:,t]/scales[t],rank=0,alpha=10.,deadline=deadline) for t in range(2)]
    def predict(ed):return np.column_stack([models[t].predict(ed)*scales[t] for t in range(2)])
    return models,predict
def corr_prediction(G,train_edges,query_edges,n):
    supported=np.bincount(train_edges.ravel(),minlength=n)>0
    pred=G[query_edges[:,0],query_edges[:,1]].copy()
    pred[~(supported[query_edges[:,0]]&supported[query_edges[:,1]])]=0
    return pred
def diag(model):return {k:v for k,v in model.items() if k!='G'}
def fit_one(family,n,ed,r,tau,nu,plan,deadline,initial=None):
    fun=fit_spectral if family=='spectral' else fit_rank1
    return fun(n,ed,r,tau,nu,max_iter=plan['optimizer_maxiter'],tol=plan['optimizer_tolerance'],deadline=deadline,initial=initial)
def run_budget(name,data,query_edges,query_ids,plan,deadline):
    out=ROOT/'results'/name;out.mkdir(parents=True,exist_ok=True)
    pos=data['low_positions'] if name=='low' else np.arange(len(data['edges']))
    edges,y,folds=data['edges'][pos],data['y'][pos],data['inner_fold'][pos]
    n=len(data['genes']);m=len(edges)
    base_oof=np.zeros_like(y);grid={};probe=[]
    configs=[(fraction,nu) for nu in plan['nu'] for fraction in plan['tau_fractions']]
    for family in ('rank1','spectral'):
        for fraction,nu in configs:
            key=f'{family}_t{fraction:g}_n{nu:g}'
            grid[key]={'family':family,'fraction':fraction,'nu':nu,'pred':np.zeros_like(y),'diagnostics':[],'eligible':True}
    for f in range(3):
        train=folds!=f;val=~train
        models,predict=baseline(data['singles'],edges[train],y[train],deadline)
        base_oof[val]=predict(edges[val]);r=y[train]-predict(edges[train])
        dump(out/f'inner{f}_baseline.json',[x.to_dict() for x in models])
        tz=float(tau_zero(n,edges[train],r))
        for nu in plan['nu']:
            warm=None
            for fraction in plan['tau_fractions']:
                pair={}
                for family in ('rank1','spectral'):
                    key=f'{family}_t{fraction:g}_n{nu:g}';row=grid[key]
                    model=fit_one(family,n,edges[train],r,tz*fraction,nu,plan,deadline,warm if family=='spectral' else None)
                    pair[family]=model
                    row['eligible'] &= bool(model['converged'])
                    row['diagnostics'].append({'inner_fold':f,'tau_zero':tz,'tau':tz*fraction,**diag(model)})
                    if family=='spectral':warm=model['G']
                    correction=corr_prediction(model['G'],edges[train],edges[val],n)
                    row['pred'][val]=base_oof[val]+correction[:,None]
                    np.savez_compressed(out/'inner_models'/f'{key}_f{f}.npz',G=model['G']) if (out/'inner_models').exists() else None
                gap=float(pair['spectral']['objective']-pair['rank1']['objective'])
                tolerance=1e-7+1e-5*max(abs(pair['rank1']['objective']),abs(pair['spectral']['objective']))
                probe.append({'inner_fold':f,'fraction':fraction,'nu':nu,'spectral_minus_rank1_objective':gap,'tolerance':tolerance,'pass':gap<=tolerance})
                dump(out/'progress.json',{'phase':'inner','fold':f,'fraction':fraction,'nu':nu,'objective_probe':probe[-1]})
                print(json.dumps({'budget':name,'phase':'inner','fold':f,'fraction':fraction,'nu':nu,
                    'seconds':[pair[k].get('seconds') for k in ('rank1','spectral')],
                    'converged':[pair[k]['converged'] for k in ('rank1','spectral')],
                    'objective_order_ok':probe[-1]['pass']}),flush=True)
    base_sse=np.sum((y-base_oof)**2,axis=0)
    rows=[]
    for key,row in grid.items():
        pred=row.pop('pred');ss=np.sum((y-pred)**2,axis=0)
        row.update({'key':key,'sse':ss,'total':float(ss.sum()),'gain_each':1-ss/base_sse})
        rows.append(row)
        np.savez_compressed(out/f'{key}_inner_predictions.npz',edges=edges,fold=folds,y=y,pred=pred,baseline=base_oof)
    dump(out/'inner_summary.json',{'M0_sse':base_sse,'grid':rows,'objective_probes':probe})
    if not all(x['pass'] for x in probe):raise RuntimeError('Convex objective order probe failed; stop before scoring')
    models,predict=baseline(data['singles'],edges,y,deadline)
    dump(out/'M0_models.json',[x.to_dict() for x in models])
    b=predict(query_edges);r=y-predict(edges);tz=float(tau_zero(n,edges,r))
    payload={'edges':query_edges,'global_ids':query_ids,'M0':b,'M0_repeat_mean':np.repeat(b.mean(axis=1)[:,None],2,axis=1)}
    selections={}
    for family in ('rank1','spectral'):
        candidates=[x for x in rows if x['family']==family and x['eligible']]
        if not candidates:raise RuntimeError('No numerically eligible candidate in '+family)
        best=min(candidates,key=lambda x:(x['total'],-x['fraction'],-x['nu']))
        use=best['total']<float(base_sse.sum())-1e-12*max(1.,float(base_sse.sum()))
        model=fit_one(family,n,edges,r,tz*best['fraction'],best['nu'],plan,deadline)
        if not model['converged']:raise RuntimeError('Selected final fit not converged: '+family)
        correction=corr_prediction(model['G'],edges,query_edges,n)
        payload[family+'_regularized']=b+correction[:,None]
        payload[family+'_selected']=payload[family+'_regularized'] if use else b.copy()
        payload[family+'_G']=model['G']
        dump(out/(family+'_model.json'),{'tau':tz*best['fraction'],'nu':best['nu'],**diag(model)})
        selections[family]={'selected':best['key'] if use else 'G0','best_regularized':best['key'],
            'inner_gain_each':best['gain_each'],'inner_gain_pooled':1-best['total']/float(base_sse.sum()),
            'regularized_result_role':'chosen predictor' if use else 'diagnostic only; not selected',
            'final_converged':True}
    np.savez_compressed(out/'predictions.npz',**payload)
    dump(out/'selection.json',selections)
    dump(out/'split.json',{'train_global_ids':data['global_ids'][pos],'evaluation_global_ids':query_ids,
        'inner_fold':folds,'train_count':m,'inner_sizes':np.bincount(folds,minlength=3)})
    print(json.dumps(clean({'budget':name,'phase':'refit_complete_no_outer_scoring','selection':selections})),flush=True)
def main():
    plan=json.loads((ROOT/'PLAN.json').read_text());start=time.monotonic();deadline=start+plan['real_fit_budget_seconds']
    with np.load(ROOT/'inputs/training.npz',allow_pickle=False) as z:data={k:z[k] for k in z.files}
    with np.load(ROOT/'inputs/evaluation.npz',allow_pickle=False) as z:query_edges=z['edges'];query_ids=z['global_ids']
    assert not(set(data['global_ids'])&set(query_ids))
    for name in ('low','high'):(ROOT/'results'/name/'inner_models').mkdir(parents=True,exist_ok=True)
    try:
        for name in ('low','high'):run_budget(name,data,query_edges,query_ids,plan,deadline)
        paths=[p for name in ('low','high') for p in [ROOT/'results'/name/'predictions.npz',ROOT/'results'/name/'selection.json']]
        dump(ROOT/'results/PREDICTIONS_LOCKED.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'evaluation_labels_scored':False,'files':{str(p.relative_to(ROOT)):digest(p) for p in paths},
            'code_sha256':{p.name:digest(p) for p in ROOT.glob('*.py')},'plan_sha256':digest(ROOT/'PLAN.json')})
    finally:
        dump(ROOT/'results/fit_timing.json',{'seconds':time.monotonic()-start,'limit_seconds':plan['real_fit_budget_seconds']})
if __name__=='__main__':main()
