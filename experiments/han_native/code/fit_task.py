from common import *
import sys,time,argparse,traceback
sys.path.insert(0,str(Path(__file__).parent/'vendor'))
import spectral_models as sm
import baselines

FRACTIONS=[.3,.1,.03,.01];NUS=[.1,1.]
def key(family,t,n):return f'{family}_t{t:g}_n{n:g}'
def model_fit(family,n,edges,r,tau,nu,path,initial=None):
    # Candidate stems contain decimal penalties. Append, never replace suffixes.
    jf=Path(str(path)+'.json');nf=Path(str(path)+'.npz')
    if jf.exists() and nf.exists():
        m=json.loads(jf.read_text());m.update(loadnpz(nf));return m
    fun=sm.fit_rank1 if family=='rank1' else sm.fit_spectral
    m=fun(n,edges,r,tau,nu,max_iter=3000,tol=1e-6,deadline=None,initial=initial)
    arrays={k:v for k,v in m.items() if isinstance(v,np.ndarray)};meta={k:v for k,v in m.items() if not isinstance(v,np.ndarray)}
    savez(nf,**arrays);dump(jf,meta);return m
def cp(G,tr,qu,n):
    supported=np.bincount(tr.ravel(),minlength=n)>0;v=G[qu[:,0],qu[:,1]].copy();v[~(supported[qu[:,0]]&supported[qu[:,1]])]=0;return v

def run(task):
    panel,seed,fold,baseline=[task[x] for x in ('panel','seed','fold','baseline')]
    out=ROOT/'runs'/panel/baseline/str(seed)/str(fold);out.mkdir(parents=True,exist_ok=True)
    if (out/'LOCK.json').exists():return {'status':'already_locked',**task}
    begun=time.monotonic();started=now();d=loadnpz(ROOT/'prepared'/panel/str(seed)/str(fold)/'training.npz')
    ed,y,single,qu=d['edges'],d['y'],d['singles'],d['query_edges'];n=len(single);folds=d['inner_fold'];base_oof=np.zeros_like(y);grid={};probes=[]
    for fam in ('rank1','spectral'):
        for nu in NUS:
            for f in FRACTIONS:grid[key(fam,f,nu)]=dict(family=fam,fraction=f,nu=nu,pred=np.zeros_like(y),eligible=True)
    for k in range(3):
        tr=folds!=k;va=~tr;st=baselines.fit(baseline,single,ed[tr],y[tr],d['pairs'][tr],seed,out/f'inner{k}_baseline.json')
        base_oof[va]=baselines.predict(st,single,ed[va]);r=y[tr]-baselines.predict(st,single,ed[tr]);tz=sm.tau_zero(n,ed[tr],r)
        for nu in NUS:
            warm=None
            for f in FRACTIONS:
                ms={}
                for fam in ('rank1','spectral'):
                    name=key(fam,f,nu);m=model_fit(fam,n,ed[tr],r,tz*f,nu,out/'inner_models'/f'{name}_f{k}',warm if fam=='spectral' else None);ms[fam]=m
                    if fam=='spectral':warm=m['G']
                    grid[name]['eligible'] &= bool(m['converged']) and not m.get('timed_out',False)
                    grid[name]['pred'][va]=base_oof[va]+cp(m['G'],ed[tr],ed[va],n)[:,None]
                    dump(out/'PROGRESS.json',dict(utc=now(),phase='inner',inner_fold=k,config=name,fit_seconds=m['seconds'],converged=m['converged'],seconds=time.monotonic()-begun))
                gap=ms['spectral']['objective']-ms['rank1']['objective'];tol=1e-7+1e-5*max(abs(ms['rank1']['objective']),abs(ms['spectral']['objective']))
                probes.append(dict(inner_fold=k,nu=nu,fraction=f,objective_gap=gap,tolerance=tol,passed=gap<=tol))
    if not all(p['passed'] for p in probes):raise RuntimeError('Native convex/nonconvex objective-order check failed')
    base_sse=float(np.sum((y-base_oof)**2));summary=[]
    for name,row in grid.items():
        row['total']=float(np.sum((y-row['pred'])**2));summary.append({k:v for k,v in row.items() if k!='pred'}|{'key':name})
    dump(out/'INNER_SELECTION_SCORES.json',dict(baseline_sse=base_sse,candidates=summary,objective_probes=probes))
    state=baselines.fit(baseline,single,ed,y,d['pairs'],seed,out/'baseline.json');bq=baselines.predict(state,single,qu);r=y-baselines.predict(state,single,ed);tz=sm.tau_zero(n,ed,r)
    pq=Project(qu,n);pin={k:Project(ed[folds==k],n) for k in range(3)}
    payload=dict(query_ids=d['query_ids'],query_edges=qu,B=bq);selection={}
    for fam in ('rank1','spectral'):
        choices=[x for x in summary if x['family']==fam and x['eligible']]
        if not choices:raise RuntimeError('No numerically eligible candidate: '+fam)
        best=min(choices,key=lambda x:(x['total'],-x['fraction'],-x['nu']));use=best['total']<base_sse-1e-12*max(1.,base_sse)
        m=model_fit(fam,n,ed,r,tz*best['fraction'],best['nu'],out/f'{fam}_final')
        if not m['converged'] or m.get('timed_out',False):raise RuntimeError('Final selected fit failed: '+fam)
        delta=np.repeat(cp(m['G'],ed,qu,n)[:,None],y.shape[1],axis=1) if use else np.zeros_like(bq)
        doof=grid[best['key']]['pred']-base_oof if use else np.zeros_like(y);e=y-base_oof;pe=np.zeros_like(e);pd=np.zeros_like(e)
        for k in range(3):
            va=folds==k;pe[va]=pin[k].p(e[va]);pd[va]=pin[k].p(doof[va])
        A=float(np.sum(pe*pd));H=float(np.sum(pd*pd));alpha=float(np.clip(A/H,0,1)) if H>1e-14 else 0.
        Af=float(np.sum(e*doof));Hf=float(np.sum(doof*doof));gamma=float(np.clip(Af/Hf,0,1)) if Hf>1e-14 else 0.
        part=pq.p(delta);payload[fam+'_D']=delta;payload[fam+'_PD']=part
        payload[fam+'_original']=bq+delta;payload[fam+'_whole_same']=bq+alpha*delta;payload[fam+'_off']=bq+delta-part;payload[fam+'_controlled']=bq+delta+(alpha-1)*part;payload[fam+'_whole_full']=bq+gamma*delta
        selection[fam]=dict(selected=best['key'] if use else 'G0',best_regularized=best['key'],inner_gain_pct=100*(1-best['total']/base_sse),alpha=alpha,gamma=gamma,partner_cross=A,partner_energy=H,full_cross=Af,full_energy=Hf,final_converged=True)
        savez(out/f'{fam}_selected_inner.npz',edges=ed,folds=folds,baseline=base_oof,D=doof)
    savez(out/'ALL_INNER_PREDICTIONS.npz',baseline=base_oof,**{k:r['pred'] for k,r in grid.items()})
    savez(out/'predictions.npz',**payload);dump(out/'SELECTION.json',selection)
    files=[out/'predictions.npz',out/'SELECTION.json',out/'baseline.json',out/'ALL_INNER_PREDICTIONS.npz',out/'INNER_SELECTION_SCORES.json',*out.glob('*_final.*'),*out.glob('*_selected_inner.npz')]
    dump(out/'LOCK.json',dict(task=task,started_utc=started,completed_utc=now(),seconds=time.monotonic()-begun,query_outcomes_opened_by_worker=False,training_sha256=digest(ROOT/'prepared'/panel/str(seed)/str(fold)/'training.npz'),protocol_sha256=digest(ROOT/'PROTOCOL.md'),source_hashes={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'code').rglob('*.py'))},files={str(p.relative_to(out)):digest(p) for p in files}))
    print(json.dumps(dict(status='locked',**task,seconds=time.monotonic()-begun)),flush=True)
    return dict(status='locked',**task)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--task',type=int,required=True);args=ap.parse_args();task=json.loads((ROOT/'TASKS.json').read_text())[args.task]
    try:run(task)
    except Exception:
        dump(ROOT/'failures'/f'task{args.task}.json',dict(task=task,utc=now(),traceback=traceback.format_exc()));raise
if __name__=='__main__':main()
