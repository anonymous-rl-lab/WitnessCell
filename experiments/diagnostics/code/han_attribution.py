"""Fixed retrospective attribution; frozen source weights, pairwise target calibration."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,argparse
import numpy as np
import pandas as pd

class Project:
    def __init__(self, edges,n):
        self.X=np.zeros((len(edges),n));self.X[np.arange(len(edges)),edges[:,0]]=1;self.X[np.arange(len(edges)),edges[:,1]]=1
        v,u=np.linalg.eigh(self.X.T@self.X);keep=v>max(self.X.shape)*np.finfo(float).eps*max(v.max(),1)
        self.inv=(u[:,keep]/v[keep])@u[:,keep].T
    def coef(self,a):return self.inv@(self.X.T@a.reshape(len(a),-1))
    def p(self,a):return (a.reshape(len(a),-1)-self.X@self.coef(a)).reshape(a.shape)

def corr(a,b):return float(np.corrcoef(a.ravel(),b.ravel())[0,1])
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def savej(p,x):p.write_text(json.dumps(x,indent=2)+'\n')

def features(singles,edges):
    n,r=singles.shape;B=singles[edges[:,0]]+singles[edges[:,1]]
    scale=np.sqrt(np.mean(singles**2,axis=0));lam=np.array([0,.05,.1,.2,.5,1,2,5])
    sat=B[:,:,None]/(1+lam[None,None,:]*abs(B[:,:,None])/scale[None,:,None])-B[:,:,None]
    t=np.tanh(singles/scale);a,b=t[edges[:,0]],t[edges[:,1]]
    ap,an=np.maximum(a,0),np.maximum(-a,0);bp,bn=np.maximum(b,0),np.maximum(-b,0)
    op=ap*bn+an*bp
    local=np.stack([ap*bp,an*bn,op,op*np.tanh(a+b)],axis=-1)*scale[None,:,None]
    return B,sat,local

def fit_residual(P,F,e,ridge):
    pf=P.p(F);pe=P.p(e);r=e.shape[1];co=[]
    for j in range(r):
        g=pf[:,j].T@pf[:,j];b=pf[:,j].T@pe[:,j];sc=np.sqrt(np.maximum(np.diag(g),0));sc=np.maximum(sc,max(sc.max()*1e-10,1e-15))
        eta=np.linalg.pinv(g/sc[:,None]/sc[None,:]+ridge*np.eye(len(sc)),rcond=1e-12)@(b/sc)
        co.append(eta/sc)
    co=np.asarray(co);raw=np.einsum('nrf,rf->nr',F,co);end=P.coef(e-raw)
    return co,end
def pred_residual(P,F,co,end):return np.einsum('nrf,rf->nr',F,co)+P.X@end

def fit_baselines(edges,pairs,n,B,Y,Sat,Loc,cal,query,seed):
    order=sorted(cal,key=lambda i:hashlib.sha256(('v37-inner|'+str(seed)+'|'+'|'.join(sorted(pairs[i]))).encode()).digest())
    blocks=np.array_split(order,3);candidates=[('saturation',i,0.) for i in range(8)]+[('local4',0,l) for l in (10.,1.,.1,.01)]
    errors=np.zeros(len(candidates))
    for va in blocks:
        tr=np.setdiff1d(cal,va);pt,pv=Project(edges[tr],n),Project(edges[va],n)
        for k,(kind,i,reg) in enumerate(candidates):
            F=Sat[:,:,i:i+1] if kind=='saturation' else Loc
            co,en=fit_residual(pt,F[tr],Y[tr]-B[tr],reg)
            pred=B[va]+pred_residual(pv,F[va],co,en)
            errors[k]+=np.sum((Y[va]-pred)**2)
    pc,pq=Project(edges[cal],n),Project(edges[query],n);out={};choices=[]
    for kind in ('saturation','local4'):
        ix=[i for i,c in enumerate(candidates) if c[0]==kind];best=min(ix,key=lambda j:errors[j]);_,i,reg=candidates[best]
        F=Sat[:,:,i:i+1] if kind=='saturation' else Loc
        co,en=fit_residual(pc,F[cal],Y[cal]-B[cal],reg)
        out[kind]=(B[cal]+pred_residual(pc,F[cal],co,en),B[query]+pred_residual(pq,F[query],co,en))
        choices.append(dict(kind=kind,saturation_lambda=[0,.05,.1,.2,.5,1,2,5][i] if kind=='saturation' else None,ridge=reg,inner_sse=float(errors[best]),coefficients=co.tolist(),endpoint_coefficients=en.tolist()))
    return out,choices,pc,pq

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--han',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--draws',type=int,default=999);args=ap.parse_args();H=args.han.resolve();O=args.out;O.mkdir(parents=True,exist_ok=True)
    prepared=(H/'PREPARED.json').exists()
    if prepared:
        cache=np.load(H/'FROZEN_SF12_FEATURES.npz');ci=json.loads((H/'SOURCE_COEFFICIENTS.json').read_text());sources={name:None for name in ('HeLa','293T')}
    else:
        sys.path.insert(0,str(H/'frozen_engine/code'));import predict as engine
        sources={name:engine.Source(name) for name in ('HeLa','293T')}
        K=np.load(H/'data/target_knowledge.npz');kl={g:i for i,g in enumerate(K['genes'])}
    randoms={};coefinfo={};source_hashes={str(p.relative_to(H)):digest(p) for p in ((H/'frozen_engine/pretrained').rglob('*') if not prepared else H.rglob('*')) if p.is_file()}
    rng=np.random.default_rng(20261109)
    for name,s in sources.items():
        if prepared:theta=np.array(ci[name]['theta']);scale=np.array(ci[name]['source_feature_scale'])
        else:
            theta=s.sf_state[s.phase_sel['sf12_selected']['key']+'__theta'];f=engine.sf_features(s.smean,s.rms,s.x,s.center,s.spread,s.train_edges)
            pf=Project(s.train_edges,s.n).p(f);scale=np.sqrt(np.sum(pf**2,axis=0));scale=np.maximum(scale,max(scale.max()*1e-10,1e-15))
        u=rng.normal(size=(12,args.draws));u/=np.linalg.norm(u,axis=0);u*=np.linalg.norm(theta*scale)
        randoms[name]=u/scale[:,None];coefinfo[name]=dict(theta=theta.tolist(),source_feature_scale=scale.tolist())
    np.savez_compressed(O/'RANDOM_COEFFICIENTS.npz',**randoms);savej(O/'SOURCE_COEFFICIENTS.json',coefinfo)
    allscores=[];randomrows=[];oracle=[];qrows=[];checks=[];choicesall=[];featurecache={}
    for panel in ('full_raw','full_published'):
        inp=np.load(H/f'data/{panel}_inputs.npz');frozen=np.load(H/f'data/{panel}_frozen_predictions.npz');Y=np.load(H/f'data/{panel}_outcomes.npz')['outcomes'];pairs=inp['pairs'];genes=inp['genes'];n=len(genes);lookup={g:i for i,g in enumerate(genes)};edges=np.array([[lookup[a],lookup[b]] for a,b in pairs]);B,Sat,Loc=features(inp['singles'],edges)
        assert np.max(abs(B-frozen['B']))<1e-12
        Fsrc={};theta={};Dsrc={}
        for name,s in sources.items():
            if prepared:Fsrc[name]=cache[panel+'__'+name]
            else:
                source_mean=s.knowledge.astype(float).mean(0);w=np.linalg.lstsq(s.knowledge.astype(float)-source_mean,s.x,rcond=None)[0];kt=np.array([K['embedding'][kl[g]] for g in genes]);z=(kt-source_mean)@w
                beta=s.meta['sf12']['native_beta'];Fs=[]
                for r in range(Y.shape[1]):
                    scale=float(np.std(inp['singles'][:,r]));su=inp['singles'][:,r]*s.h/scale
                    Fs.append(engine.sf_features(su,s.rms,z,s.center,s.spread,edges)*beta*scale/s.h)
                Fsrc[name]=np.stack(Fs,axis=1)
            theta[name]=np.array(coefinfo[name]['theta']);Dsrc[name]=np.einsum('nrf,f->nr',Fsrc[name],theta[name]);m=list(frozen['models']).index(name+'__sf12')
            err=float(np.max(abs(Dsrc[name]-frozen['D'][:,:,m])));assert err<1e-10;checks.append(dict(check='source_direction_replay',panel=panel,source=name,max_abs=err));featurecache[panel+'__'+name]=Fsrc[name]
        for p in sorted((H/'data/splits'/panel).glob('*/*/identities.json')):
            ids=json.loads(p.read_text());cal=np.array(ids['train']);query=np.array(ids['query']);seed=int(p.parent.parent.name);fold=int(p.parent.name)
            assert not set(cal)&set(query);assert len(cal)+len(query)==len(pairs)
            bases,ch,pc,pq=fit_baselines(edges,pairs,n,B,Y,Sat,Loc,cal,query,seed)
            bases={'additive':(B[cal],B[query]),**bases}
            # Reproduce the original two calibration projection blocks exactly.
            order=sorted(cal,key=lambda i:hashlib.sha256((str(seed+1000+37*fold)+'|'+'|'.join(sorted(pairs[i]))).encode()).hexdigest())
            positions={i:j for j,i in enumerate(cal)}
            cb=[(np.array([positions[i] for i in block]),Project(edges[block],n)) for block in np.array_split(order,2)]
            def calibrate_project(a):
                out=np.empty_like(a)
                for loc,pr in cb:out[loc]=pr.p(a[loc])
                return out
            for v in ch:choicesall.append(dict(panel=panel,seed=seed,fold=fold,**v))
            if Y.shape[1]==2:
                qrows.append(dict(panel=panel,seed=seed,fold=fold,training_pairs=len(cal),query_pairs=len(query),training_q_raw=corr(*pc.p(Y[cal]).T),query_repeat_r=corr(*pq.p(Y[query]).T)))
            for baseline,(bc,bq) in bases.items():
                erc=Y[cal]-bc;erq=Y[query]-bq;rc,rq=calibrate_project(erc),pq.p(erq);E=float(np.sum(rq*rq));EF=float(np.sum(erq*erq))
                allscores.append(dict(panel=panel,seed=seed,fold=fold,baseline=baseline,source='none',policy='baseline',alpha=0.,partner_sse=E,full_sse=EF,base_partner_sse=E,base_full_sse=EF))
                for name in sources:
                    fc,fq=calibrate_project(Fsrc[name][cal]),pq.p(Fsrc[name][query]);th=theta[name];dc=np.einsum('nrf,f->nr',fc,th);dq=np.einsum('nrf,f->nr',fq,th)
                    A=float(np.sum(rc*dc));BB=float(np.sum(dc*dc));a=float(np.clip(A/BB,0,1)) if BB>1e-14 else 0.
                    for pol,aa in [('original',1.),('controlled',a)]:
                        sp=float(np.sum((rq-aa*dq)**2));sf=float(np.sum((erq-aa*dq)**2))
                        allscores.append(dict(panel=panel,seed=seed,fold=fold,baseline=baseline,source=name,policy=pol,alpha=aa,partner_sse=sp,full_sse=sf,base_partner_sse=E,base_full_sse=EF))
                        assert abs((EF-sf)-(E-sp))<1e-8
                    if baseline!='additive':continue
                    original_params=json.loads((H/'results'/panel/str(seed)/str(fold)/'PARAMETERS.json').read_text());m=list(frozen['models']).index(name+'__sf12');assert abs(a-original_params['alpha'][m])<1e-10
                    aq=float(np.sum(rq*dq));bqenergy=float(np.sum(dq*dq));au=aq/bqenergy;ab=float(np.clip(au,0,1));an=max(0.,au)
                    oracle.append(dict(panel=panel,seed=seed,fold=fold,source=name,E=E,A=aq,B=bqenergy,actual_alpha=a,unconstrained_alpha=au,bounded_alpha=ab,nonnegative_alpha=an,actual_improvement=2*a*aq-a*a*bqenergy,bounded_improvement=2*ab*aq-ab*ab*bqenergy,nonnegative_improvement=2*an*aq-an*an*bqenergy))
                    # Gram arithmetic avoids materializing 999 complete prediction panels.
                    c=fc.reshape(-1,12);q=fq.reshape(-1,12);G=c.T@c;h=c.T@rc.ravel();Gq=q.T@q;hq=q.T@rq.ravel();U=randoms[name]
                    en=np.einsum('ij,ij->j',U,G@U);mult=np.sqrt(BB/np.maximum(en,1e-300));V=U*mult
                    ar=np.clip((V.T@h)/BB,0,1);arq=V.T@hq;brq=np.einsum('ij,ij->j',V,Gq@V);improve=2*ar*arq-ar*ar*brq
                    if fold==0:
                        for j in (0,args.draws-1):
                            direct=q@(V[:,j]*ar[j]);sse=np.sum((rq.ravel()-direct)**2);assert abs(sse-(E-improve[j]))<1e-8
                    for j in range(args.draws):randomrows.append(dict(panel=panel,seed=seed,fold=fold,source=name,draw=j,alpha=ar[j],scale=mult[j],improvement=float(improve[j]),base_partner_sse=E))
            print(f'{panel} {seed} {fold} complete',flush=True)
    scores=pd.DataFrame(allscores);scores.to_csv(O/'SCORES.csv',index=False)
    groups=['panel','baseline','source','policy'];sums=['partner_sse','full_sse','base_partner_sse','base_full_sse']
    pooled=scores.groupby(groups)[sums].sum().reset_index()
    for metric in ['partner','full']:pooled[metric+'_gain_pct']=100*(1-pooled[metric+'_sse']/pooled['base_'+metric+'_sse'])
    pooled.to_csv(O/'POOLED_SCORES.csv',index=False)
    rnd=pd.DataFrame(randomrows);rnd.to_csv(O/'RANDOM_FOLD_SCORES.csv.gz',index=False)
    rp=rnd.groupby(['panel','source','draw'])[['improvement','base_partner_sse']].sum().reset_index();rp['gain_pct']=100*rp.improvement/rp.base_partner_sse;rp.to_csv(O/'RANDOM_POOLED.csv',index=False)
    summaries=[]
    for (panel,name),t in rp.groupby(['panel','source']):
        obs=float(pooled[(pooled.panel==panel)&(pooled.source==name)&(pooled.baseline=='additive')&(pooled.policy=='controlled')].partner_gain_pct.iloc[0]);v=t.gain_pct.to_numpy();exc=int((v>=obs).sum())
        summaries.append(dict(panel=panel,source=name,observed_gain_pct=obs,draws=len(v),random_median=float(np.median(v)),random_q025=float(np.quantile(v,.025)),random_q95=float(np.quantile(v,.95)),random_q975=float(np.quantile(v,.975)),random_max=float(v.max()),exceedances=exc,empirical_tail_fraction=(1+exc)/(1+len(v))))
    pd.DataFrame(summaries).to_csv(O/'RANDOM_SUMMARY.csv',index=False)
    pd.DataFrame(oracle).to_csv(O/'HAN_ORACLE_FOLDS.csv',index=False);pd.DataFrame(qrows).to_csv(O/'HAN_REPEATABILITY.csv',index=False)
    savej(O/'BASELINE_SELECTIONS.json',choicesall);savej(O/'CHECKS.json',checks)
    np.savez_compressed(O/'FROZEN_SF12_FEATURES.npz',**featurecache)
    assert source_hashes=={str(p.relative_to(H)):digest(p) for p in ((H/'frozen_engine/pretrained').rglob('*') if not prepared else H.rglob('*')) if p.is_file()}
    savej(O/'RUN_RECORD.json',dict(status='PASS',source_predictor_fits=0,random_draws_per_source=args.draws,folds=30,source_hashes=source_hashes,plan_sha256=digest(Path(__file__).parents[1]/'ANALYSIS_PLAN.md'),script_sha256=digest(__file__),status_scope='retrospective attribution; comparator baselines and scalar calibration fitted on permitted target training pairs'))
    print(pd.DataFrame(summaries).to_string(index=False),flush=True)
    print(pooled[pooled.policy.isin(['baseline','controlled'])].to_string(index=False),flush=True)

if __name__=='__main__':main()
