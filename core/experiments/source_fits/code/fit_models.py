"""Finite functional adaptation of Next2 partner projection/calibration.

Reads training responses and query identities only. Never imports outer truth.
All reported selection losses are development cross-validation quantities.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[k]='1'
import argparse, ast, hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'vendor'))
import run_smoke as native
import spectral_models as spectral

# Use the exact recovered Next2 constrained scalar optimizer, without importing
# its RNA-only data interface. Preserve the source and extract its AST verbatim.
tree=ast.parse((ROOT/'vendor/next2_baseline_core.py').read_text())
node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='solve_beta')
space={'np':np,'CONFIG':{'numerical_relative_tolerance':float(64*np.finfo(float).eps)}}
exec(compile(ast.Module(body=[node],type_ignores=[]),'next2_solve_beta','exec'),space)
solve_beta=space['solve_beta']


def dump(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(native.clean(data),indent=2,allow_nan=False)+'\n');tmp.replace(path)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Projection:
    def __init__(self,edges,n):
        self.edges=np.asarray(edges,int);self.n=n;self.m=len(edges)
        self.X=np.zeros((self.m,n));self.X[np.arange(self.m)[:,None],self.edges]=1
        self.inv=np.linalg.pinv(self.X.T@self.X, rcond=1e-12)
    def coeff(self,value):
        v=np.asarray(value);out=np.zeros((self.n,)+v.shape[1:])
        np.add.at(out,self.edges[:,0],v);np.add.at(out,self.edges[:,1],v)
        return self.inv@out
    def apply(self,value):
        b=self.coeff(value)
        return value-b[self.edges[:,0]]-b[self.edges[:,1]]
    def extend(self,value,edges):
        b=self.coeff(value);return b[edges[:,0]]+b[edges[:,1]]


def agreement(z):
    return float(np.clip(np.sum(z[:,0]*z[:,1])/(.5*np.sum(z*z)+1e-300),0,1))


def fit_projected_rank1(n,edges,resid,tau,nu,projection,maxiter=1500,tol=1e-6):
    r=projection.apply(np.mean(resid,axis=1));i,j=edges.T;m=len(edges)
    scale=max(float(np.mean(r*r)),1e-15)
    mat=np.zeros((n,n));mat[i,j]=r;mat[j,i]=r
    ev,V=np.linalg.eigh(mat);order=np.argsort(-np.abs(ev))[:4]
    rng=np.random.default_rng(20261018)
    starts=[V[:,k] for k in order]+[rng.normal(size=n) for _ in range(2)]
    def objective(w):
        norm=np.linalg.norm(w);u=w/norm;v=projection.apply(u[i]*u[j])
        a=np.mean(v*v)+nu/n**2;b=np.mean(r*v)
        lam=np.sign(b)*max(abs(b)-tau/(2*n),0)/a if a>0 else 0.
        error=lam*v-r
        val=np.mean(error*error)+tau/n*abs(lam)+nu/n**2*lam*lam
        coeff=2*lam*projection.apply(error)/m
        gu=np.bincount(i,weights=coeff*u[j],minlength=n)+np.bincount(j,weights=coeff*u[i],minlength=n)
        grad=(gu-u*np.dot(u,gu))/norm
        return float(val/scale),grad/scale,u,float(lam)
    best=(float(np.mean(r*r)/scale),np.zeros(n),0.,True,0.)
    diagnostics=[]
    for start in starts:
        fit=minimize(lambda w:objective(w)[:2],start,jac=True,method='L-BFGS-B',
                     options={'maxiter':maxiter,'ftol':1e-14,'gtol':1e-8,'maxls':50})
        val,grad,u,lam=objective(fit.x);gn=float(np.linalg.norm(grad)*np.linalg.norm(fit.x))
        good=gn<=tol
        diagnostics.append({'objective':val*scale,'stationarity_relative':gn,'converged':good,'nit':int(fit.nit)})
        if good and val<best[0]:best=(val,u,lam,good,gn)
    _,u,lam,good,gn=best
    return dict(G=lam*np.outer(u,u),converged=good,stationarity_relative=gn,starts=diagnostics)


def symfeatures(f,edges):
    a,b=f[edges[:,0]],f[edges[:,1]];d=f.shape[1]
    return np.column_stack([a[:,i]*b[:,i] if i==j else (a[:,i]*b[:,j]+a[:,j]*b[:,i])/np.sqrt(2)
                            for i in range(d) for j in range(i,d)])


def features(family,singles,edges,resid,query,n,proj):
    all_edges=np.vstack([edges,query]);s=singles.mean(axis=1)
    z=s/max(np.sqrt(np.mean(s*s)),1e-12)
    if family=='native_product':
        f=singles/max(np.sqrt(np.mean(singles*singles)),1e-12)
        F=symfeatures(f,all_edges)
    elif family=='next2_nonlinear':
        a,b=z[all_edges[:,0]],z[all_edges[:,1]]
        cols=[]
        for c in (.5,1.,2.):
            aa,bb=c*a,c*b
            cols.extend([np.tanh(aa+bb)-np.tanh(aa)-np.tanh(bb),
                (aa+bb)/(1+np.exp(-np.clip(aa+bb,-50,50)))-aa/(1+np.exp(-np.clip(aa,-50,50)))-bb/(1+np.exp(-np.clip(bb,-50,50)))])
        cols.extend([a*b,np.maximum(a,0)*np.maximum(b,0),np.maximum(-a,0)*np.maximum(-b,0),
                     np.maximum(a,0)*np.maximum(-b,0)+np.maximum(-a,0)*np.maximum(b,0)])
        F=np.column_stack(cols)
    elif family=='single_rbf':
        centers=np.quantile(z,np.linspace(0,1,8));width=max(np.std(z),.1)
        f=np.exp(-.5*((z[:,None]-centers)/width)**2)
        F=symfeatures(f,all_edges)
    else:
        rr=proj.apply(resid);mats=[]
        degree=np.bincount(edges.ravel(),minlength=n)
        for k in range(2):
            mat=np.zeros((n,n));mat[edges[:,0],edges[:,1]]=rr[:,k];mat[edges[:,1],edges[:,0]]=rr[:,k]
            # Identity-only exposure normalization, refitted for every inner split.
            mat=mat/np.sqrt(np.maximum(degree[:,None]*degree[None,:],1))
            mats.append(mat)
        if 'consensus' in family:C=(mats[0]@mats[1].T+mats[1]@mats[0].T)/2
        else:
            mean=(mats[0]+mats[1])/2;C=mean@mean.T
        vals,V=np.linalg.eigh((C+C.T)/2);k=int(family[-1]);ix=np.argsort(vals)[::-1][:k]
        f=V[:,ix]*np.sqrt(np.maximum(vals[ix],0))
        f=f/max(np.sqrt(np.mean(f*f)),1e-12)
        F=symfeatures(f,all_edges)
    Ft,Fq=F[:len(edges)],F[len(edges):]
    b=proj.coeff(Ft);Ft=proj.apply(Ft);Fq=Fq-b[query[:,0]]-b[query[:,1]]
    trace=np.sum(Ft*Ft);scale=np.sqrt(trace/len(edges)) if trace>1e-24 else 1.
    return Ft/scale,Fq/scale


def one_fit(d,idx,query,plan):
    n=len(d['genes']);edges=d['edges'][idx];y=d['y'][idx];singles=d['singles']
    _,predict=native.baseline(singles,edges,y,time.monotonic()+1800)
    base=predict(query);resid=y-predict(edges);proj=Projection(edges,n)
    w=agreement(proj.apply(y));out={};audit={}
    for family in plan['feature_families']:
        Ft,Fq=features(family,singles,edges,resid,query,n,proj)
        gram=Ft.T@Ft;rhs=Ft.T@proj.apply(resid.mean(axis=1))
        vals,V=np.linalg.eigh((gram+gram.T)/2)
        for lam in plan['ridge']:
            theta=V@((V.T@rhs)/(np.maximum(vals,0)+lam))
            corr=Fq@theta
            out[f'{family}_l{lam:g}']=(np.zeros_like(base),np.repeat(corr[:,None],2,axis=1))
    tz=spectral.tau_zero(n,edges,resid)
    pz=spectral.tau_zero(n,edges,proj.apply(resid))
    for nu in plan['rank1_nu']:
        for frac in plan['rank1_tau_fractions']:
            key=f'native_decomp_t{frac:g}_n{nu:g}'
            model=spectral.fit_rank1(n,edges,resid,tz*frac,nu,max_iter=3000,tol=1e-6)
            G=model['G'];train_g=G[edges[:,0],edges[:,1]];end=proj.extend(train_g,query)
            corr=G[query[:,0],query[:,1]]-end
            out[key]=(np.repeat(end[:,None],2,axis=1),np.repeat(corr[:,None],2,axis=1))
            audit[key]={'converged':model['converged'],'stationarity_relative':model['stationarity_relative']}
            key=f'projected_rank1_t{frac:g}_n{nu:g}'
            model=fit_projected_rank1(n,edges,resid,pz*frac,nu,proj)
            G=model['G'];train_g=G[edges[:,0],edges[:,1]]
            corr=G[query[:,0],query[:,1]]-proj.extend(train_g,query)
            out[key]=(np.zeros_like(base),np.repeat(corr[:,None],2,axis=1))
            audit[key]={k:v for k,v in model.items() if k!='G'}
    return base,out,w,audit


def stats(base,y,add,corr,proj):
    e=base+add-y;pe=proj.apply(e);pc=proj.apply(corr)
    return np.array([np.sum(pe*pe),np.sum(pe*pc),np.sum(pc*pc),
                     np.sum(e*e),np.sum(e*corr),np.sum(corr*corr)])


def load_input(path):
    v=json.loads(path.read_text());genes=v['genes'];lookup={g:i for i,g in enumerate(genes)}
    edges=np.array([[lookup[g] for g in r['pair']] for r in v['training']]);ids=['|'.join(sorted(r['pair'])) for r in v['training']]
    queries=np.array([[lookup[g] for g in p] for p in v['queries']]);qids=['|'.join(sorted(p)) for p in v['queries']]
    assert not(set(ids)&set(qids))
    folds=np.array([int.from_bytes(hashlib.sha256(('wc_spectral_smoke_v1|'+s).encode()).digest()[:8],'big')%3 for s in ids])
    return dict(genes=np.array(genes),singles=np.array(v['singles']),edges=edges,y=np.array([r['response'] for r in v['training']]),fold=folds,query=queries,ids=np.array(ids),query_ids=np.array(qids))


def run(cell):
    plan=json.loads((ROOT/'PLAN.json').read_text());out=ROOT/'results'/cell
    if (out/'FREEZE.json').exists():print(cell,'already frozen');return
    out.mkdir(exist_ok=True);d=load_input(ROOT/'inputs'/f'{cell}.json');n=len(d['genes']);m=len(d['edges'])
    acc={};cache=[];base_sse=0.;base_partner=0.;elig={};begin=time.monotonic()
    for f in range(3):
        tr=np.flatnonzero(d['fold']!=f);va=np.flatnonzero(d['fold']==f);q=d['edges'][va];y=d['y'][va]
        base,candidates,w,audit=one_fit(d,tr,q,plan);P=Projection(q,n)
        base_sse+=np.sum((base-y)**2);base_partner+=np.sum(P.apply(base-y)**2)
        np.savez_compressed(out/f'inner_{f}.npz',edges=q,ids=d['ids'][va],y=y,base=base,
            **{k+'__endpoint':a for k,(a,c) in candidates.items()},**{k+'__partner':c for k,(a,c) in candidates.items()})
        dump(out/f'inner_{f}_audit.json',dict(training_ids=d['ids'][tr],validation_ids=d['ids'][va],reliability=w,models=audit))
        for k,(a,c) in candidates.items():
            acc.setdefault(k,np.zeros(6));acc[k]+=stats(base,y,a,c,P)
            elig[k]=elig.get(k,True) and audit.get(k,{}).get('converged',True)
        cache.append((base,y,candidates,P,w))
        print(json.dumps({'cell':cell,'phase':'inner_complete','fold':f,'seconds':time.monotonic()-begin}),flush=True)
    rows=[]
    for k,st in acc.items():
        s=solve_beta(st,float(base_sse));s['eligible']=s['eligible'] and elig[k]
        rows.append(dict(key=k,**s))
    zero=dict(key='M0',beta=0.,eligible=True,partner_sse=float(base_partner),rna_sse=float(base_sse))
    chosen=min([zero]+[r for r in rows if r['eligible']],key=lambda r:(r['partner_sse'],r['rna_sse'],r['key']))
    groups={'selected':chosen,'next2_selected':min([zero]+[r for r in rows if r['eligible'] and not any(r['key'].startswith(p) for p in ['native_decomp','projected_rank1'])],key=lambda r:(r['partner_sse'],r['rna_sse'],r['key'])),
        'projected_selected':min([zero]+[r for r in rows if r['eligible'] and r['key'].startswith('projected_rank1')],key=lambda r:(r['partner_sse'],r['rna_sse'],r['key'])),
        'native_decomp_selected':min([zero]+[r for r in rows if r['eligible'] and r['key'].startswith('native_decomp')],key=lambda r:(r['partner_sse'],r['rna_sse'],r['key']))}
    dump(out/'SELECTION_BEFORE_OUTER.json',dict(cell=cell,baseline_sse=base_sse,baseline_partner_sse=base_partner,selected=groups,candidates=rows,outer_truth_read=False))
    base,candidates,w,audit=one_fit(d,np.arange(m),d['query'],plan)
    predictions={'M0':base};states={}
    for label,row in groups.items():
        if row['key']=='M0':predictions[label]=base.copy();continue
        a,c=candidates[row['key']];predictions[label]=base+a+row['beta']*c
        states[label+'__endpoint']=a;states[label+'__partner']=c
    # Preserve all finite candidates at beta=1 for diagnostic comparisons;
    # family selection remains exclusively in SELECTION_BEFORE_OUTER.json.
    np.savez_compressed(out/'PREDICTIONS.npz',edges=d['query'],ids=d['query_ids'],**predictions,**states)
    dump(out/'final_audit.json',dict(reliability=w,models=audit,seconds=time.monotonic()-begin))
    files=[p for p in out.iterdir() if p.is_file()]
    dump(out/'FREEZE.json',dict(cell=cell,code_sha256=sha(__file__),plan_sha256=sha(ROOT/'PLAN.json'),input_sha256=sha(ROOT/'inputs'/f'{cell}.json'),files={p.name:sha(p) for p in files},outer_truth_read=False))
    print(json.dumps({'cell':cell,'phase':'frozen','selected':chosen,'seconds':time.monotonic()-begin}),flush=True)


def smoke():
    rng=np.random.default_rng(64);n=9;edges=np.array([(i,j) for i in range(n) for j in range(i+1,n)])
    P=Projection(edges,n);r=rng.normal(size=(len(edges),2))
    assert np.max(abs(P.X.T@P.apply(r)))<1e-12
    direct=r-P.X@np.linalg.lstsq(P.X,r,rcond=None)[0]
    assert np.max(abs(P.apply(r)-direct))<1e-12
    f=rng.normal(size=(n,3));F=symfeatures(f,edges);K=f@f.T
    a,b=edges.T;expected=.5*(K[np.ix_(a,a)]*K[np.ix_(b,b)]+K[np.ix_(a,b)]*K[np.ix_(b,a)])
    assert np.max(abs(F@F.T-expected))<1e-11
    # Equivalent feature implementation of the archived train/query projection.
    tr=edges[:25];q=edges[25:];pt=Projection(tr,n)
    ft,fq=F[:25],F[25:];pr=pt.apply(ft);pq=fq-pt.extend(ft,q)
    A=np.linalg.pinv(pt.X);pm=np.eye(len(tr))-pt.X@A;xq=Projection(q,n).X
    raw=F@F.T;qt=(raw[25:,:25]-xq@A@raw[:25,:25])@pm
    assert np.max(abs(pq@pr.T-qt))<1e-11
    st=np.array([9.,-3.,1.,9.,-3.,1.]);assert abs(solve_beta(st,1.)['beta']-3)<1e-12
    u=rng.normal(size=n);u/=np.linalg.norm(u);truth=(np.outer(u,u))[edges[:,0],edges[:,1]]
    fit=fit_projected_rank1(n,edges,truth[:,None].repeat(2,axis=1),0.,0.,P)
    err=P.apply(fit['G'][edges[:,0],edges[:,1]]-truth)
    assert np.mean(err*err)<1e-12
    report={'projection_lstsq_equivalence':True,'symmetric_tensor_native_kernel_equivalence':True,'next2_train_query_projection_equivalence':True,'next2_scalar_solver_reused':True,'projected_rank1_positive_control_mse':float(np.mean(err*err))}
    dump(ROOT/'diagnostics/SMOKE.json',report);print(report)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--cell');ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    if args.smoke:smoke()
    else:run(args.cell)
