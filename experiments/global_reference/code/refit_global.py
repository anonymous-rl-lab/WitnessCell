"""Matched global-product reference; outcomes remain in a separate file."""
import os
for _k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[_k]='1'
from pathlib import Path
import sys,json,hashlib,time,argparse
import numpy as np
from scipy import linalg
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'vendor'))
import frozen_models as fm
import spectral_models as sm

def clean(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [clean(x) for x in v]
    return v
def dump(p,v):p.write_text(json.dumps(clean(v),indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class Project:
    def __init__(self,ed,n):
        self.X=np.zeros((len(ed),n));self.X[np.arange(len(ed))[:,None],ed]=1
        self.inv=np.linalg.pinv(self.X.T@self.X,hermitian=True)
    def p(self,a):return a-self.X@(self.inv@(self.X.T@a))

def fit_base(singles,ed,y,product=True):
    scale=np.std(singles,axis=0);scale[scale==0]=1;models=[]
    for j in range(y.shape[1]):
        a=singles[:,j]/scale[j];de=fm._Design(a,ed,y[:,j]/scale[j],10.,0,0.)
        X=np.column_stack((de.X,a[ed[:,0]]*a[ed[:,1]])) if product else de.X
        pen=np.r_[de.penalty,0.] if product else de.penalty
        ix=np.flatnonzero(pen);A=np.vstack((X/np.sqrt(len(y)),np.diag(np.sqrt(pen))[ix]));b=np.r_[(y[:,j]/scale[j]-de.offset)/np.sqrt(len(y)),np.zeros(len(ix))]
        coef=linalg.lstsq(A,b,lapack_driver='gelsy',check_finite=False)[0]
        native=de.make_model(coef[:-1] if product else coef).to_dict()
        models.append(dict(native=native,gamma=float(coef[-1]) if product else 0.,scale=float(scale[j])))
    return models
def predict(models,ed):
    cols=[]
    for st in models:
        m=fm.FittedModel.from_dict(st['native']);cols.append((m.predict(ed)+st['gamma']*m.a[ed[:,0]]*m.a[ed[:,1]])*st['scale'])
    return np.column_stack(cols)
def fit(label):
    O=R/'fitted'/label;O.mkdir(parents=True,exist_ok=True)
    if (O/'LOCK.json').exists():return
    start=time.monotonic();d=dict(np.load(R/'inputs'/label/'training.npz'));ed,y,sing,folds,qu=[d[k] for k in ('edges','y','singles','folds','query_edges')];n=len(sing)
    orig=dict(np.load(R/'inputs'/label/'original.npz'))
    native=fit_base(sing,ed,y,False);gap=float(np.max(abs(predict(native,qu)-orig['M0'])));assert gap<1e-10,(label,gap)
    oo=np.zeros_like(y);grid={};probes=[]
    for fam in ['rank1','spectral']:
        for nu in [.1,1.]:
            for frac in [.3,.1,.03,.01]:grid[f'{fam}_t{frac:g}_n{nu:g}']=dict(family=fam,nu=nu,fraction=frac,D=np.zeros_like(y),eligible=True)
    for k in range(3):
        tr=folds!=k;va=~tr;st=fit_base(sing,ed[tr],y[tr]);oo[va]=predict(st,ed[va]);res=y[tr]-predict(st,ed[tr]);tz=sm.tau_zero(n,ed[tr],res)
        dump(O/f'inner{k}_baseline.json',st)
        for nu in [.1,1.]:
            warm=None
            for frac in [.3,.1,.03,.01]:
                ms={}
                for fam in ['rank1','spectral']:
                    fun=sm.fit_rank1 if fam=='rank1' else sm.fit_spectral
                    m=fun(n,ed[tr],res,tz*frac,nu,max_iter=3000,tol=1e-6,initial=warm if fam=='spectral' else None)
                    if fam=='spectral':warm=m['G']
                    key=f'{fam}_t{frac:g}_n{nu:g}';row=grid[key];row['eligible'] &= bool(m['converged']) and not m.get('timed_out',False)
                    supported=np.bincount(ed[tr].ravel(),minlength=n)>0;dd=m['G'][ed[va,0],ed[va,1]].copy();dd[~supported[ed[va]].all(axis=1)]=0;row['D'][va]=dd[:,None];ms[fam]=m
                gap_obj=ms['spectral']['objective']-ms['rank1']['objective'];tol=1e-7+1e-5*max(abs(ms['rank1']['objective']),abs(ms['spectral']['objective']));assert gap_obj<=tol
                probes.append(dict(fold=k,nu=nu,fraction=frac,objective_gap=gap_obj,tolerance=tol))
        print(label,'inner',k,'seconds',round(time.monotonic()-start,1),flush=True)
    bsse=float(np.sum((y-oo)**2));rows=[]
    for key,row in grid.items():rows.append({k:v for k,v in row.items() if k!='D'}|dict(key=key,sse=float(np.sum((y-oo-row['D'])**2))))
    st=fit_base(sing,ed,y);bq=predict(st,qu);r=y-predict(st,ed);tz=sm.tau_zero(n,ed,r);pay={'B':bq,'B_oof':oo,'query_edges':qu};sel={}
    dump(O/'baseline.json',st)
    for fam in ['rank1','spectral']:
        candidates=[v for v in rows if v['family']==fam and v['eligible']];assert candidates
        best=min(candidates,key=lambda v:(v['sse'],-v['fraction'],-v['nu']));use=best['sse']<bsse-1e-12*max(1.,bsse)
        fun=sm.fit_rank1 if fam=='rank1' else sm.fit_spectral;m=fun(n,ed,r,tz*best['fraction'],best['nu'],max_iter=3000,tol=1e-6);assert m['converged']
        G=m['G'] if use else np.zeros((n,n));D=np.repeat(G[qu[:,0],qu[:,1]][:,None],y.shape[1],axis=1);Do=grid[best['key']]['D'] if use else np.zeros_like(y)
        aa=bb=0.
        for k in range(3):
            va=folds==k;p=Project(ed[va],n);er=p.p(y[va]-oo[va]);dr=p.p(Do[va]);aa+=np.sum(er*dr);bb+=np.sum(dr*dr)
        alpha=float(np.clip(aa/bb,0,1)) if bb>1e-14 else 0.
        pay.update({fam+'_D':D,fam+'_G':G,fam+'_D_oof':Do,fam+'_controlled':bq+D+(alpha-1)*Project(qu,n).p(D)})
        sel[fam]=dict(selected=best['key'] if use else 'G0',best_regularized=best['key'],alpha=alpha,inner_sse=best['sse'],baseline_inner_sse=bsse,final_converged=bool(m['converged']))
        dump(O/(fam+'_fit.json'),{k:v for k,v in m.items() if k!='G'})
    np.savez_compressed(O/'PREDICTIONS.npz',**pay);dump(O/'SELECTION.json',sel);dump(O/'INNER_GRID.json',dict(candidates=rows,probes=probes))
    dump(O/'LOCK.json',dict(configuration=label,seconds=time.monotonic()-start,query_outcomes_opened=False,native_baseline_replay_max_abs=gap,plan_sha256=sha(R/'ANALYSIS_PLAN.json'),training_sha256=sha(R/'inputs'/label/'training.npz'),predictions_sha256=sha(O/'PREDICTIONS.npz'),selection_sha256=sha(O/'SELECTION.json'),code_sha256=sha(Path(__file__))))
    print(label,'LOCKED',sel,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--label');a=ap.parse_args()
    labels=[a.label] if a.label else [x['label'] for x in json.loads((R/'CONFIGURATIONS.json').read_text())]
    for label in labels:fit(label)
