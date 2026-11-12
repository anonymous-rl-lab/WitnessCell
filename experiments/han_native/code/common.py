import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,hashlib,datetime
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    return x
def dump(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(clean(x),indent=2,allow_nan=False)+'\n');tmp.replace(p)
def savez(p,**x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.tmp.npz');np.savez_compressed(tmp,**x);tmp.replace(p)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def loadnpz(p):
    with np.load(p,allow_pickle=False) as f:return {k:f[k] for k in f.files}
def innerfold(pairs):
    return np.array([int.from_bytes(hashlib.sha256(('wc_spectral_smoke_v1|'+'|'.join(sorted(p))).encode()).digest()[:8],'big')%3 for p in pairs])
class Project:
    def __init__(self,edges,n):
        self.X=np.zeros((len(edges),n));self.X[np.arange(len(edges)),edges[:,0]]=1;self.X[np.arange(len(edges)),edges[:,1]]=1
        v,u=np.linalg.eigh(self.X.T@self.X);keep=v>max(self.X.shape)*np.finfo(float).eps*max(v.max(),1)
        self.inv=(u[:,keep]/v[keep])@u[:,keep].T
    def coef(self,a):return self.inv@(self.X.T@np.asarray(a).reshape(len(a),-1))
    def p(self,a):return (a.reshape(len(a),-1)-self.X@self.coef(a)).reshape(a.shape)
