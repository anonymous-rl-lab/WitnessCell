import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import hashlib,json,datetime
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def load(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def split(pairs,ids,folds,seed):
    order=sorted(ids,key=lambda i:hashlib.sha256((str(seed)+'|'+'|'.join(sorted(pairs[i]))).encode()).hexdigest())
    return [np.asarray(a,dtype=int) for a in np.array_split(order,folds)]
class Projection:
    def __init__(self,pairs):
        self.genes=sorted(set(pairs.ravel()));ix={g:i for i,g in enumerate(self.genes)}
        self.edges=np.array([[ix[a],ix[b]] for a,b in pairs]);n=len(self.genes)
        self.X=np.zeros((len(pairs),n))
        self.X[np.arange(len(pairs)),self.edges[:,0]]=1
        self.X[np.arange(len(pairs)),self.edges[:,1]]=1
        gram=self.X.T@self.X
        w,v=np.linalg.eigh(gram);keep=w>max(gram.shape)*np.finfo(float).eps*max(float(w.max()),1.)
        self.rank=int(keep.sum());self.dim=len(pairs)-self.rank
        self.inv=(v[:,keep]/w[keep])@v[:,keep].T
    def H(self,a):
        shape=a.shape;b=a.reshape(len(a),-1)
        return (self.X@(self.inv@(self.X.T@b))).reshape(shape)
    def P(self,a):return a-self.H(a)
