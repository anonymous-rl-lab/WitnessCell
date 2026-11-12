from common import *
import sys
sys.path.insert(0,str(Path(__file__).parent/'vendor'))
from frozen_models import fit_model,FittedModel

LAM=np.array([0,.05,.1,.2,.5,1,2,5])
def design(single,edges):
    b=single[edges[:,0]]+single[edges[:,1]];scale=np.sqrt(np.mean(single**2,axis=0));scale[scale==0]=1
    return b,b[:,:,None]/(1+abs(b[:,:,None])/scale[None,:,None]*LAM[None,None,:])-b[:,:,None]
def satfit(P,f,e):
    pe,pf=P.p(e),P.p(f);c=[]
    for j in range(e.shape[1]):
        g=float(pf[:,j]@pf[:,j]);h=float(pf[:,j]@pe[:,j]);c.append(h/g if g>1e-30 else 0.)
    c=np.array(c);en=P.coef(e-f*c);return c,en
def fit(kind,singles,edges,y,pairs,seed,path):
    path=Path(path)
    if path.exists():return json.loads(path.read_text())
    if kind=='native':
        scale=np.std(singles,axis=0);scale[scale==0]=1
        models=[fit_model(singles[:,j]/scale[j],edges,y[:,j]/scale[j],rank=0,alpha=10.,deadline=None).to_dict() for j in range(y.shape[1])]
        state=dict(kind=kind,scale=scale,models=models)
    else:
        B,F=design(singles,edges);n=len(singles)
        order=sorted(range(len(edges)),key=lambda i:hashlib.sha256(('v37-inner|'+str(seed)+'|'+'|'.join(sorted(pairs[i]))).encode()).digest())
        scores=np.zeros(8)
        for va in np.array_split(order,3):
            tr=np.setdiff1d(np.arange(len(edges)),va);pt,pv=Project(edges[tr],n),Project(edges[va],n)
            for i in range(8):
                c,en=satfit(pt,F[tr,:,i],y[tr]-B[tr]);pred=B[va]+F[va,:,i]*c+pv.X@en;scores[i]+=np.sum((y[va]-pred)**2)
        ix=int(np.argmin(scores));c,en=satfit(Project(edges,n),F[:,:,ix],y-B)
        state=dict(kind=kind,index=ix,lambda_value=float(LAM[ix]),coef=c,endpoint=en,inner_sse=scores)
    dump(path,state);return clean(state)
def predict(state,singles,edges):
    if state['kind']=='native':return np.column_stack([FittedModel.from_dict(m).predict(edges)*state['scale'][j] for j,m in enumerate(state['models'])])
    B,F=design(singles,edges);en=np.array(state['endpoint']);return B+F[:,:,state['index']]*np.array(state['coef'])+en[edges[:,0]]+en[edges[:,1]]
