#!/usr/bin/env python3
"""Native conditional-RNA predictor without an unperturbed condition.

The mean of permitted single profiles is an origin, not a fabricated control.
All fitted operations receive explicit training indices. Batch-independent
partner corrections use a fixed identity-pair universe.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
import hashlib, itertools, json
from pathlib import Path
import numpy as np

ALPHAS=(0.,.02,.05,.1,.2,.5)
NOISES=(.01,.03,.1,.3,1.)
LAMBDAS=(.1,1.,10.)

def projection(x):
    p=np.eye(len(x))-x@np.linalg.pinv(x)
    return (p+p.T)/2

def sat(z,a):return z/(1+a*np.abs(z))

def incidence(pairs,nodes):
    x=np.zeros((len(pairs),len(nodes)))
    for i,p in enumerate(pairs):
        for g in p.split('+'):x[i,nodes.index(g)]=1.
    return x

def label_folds(allowed,pairs,k=3):
    """Deterministic geometry-only cycle allocation, no expression values."""
    ids=list(map(int,allowed));nodes=sorted({g for i in ids for g in pairs[i].split('+')})
    lut={pairs[i]:i for i in ids};rects=[]
    for ns in itertools.combinations(nodes,4):
        a,b,c,d=ns
        for ee in (((a,b),(a,d),(c,b),(c,d)),((a,c),(a,d),(b,c),(b,d)),((a,b),(a,c),(d,b),(d,c))):
            ps=['+'.join(sorted(e)) for e in ee]
            if all(p in lut for p in ps):rects.append(tuple(sorted(lut[p] for p in ps)))
    digest=lambda group:hashlib.sha256(('originfree-v1|'+'|'.join(sorted(pairs[i] for i in group))).encode()).hexdigest()
    rects=sorted(set(rects),key=digest);used=set();groups=[]
    for r in rects:
        if not used.intersection(r):groups.append(list(r));used.update(r)
    groups += [[i] for i in sorted(set(ids)-used,key=lambda i:digest([i]))]
    folds=[[] for _ in range(min(k,len(ids)))]
    for g in groups:
        j=min(range(len(folds)),key=lambda j:(len(folds[j]),j));folds[j].extend(g)
    return [np.array(sorted(f),int) for f in folds if f]

class Data:
    def __init__(self,pairs,nodes,singles,truth,halves,genes=None):
        self.pairs=list(pairs);self.nodes=list(nodes);self.genes=genes
        self.single=np.asarray(singles,float);self.origin=self.single.mean(axis=0)
        self.x=incidence(pairs,nodes);self.e=self.single-self.origin
        self.y=np.asarray(truth,float)-self.origin
        self.halves=np.asarray(halves,float)-self.origin
        self.additive=.5*self.x@self.e
        self.gram=self.e@self.e.T/self.e.shape[1]
        self.edges=np.array([[nodes.index(g) for g in p.split('+')] for p in pairs],int)
        self.universe=['+'.join(p) for p in itertools.combinations(nodes,2)]
        self.xref=incidence(self.universe,nodes)
        self.eref=np.array([[nodes.index(g) for g in p.split('+')] for p in self.universe],int)

def backbone(d,q,alpha):
    return np.zeros_like(d.y[q]) if alpha is None else sat(d.additive[q],alpha)

def baseline(d,tr,q,pars):
    tr=np.asarray(tr,int);q=np.asarray(q,int)
    xt=d.x[tr];xm=xt.mean(axis=0);xc=xt-xm
    r=d.y[tr]-backbone(d,tr,pars['alpha']);u=r.mean(axis=0)
    beta=np.linalg.solve(xc.T@xc+2*pars['noise']*np.eye(xc.shape[1]),xc.T@(r-u))
    return backbone(d,q,pars['alpha'])+u+pars['gamma']*((d.x[q]-xm)@beta)

def saturation_predict(d,tr,q,alpha):
    return backbone(d,q,alpha)+(d.y[tr]-backbone(d,tr,alpha)).mean(axis=0)

def tune_baseline(d,tr,folds):
    tr=np.asarray(tr,int);allowed=set(tr)
    assert set(np.concatenate(folds))==allowed and sum(map(len,folds))==len(tr)
    alphas=(*ALPHAS,None)
    dot=np.zeros((len(alphas),len(NOISES)));raw2=np.zeros_like(dot);truth2=np.zeros(len(alphas))
    for va in folds:
        fit=np.array(sorted(allowed-set(va)),int);xt=d.x[fit];xm=xt.mean(axis=0);xc=xt-xm;xv=d.x[va]-xm
        xx=xv.T@xv;invs=[np.linalg.inv(xc.T@xc+2*eta*np.eye(xc.shape[1])) for eta in NOISES]
        for ai,alpha in enumerate(alphas):
            r=d.y[fit]-backbone(d,fit,alpha);u=r.mean(axis=0);r=r-u
            target=d.y[va]-backbone(d,va,alpha)-u
            q=xc.T@r;t=xv.T@target;cross=t@q.T;qq=q@q.T
            truth2[ai]+=float(np.sum(target*target))
            for ni,a in enumerate(invs):
                dot[ai,ni]+=float(np.sum(a*cross));raw2[ai,ni]+=float(np.sum((a@qq@a.T)*xx))
    gamma=np.divide(dot,raw2,out=np.zeros_like(dot),where=raw2>1e-12).clip(0,1)
    loss=truth2[:,None]-2*gamma*dot+gamma*gamma*raw2
    def pars(ai,ni):return {'alpha':alphas[ai],'noise':NOISES[ni],'gamma':float(gamma[ai,ni]),'oof_sse':float(loss[ai,ni]),'family':'free_endpoint' if alphas[ai] is None else 'single_mean_saturation_endpoint'}
    best=min(np.ndindex(loss.shape),key=lambda ij:(loss[ij],-NOISES[ij[1]],ij[0]))
    free=min(range(len(NOISES)),key=lambda ni:(loss[-1,ni],-NOISES[ni]))
    anchor=min(itertools.product(range(len(ALPHAS)),range(len(NOISES))),key=lambda ij:(loss[ij],-NOISES[ij[1]],ij[0]))
    sa=int(np.argmin(truth2[:len(ALPHAS)]))
    return {'selected':pars(*best),'free_endpoint':pars(len(alphas)-1,free),'anchored_M0':pars(*anchor),
            'saturation_alpha':ALPHAS[sa],'saturation_oof_sse':float(truth2[sa]),
            'losses':loss.tolist(),'gammas':gamma.tolist(),'folds':[v.tolist() for v in folds]}

def cross(a,b,gram):
    x,y=a.T;u,v=b.T
    return .5*(gram[np.ix_(x,u)]*gram[np.ix_(y,v)]+gram[np.ix_(x,v)]*gram[np.ix_(y,u)])

class Operator:
    def __init__(self,d,tr):
        self.d=d;self.tr=np.asarray(tr,int);self.xt=d.x[tr];self.et=d.edges[tr]
        self.plus=np.linalg.pinv(self.xt);self.p=projection(self.xt)
        self.ktt=cross(self.et,self.et,d.gram)
        k=self.p@self.ktt@self.p;k=(k+k.T)/2
        self.scale=float(np.trace(k)/len(tr)) if np.trace(k)>1e-12*max(float(np.trace(self.ktt)),1e-300) else 0.
        self.eig,self.v=np.linalg.eigh(k/self.scale if self.scale else k*0)
        self.refbeta={}
        for lam in LAMBDAS:
            self.refbeta[lam]=np.linalg.pinv(d.xref)@self.raw(d.xref,d.eref,lam)
    def raw(self,xq,eq,lam):
        if not self.scale:return np.zeros((len(eq),len(self.tr)))
        k=(cross(eq,self.et,self.d.gram)-xq@self.plus@self.ktt)@self.p/self.scale
        return ((k@self.v)/(self.eig+lam))@self.v.T
    def predict_operator(self,q,lam):
        q=np.asarray(q,int);return self.raw(self.d.x[q],self.d.edges[q],lam)-self.d.x[q]@self.refbeta[lam]

def weights(d,tr):
    p=projection(d.x[tr]);a,b=[p@v[tr] for v in d.halves]
    den=np.sum(a*a+b*b,axis=0);num=2*np.sum(a*b,axis=0)
    return np.divide(num,den,out=np.zeros_like(num),where=den>0).clip(0,1)

def fit(d,tr,q,inners):
    tr=np.asarray(tr,int);q=np.asarray(q,int)
    pars=tune_baseline(d,tr,inners);base=baseline(d,tr,q,pars['selected'])
    scores={k:{lam:0. for lam in LAMBDAS} for k in ('R','M1')};dimensions=0;info=[]
    for va in inners:
        train=np.array(sorted(set(tr)-set(va)),int);aux=label_folds(train,d.pairs)
        ip=tune_baseline(d,train,aux);b=baseline(d,train,va,ip['selected'])
        resid=d.y[train]-baseline(d,train,train,ip['selected']);w=weights(d,train);op=Operator(d,train)
        p=projection(d.x[va]);df=len(va)-np.linalg.matrix_rank(d.x[va]);dimensions+=df
        for lam in LAMBDAS:
            g=op.predict_operator(va,lam)@resid
            for kind in scores:
                err=b+g*(w if kind=='R' else 1.)-d.y[va]
                scores[kind][lam]+=float(np.sum((p@err if df>0 else err)**2))
        info.append({'train':train.tolist(),'validation':va.tolist(),'projection_df':int(df),'fallback_full':bool(df==0),'baseline':ip})
    choices={k:min(LAMBDAS,key=lambda lam:(scores[k][lam],-lam)) for k in scores}
    op=Operator(d,tr);w=weights(d,tr);resid=d.y[tr]-baseline(d,tr,tr,pars['selected'])
    pred={'M0':base,'single_average':d.additive[q],'saturation':saturation_predict(d,tr,q,pars['saturation_alpha']),
          'free_endpoint':baseline(d,tr,q,pars['free_endpoint']),'anchored_M0':baseline(d,tr,q,pars['anchored_M0'])}
    state={'train':tr,'query':q,'weights':w,'origin':d.origin,'M0':base}
    for lam in LAMBDAS:
        o=op.predict_operator(q,lam);g=o@resid;state['operator_'+str(lam)]=o
        pred['M1_lambda'+str(lam)]=base+g;pred['R_lambda'+str(lam)]=base+g*w
    for kind,lam in choices.items():pred[kind]=pred[kind+'_lambda'+str(lam)]
    audits=[]
    for lam in LAMBDAS:
        one=np.concatenate([op.predict_operator([j],lam) for j in q])
        full=op.predict_operator(q,lam);reverse=op.predict_operator(q[::-1],lam)[::-1]
        audits.append(max(float(np.max(abs(one-full))),float(np.max(abs(reverse-full)))))
    assert max(audits)<1e-10
    return pred,state,{'baseline':pars,'choices':choices,'scores':scores,'inner':info,'inner_projection_df_total':int(dimensions),
                       'query_invariance_max_abs':max(audits),'weight_mean':float(w.mean()),'weight_median':float(np.median(w))}

def array_hash(a):return hashlib.sha256(np.ascontiguousarray(a,dtype=np.float64).tobytes()).hexdigest()

def evaluate(d,pred,foldids):
    rows=[]
    for kind,y in pred.items():
        full=float(np.sum((y-d.y)**2));psse=0.;fold=[]
        for f in sorted(set(foldids)):
            q=np.flatnonzero(foldids==f);err=y[q]-d.y[q];p=projection(d.x[q]);v=float(np.sum((p@err)**2));psse+=v
            fold.append({'fold':int(f),'full_sse':float(np.sum(err**2)),'P_sse':v})
        rows.append({'method':kind,'full_sse':full,'P_sse':psse,'folds':fold})
    m0=next(r for r in rows if r['method']=='M0')
    for r in rows:
        for k in ('full_sse','P_sse'):r[k+'_gain_pct']=100*(1-r[k]/m0[k])
    return rows
