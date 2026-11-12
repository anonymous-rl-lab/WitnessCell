#!/usr/bin/env python3
"""Frozen-universe postfit program reference, bounded retrospective probe."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse,hashlib,importlib.util,itertools,json,sys,time
from pathlib import Path
import numpy as np

METHODS=('raw','query','global','complement')
LAMBDAS=(.1,1.,10.)

def mod(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def pair_cross(a,b,gram):
    x,y=a.T;u,v=b.T
    return .5*(gram[np.ix_(x,u)]*gram[np.ix_(y,v)]+gram[np.ix_(x,v)]*gram[np.ix_(y,u)])

class FixedFit:
    """Only fitted labels, measured single information and fixed node universe."""
    def __init__(self,task,fit,pi,gm,b):
        self.gm=gm;self.b=b;self.nodes=tuple(task.nodes);self.node_index={v:i for i,v in enumerate(self.nodes)}
        self.fit_pairs=tuple(task.pairs[i] for i in fit)
        self.universe=tuple('+'.join(p) for p in itertools.combinations(self.nodes,2))
        self.complement=tuple(p for p in self.universe if p not in set(self.fit_pairs))
        self.allowed=('control',*self.nodes,*self.fit_pairs)
        matrix=np.stack([task.values[v] for v in self.allowed]);var=np.var(matrix,axis=0,ddof=1)
        self.floor=float(np.quantile(var[var>0],.01));self.sigma=np.sqrt(np.maximum(var,self.floor));self.pi=pi.copy()
        self.e=np.stack([task.values[g]-task.values['control'] for g in self.nodes]);self.gram=self.e@self.e.T/self.e.shape[1]
        self.xt,self.et=self.design(self.fit_pairs);self.plus=np.linalg.pinv(self.xt);self.pt=b.projection(self.xt)
        self.ktt=pair_cross(self.et,self.et,self.gram)
        target=self.h(task.d.y[fit]);add=self.h(self.xt@self.e)
        self.alpha=min(gm.ALPHAS,key=lambda a:(float(np.sum((gm.sat(add,a)-target)**2)),a))
        self.resid=target-gm.sat(add,self.alpha)
        kt=self.pt@self.ktt@self.pt;kt=(kt+kt.T)/2
        trace=float(np.trace(kt));self.scale=trace/len(fit) if trace>1e-12*max(abs(float(np.trace(self.ktt))),1e-300) else 0.
        kt=kt/self.scale if self.scale else np.zeros_like(kt)
        self.eig,self.V=np.linalg.eigh(kt)
        self.c={lam:((self.V/(self.eig+lam))@self.V.T)@self.resid for lam in LAMBDAS}
        xf,_=self.design(self.universe);xc,_=self.design(self.complement)
        self.beta_global={lam:np.linalg.pinv(xf)@self.raw(self.universe,lam) for lam in LAMBDAS}
        self.beta_complement={lam:np.linalg.pinv(xc)@self.raw(self.complement,lam) for lam in LAMBDAS}
    def h(self,y):return np.mean(np.asarray(y)[...,self.pi]/self.sigma[self.pi],axis=-1,keepdims=True)
    def design(self,pairs):
        edges=np.asarray([[self.node_index[g] for g in p.split('+')] for p in pairs],int)
        if not len(edges):return np.zeros((0,len(self.nodes))),np.empty((0,2),int)
        if np.any(edges[:,0]==edges[:,1]):raise ValueError('Only distinct endpoint pairs supported')
        x=np.zeros((len(edges),len(self.nodes)));x[np.arange(len(edges)),edges[:,0]]=1.;x[np.arange(len(edges)),edges[:,1]]=1.
        return x,edges
    def base(self,pairs):
        x,_=self.design(pairs);return self.gm.sat(self.h(x@self.e),self.alpha)
    def additive(self,pairs):
        x,_=self.design(pairs);return self.h(x@self.e)
    def raw(self,pairs,lam):
        x,edges=self.design(pairs)
        if not self.scale:return np.zeros((len(pairs),1))
        kq=(pair_cross(edges,self.et,self.gram)-x@self.plus@self.ktt)@self.pt/self.scale
        return kq@self.c[lam]
    def correction(self,pairs,lam,method):
        x,_=self.design(pairs);g=self.raw(pairs,lam)
        if method=='raw':return g
        if method=='query':return self.b.projection(x)@g
        if method=='global':return g-x@self.beta_global[lam]
        if method=='complement':return g-x@self.beta_complement[lam]
        raise ValueError(method)
    def predict(self,pairs,lam,method):return self.base(pairs)+self.correction(pairs,lam,method)
    def save(self,path):
        np.savez_compressed(path,nodes=np.array(self.nodes),fit_pairs=np.array(self.fit_pairs),universe=np.array(self.universe),
            complement=np.array(self.complement),sigma=self.sigma,pi=self.pi,alpha=self.alpha,e=self.e,gram=self.gram,xt=self.xt,
            et=self.et,plus=self.plus,pt=self.pt,ktt=self.ktt,scale=self.scale,resid=self.resid,eig=self.eig,V=self.V,
            **{f'c_{lam}':self.c[lam] for lam in LAMBDAS},
            **{f'beta_global_{lam}':self.beta_global[lam] for lam in LAMBDAS},
            **{f'beta_complement_{lam}':self.beta_complement[lam] for lam in LAMBDAS})

def fit_outer(task,train,fold,pi,gm,b,geo):
    scores={m:{lam:0. for lam in LAMBDAS} for m in METHODS};records=[];zero=0.
    for j,held in enumerate(geo.graph_folds(train,task.pairs,fold,task.d.x,'cycle3')):
        fit=np.array(sorted(set(train)-set(held)));state=FixedFit(task,fit,pi,gm,b);qp=[task.pairs[i] for i in held]
        truth=state.h(task.d.y[held]);err=state.base(qp)-truth;zero+=float(np.sum(err**2))
        detail={m:{} for m in METHODS}
        for method in METHODS:
            for lam in LAMBDAS:
                ss=float(np.sum((state.predict(qp,lam,method)-truth)**2));scores[method][lam]+=ss;detail[method][lam]=ss
        records.append({'inner':j,'fit':fit.tolist(),'validation':held.tolist(),'allowed_normalization':list(state.allowed),
            'alpha':state.alpha,'variance_floor':state.floor,'scores':detail,'complement_reference_size':len(state.complement)})
    choices={m:min(LAMBDAS,key=lambda lam:(scores[m][lam],-lam)) for m in METHODS}
    state=FixedFit(task,train,pi,gm,b)
    return state,choices,{'scores':scores,'choices':choices,'zero_full_sse':zero,'inner':records}

def invariance(state,qp,lam):
    rng=np.random.default_rng(7719);perm=rng.permutation(len(qp));rows=[]
    for method in METHODS:
        whole=state.predict(qp,lam,method);single=np.concatenate([state.predict([p],lam,method) for p in qp])
        chunks=np.concatenate([state.predict(qp[j:j+4],lam,method) for j in range(0,len(qp),4)])
        reverse=state.predict(['+'.join(p.split('+')[::-1]) for p in qp],lam,method)
        order=state.predict([qp[i] for i in perm],lam,method)[np.argsort(perm)]
        duplicates=state.predict(qp+[qp[0]],lam,method)
        diffs={'singleton':float(np.max(abs(single-whole))),'chunks4':float(np.max(abs(chunks-whole))),
            'order':float(np.max(abs(order-whole))),'endpoint_swap':float(np.max(abs(reverse-whole))),
            'duplicate_request':float(np.max(abs(duplicates[:-1]-whole)))}
        rows.append({'method':method,'lambda':lam,**diffs,'max_difference':max(diffs.values()),
            'stable_fixed_reference_expected':method!='query','PASS':max(diffs.values())<=1e-12 if method!='query' else None})
    return rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--asset-root',type=Path);ap.add_argument('--compact-root',type=Path);ap.add_argument('--bridge-code',type=Path,required=True)
    ap.add_argument('--geometry-code',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    out=a.out;out.mkdir(parents=True,exist_ok=True)
    if (out/'STARTED.json').exists():raise FileExistsError(out)
    if a.compact_root is None and a.asset_root is None:ap.error('Provide --asset-root or --compact-root')
    b=mod('reference_bridge',a.bridge_code);geo=mod('reference_geo',a.geometry_code);r=a.asset_root
    if r is None:ap.error('--asset-root supplies the original immutable model source; --compact-root replaces scientific data only')
    gm=b.module('global_pair_model',r/'experiment_assets/WitnessCell_Partner_Component_Run_20261018/model/global_pair_model.py')
    if a.compact_root is not None:
        loader=mod('reference_compact_loader',Path(__file__).resolve().parent.parent/'compact_loader.py');task=loader.load_task(a.compact_root,gm,b)
        program_path=a.compact_root/'PREDEFINED_G2M_PROGRAM.json'
    else:
        old=b.module('reference_old',r/'source_v19/WitnessCell_Paper_v19_Full/source_data/RNA_controls/run_controls.py')
        new=b.module('reference_new',r/'work_B/run_frozen_rna_confirmation.py');task=b.load_task('Dixit',r,gm,new,old)
        program_path=r/'work_B/PREDEFINED_G2M_PROGRAM.json'
    program=json.loads(program_path.read_text());pi=np.array([int(np.flatnonzero(task.symbols==s)[0]) for s in program['genes'] if np.count_nonzero(task.symbols==s)==1]);assert len(pi)==52
    universe=tuple('+'.join(p) for p in itertools.combinations(task.nodes,2));assert set(universe)==set(task.pairs);assert len(universe)==45
    started={'utc':b.stamp(),'plan_sha256':b.sha(Path(__file__).with_name('PLAN.txt')),'script_sha256':b.sha(__file__),
        'bridge_sha256':b.sha(a.bridge_code),'geometry_sha256':b.sha(a.geometry_code),'model_source_sha256':b.sha(gm.__file__),
        'program_json_sha256':b.sha(program_path),'compact_loader_sha256':b.sha(loader.__file__) if a.compact_root is not None else None,
        'status':'retrospective development',
        'program_symbols':task.symbols[pi].tolist(),'identity_only_complete_universe':list(universe),'source_hashes':{str(p):b.sha(p) for p in task.source_paths}}
    b.dump(out/'STARTED.json',started);tick=time.monotonic();pred={};allstates=[];foldids=np.full(len(task.pairs),-1)
    lookup={p:i for i,p in enumerate(task.pairs)};prior={0:.1,1:10.,2:10.};audits=[];inv=[]
    def fill(name,q,v):
        if name not in pred:pred[name]=np.full(len(task.pairs),np.nan)
        pred[name][q]=np.asarray(v).ravel()
    for fold in task.folds:
        f=fold['fold'];tr=np.array([lookup[p] for p in fold['train_pairs']]);q=np.array([lookup[p] for p in fold['test_pairs']]);qp=[task.pairs[i] for i in q]
        state,choices,info=fit_outer(task,tr,f,pi,gm,b,geo);foldids[q]=f
        fill('scalar_saturation',q,state.base(qp));fill('measured_additive',q,state.additive(qp))
        corrections={}
        for method in METHODS:
            for regime,lam in [('matched',prior[f]),('selected',choices[method])]:fill(f'{method}_{regime}',q,state.predict(qp,lam,method))
            for lam in LAMBDAS:
                fill(f'{method}_fixed_{lam}',q,state.predict(qp,lam,method));corrections[f'{method}_{lam}']=state.correction(qp,lam,method)
        state.save(out/f'FITTED_STATE_{f}.npz')
        np.savez_compressed(out/f'CORRECTIONS_{f}.npz',q=q,train=tr,raw_universe=np.stack([state.raw(universe,lam) for lam in LAMBDAS]),**corrections)
        b.dump(out/f'FOLD_{f}.json',{'fold':f,'train':tr.tolist(),'query':q.tolist(),'prior_lambda':prior[f],**info,
            'alpha':state.alpha,'variance_floor':state.floor,'reference_universe':list(state.universe),'complement_reference':list(state.complement)})
        inv.extend({'fold':f,**row} for row in invariance(state,qp,prior[f]));allstates.append((f,tr,q,state,choices,info))
        b.emit('reference_outer_frozen',fold=f,choices=choices,seconds=time.monotonic()-tick,outer_scored=False)
    np.savez_compressed(out/'PREDICTIONS.npz',foldids=foldids,pairs=np.array(task.pairs),**pred)
    b.dump(out/'ALL_PREDICTIONS_FROZEN.json',{'utc':b.stamp(),'sha256':b.sha(out/'PREDICTIONS.npz'),
        'arrays':{k:hashlib.sha256(v.tobytes()).hexdigest() for k,v in pred.items()},'outer_scored':False})
    truth=np.full(len(task.pairs),np.nan)
    for f,tr,q,state,choices,info in allstates:truth[q]=state.h(task.d.y[q]).ravel()
    np.savez_compressed(out/'EVALUATION_LABELS.npz',truth=truth,foldids=foldids,pairs=np.array(task.pairs),x=task.d.x)
    rows=[];foldrows=[]
    for name,pp in pred.items():
        total=np.zeros(3)
        for f,tr,q,state,choices,info in allstates:
            P=b.projection(task.d.x[q]);err=pp[q]-truth[q];ss=np.array([np.sum(err**2),np.sum((P@err)**2),np.sum(((np.eye(len(q))-P)@err)**2)])
            total+=ss;foldrows.append({'method':name,'fold':f,'full_sse':ss[0],'P_sse':ss[1],'endpoint_sse':ss[2]})
        rows.append({'method':name,'full_sse':total[0],'P_sse':total[1],'endpoint_sse':total[2]})
    ref={v['method']:v for v in rows}
    for row in rows:
        for metric in ('full_sse','P_sse','endpoint_sse'):row[metric+'_gain_vs_saturation']=1-row[metric]/ref['scalar_saturation'][metric]
    fref={(v['method'],v['fold']):v for v in foldrows}
    for row in foldrows:
        for metric in ('full_sse','P_sse','endpoint_sse'):row[metric+'_gain_vs_saturation']=1-row[metric]/fref[('scalar_saturation',row['fold'])][metric]
    b.csv_write(out/'RESULTS.csv',rows);b.csv_write(out/'FOLD_RESULTS.csv',foldrows);b.dump(out/'INVARIANCE.json',inv)
    # Full outer fitting/inner selection rerun after excluded labels are poisoned.
    poison=[];rng=np.random.default_rng(83072)
    for f,tr,q,state,choices,info in allstates:
        yy=task.d.y[q].copy();vv={task.pairs[i]:task.values[task.pairs[i]].copy() for i in q};hh={task.pairs[i]:task.halves[task.pairs[i]] for i in q} if hasattr(task,'halves') else None
        task.d.y[q]=rng.normal(0,1e8,task.d.y[q].shape)
        for i in q:
            p=task.pairs[i];task.values[p]=rng.normal(0,1e8,task.values[p].shape)
            if hh is not None:task.halves[p]=tuple(rng.normal(0,1e8,z.shape) for z in hh[p])
        other,cc,ii=fit_outer(task,tr,f,pi,gm,b,geo)
        task.d.y[q]=yy;task.values.update(vv)
        if hh is not None:task.halves.update(hh)
        qp=[task.pairs[i] for i in q]
        difference=max(float(np.max(abs(state.predict(qp,lam,m)-other.predict(qp,lam,m)))) for m in METHODS for lam in LAMBDAS)
        samechoices=choices==cc;score_diff=max(abs(info['scores'][m][lam]-ii['scores'][m][lam]) for m in METHODS for lam in LAMBDAS)
        poison.append({'fold':f,'prediction_max_abs':difference,'inner_score_max_abs':score_diff,'choices_identical':samechoices,
            'excluded_query_means_poisoned':True,'excluded_query_halves_poisoned':hh is not None,
            'halves_not_used_by_learner':True,'PASS':difference==0 and samechoices and score_diff==0})
    b.dump(out/'LABEL_POISON_CHECK.json',poison)
    # Exact endpoint subspace algebra and candidate contrasts.
    algebra=[]
    for f,tr,q,state,choices,info in allstates:
        qp=[task.pairs[i] for i in q];x=task.d.x[q];P=b.projection(x);lam=prior[f];raw=state.raw(qp,lam);global_g=state.correction(qp,lam,'global');query_g=state.correction(qp,lam,'query');comp=state.correction(qp,lam,'complement')
        xf,_=state.design(universe);gf=state.correction(universe,lam,'global')
        algebra.append({'fold':f,'matched_lambda':lam,'P_global_vs_raw_max_abs':float(np.max(abs(P@global_g-P@raw))),
            'P_global_vs_query_max_abs':float(np.max(abs(P@global_g-query_g))),
            'complement_equals_outer_query_max_abs':float(np.max(abs(comp-query_g))),
            'global_universe_orthogonality_max_abs':float(np.max(abs(xf.T@gf))),
            'query_endpoint_raw_norm':float(np.linalg.norm((np.eye(len(q))-P)@raw)),
            'query_endpoint_global_norm':float(np.linalg.norm((np.eye(len(q))-P)@global_g)),
            'global_reference_beta_norm':float(np.linalg.norm(state.beta_global[lam])),
            'query_reference_beta_norm':float(np.linalg.norm(np.linalg.pinv(x)@raw)),
            'complement_reference_size':len(state.complement)})
    b.dump(out/'ALGEBRA.json',algebra)
    b.dump(out/'COMPLETE.json',{'utc':b.stamp(),'seconds':time.monotonic()-tick,'status':'PASS' if all(v['PASS'] for v in poison) and all(v['PASS'] for v in inv if v['stable_fixed_reference_expected']) else 'FAIL',
        'new_native_RNA_primary_claim':False,'independent_confirmation':False,'all_predictions_scored_after_freeze':True})
    b.emit('reference_complete',seconds=time.monotonic()-tick,results=rows[:10])

if __name__=='__main__':main()
