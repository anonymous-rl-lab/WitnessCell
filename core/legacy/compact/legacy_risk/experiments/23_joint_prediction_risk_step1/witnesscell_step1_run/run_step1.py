#!/usr/bin/env python3
"""Isolated training/scoring runner. Evaluation is a separate subcommand."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
import argparse,csv,hashlib,json,pickle,sys,time,resource
from pathlib import Path
from collections import Counter
import numpy as np

HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE.parent/'witnesscell_step1_protocol_v1/protocol.json').read_text())
PID=CFG['protocol_id']
REPO=Path(os.environ.get('WITNESS_REPO',str(HERE.parents[1]/'grn_extracted/grn_repo_v18')))
SRC=REPO/'experiments/19_v14_incremental_amplitude_gate/src'
sys.path.insert(0,str(SRC))
import witnesscell_v4_core as core
import scipy
from fast_amplitude_bank import install as install_fast_bank

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def h(*parts):return hashlib.sha256('|'.join(map(str,parts)).encode()).hexdigest()
def pair(c):return '+'.join(sorted(c.split('+')))
def emit(**x):print(json.dumps(x,allow_nan=False),flush=True)
def dump(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=True)+'\n')
def csvwrite(path,rows):
    rows=list(rows)
    if not rows:raise ValueError('refuse empty records')
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def weighted_q(x,w,q):
    order=np.argsort(x,kind='stable');a=np.asarray(x)[order];b=np.asarray(w)[order]
    return float(a[min(np.searchsorted(np.cumsum(b),q*b.sum(),side='left'),len(a)-1)])

def prepare():
    """Partition labels into allowed fitting inputs and separate evaluation files."""
    for ds in ['Norman','Wessels']:
        with np.load(REPO/f'experiments/18_dual_head_evidence_gate/cache/{ds}.condition_moments.npz',allow_pickle=False) as z:
            z={k:z[k] for k in z.files}
        names=z['conditions'].astype(str).tolist();idx={c:i for i,c in enumerate(names)}
        ledger=json.loads((REPO/f'experiments/22_metric_calibration_stress/assets/formal_contracts/{ds}.split_ledger.json').read_text())['split_contract']
        with np.load(REPO/f'experiments/22_metric_calibration_stress/assets/formal_contracts/{ds}.phase_m_contract.npz',allow_pickle=False) as zc:
            weight_names=zc['conditions'].astype(str).tolist();weights=zc['weights'];weight_genes=zc['genes']
        assert np.array_equal(weight_genes,z['genes'])
        wi={c:i for i,c in enumerate(weight_names)}
        for seed in [1,2,3]:
            out=HERE/'results'/f'{ds}_seed{seed}';out.mkdir(parents=True,exist_ok=True)
            sp=ledger[str(seed)];train=sp['train'];val=sp['validation'];alltest=sp['test'];test=[c for c in alltest if '+' in c]
            allowed=sorted(set(['control']+train+val));assert not set(alltest)&set(allowed)
            assert len(test)==CFG['test_double_counts'][ds][seed-1]
            with np.load(REPO/f'experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo/{ds}/seed{seed}/deploy_predictions.npz',allow_pickle=False) as d:
                assert np.array_equal(d['genes'],z['genes']);di={c:i for i,c in enumerate(d['conditions'].astype(str))}
                pred=np.stack([d['prediction'][di[c]] for c in test]);fact=np.stack([d['factorized_prediction'][di[c]] for c in test])
            arows=[idx[c] for c in allowed]
            np.savez_compressed(out/'allowed_input.npz',conditions=np.array(allowed),means=z['means'][arows],variances=z['variances'][arows],counts=z['counts'][arows],genes=z['genes'],train=np.array(train),validation=np.array(val),targets=np.array(test),frozen_prediction=pred,frozen_factorized=fact)
            # No evaluation statistic is calculated here; this file is excluded from fitting/scoring.
            np.savez_compressed(out/'evaluation_only.npz',conditions=np.array(test),truth=z['means'][[idx[c] for c in test]],weights=np.stack([weights[wi[c]] for c in test]),genes=z['genes'])
            dump(out/'split_manifest.json',{'dataset':ds,'seed':seed,'split':sp,'allowed_input_sha256':sha(out/'allowed_input.npz'),'evaluation_only_sha256':sha(out/'evaluation_only.npz'),'risk_forbidden_fields':['test_cell_counts','test_variances','test_truth','test_DE_weights']})
    emit(event='prepared',splits=6,test_truth_scored=False)

_go=None;_go_sha=None;_orig_head=core.fit_incremental_amplitude_head
_projection_mode=False
_head_cache_dir=HERE/'identity_cache'
def head_key(kwargs):
    d=hashlib.sha256();d.update(_go_sha.encode())
    for p in sorted(SRC.glob('*.py')):d.update(sha(p).encode())
    d.update(sha(HERE/'fast_amplitude_bank.py').encode())
    d.update((sys.version+np.__version__+scipy.__version__).encode())
    for k,v in kwargs.items():
        if k=='gene2go':continue
        d.update(k.encode())
        if isinstance(v,dict):
            for name,value in v.items():
                d.update(str(name).encode());a=np.asarray(value);d.update(str(a.dtype).encode());d.update(str(a.shape).encode());d.update(a.tobytes())
        elif k=='genes':d.update(json.dumps(list(v)).encode())
        else:
            a=np.asarray(v);d.update(str(a.dtype).encode());d.update(str(a.shape).encode());d.update(a.tobytes())
    return d.hexdigest()
def cached_head(**kwargs):
    key=head_key(kwargs)
    if _projection_mode:key='projected_'+h(key,sha(HERE/'run_step1.py'))
    p=_head_cache_dir/f'{key}.pkl'
    if p.exists():
        with p.open('rb') as f:return pickle.load(f)
    t=time.monotonic()
    if _projection_mode:
        # Exact projection: v14's incremental correction changes only each
        # query endpoint coordinate. Hyperparameter selection and residual
        # fitting already use baseline_head in the original v14 core.
        from multihead_v14 import fit_frozen_v13,IncrementalAmplitudeHead
        baseline=fit_frozen_v13(kwargs['single_effects'],kwargs['single_means'],kwargs['single_variances'],kwargs['single_counts'],kwargs['control_mean'],kwargs['control_variance'],kwargs['control_count'],kwargs['genes'],kwargs['gene2go'])
        obj=IncrementalAmplitudeHead(baseline,False,(0.,0.,0.),0.,float('-inf'),0.,float('-inf'),0.,float('-inf'),len(kwargs['single_effects']),())
    else:obj=_orig_head(**kwargs)
    _head_cache_dir.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(f'.{os.getpid()}.tmp')
    with tmp.open('wb') as f:pickle.dump(obj,f)
    tmp.replace(p)
    emit(event='identity_fit',seconds=round(time.monotonic()-t,3),known_singles=len(kwargs['single_effects']),projection=_projection_mode)
    return obj
def load_go():
    global _go,_go_sha
    if _go is None:
        p=REPO/'experiments/15_scperturbench_sota/module/data/gears_assets/gene2go_all.pkl'
        _go_sha=sha(p)
        with p.open('rb') as f:_go=pickle.load(f)
        install_fast_bank()
        core.fit_incremental_amplitude_head=cached_head
    return _go

class Split:
    def __init__(self,ds,seed):
        self.ds=ds;self.seed=seed;self.out=HERE/'results'/f'{ds}_seed{seed}';self.out.mkdir(parents=True,exist_ok=True)
        with np.load(self.out/'allowed_input.npz',allow_pickle=False) as z:
            z={k:z[k] for k in z.files}
            names=z['conditions'].astype(str).tolist();self.genes=z['genes'].astype(str).tolist()
            self.means={c:z['means'][i].astype(float) for i,c in enumerate(names)}
            self.variances={c:z['variances'][i].astype(float) for i,c in enumerate(names)}
            self.counts={c:int(z['counts'][i]) for i,c in enumerate(names)}
            self.train=z['train'].astype(str).tolist();self.val=z['validation'].astype(str).tolist();self.targets=z['targets'].astype(str).tolist()
            self.pred=z['frozen_prediction'].astype(float);self.fact=z['frozen_factorized'].astype(float)
        self.ctrl=self.means['control'];self.gi={g:i for i,g in enumerate(self.genes)};self.go=load_go()
        self.source_digest=h(*[sha(p) for p in sorted(SRC.glob('*.py'))],*[sha(HERE/f) for f in ['run_step1.py','fast_amplitude_bank.py','baselines.py']],sys.version,np.__version__,scipy.__version__)
        self.input_digest=sha(self.out/'allowed_input.npz')
        assert not set(self.targets)&set(self.means)
    def fit(self,train,val,targets,tag,want_object=False,projected=True):
        global _projection_mode
        _projection_mode=projected
        allowed=set(['control']+train+val)
        assert not allowed&set(targets),(tag,'label_overlap')
        keys=sorted(allowed);fm={c:self.means[c] for c in keys};fv={c:self.variances[c] for c in keys};fc={c:self.counts[c] for c in keys}
        digest=h(self.source_digest,self.input_digest,json.dumps([train,val,targets,projected]))
        p=self.out/'fits'/f'{tag}.npz';meta=p.with_suffix('.json')
        if p.exists() and not want_object:
            md=json.loads(meta.read_text());assert md['input_digest']==digest
            with np.load(p,allow_pickle=False) as z:return z['prediction'].astype(float),None
        emit(event='fit_start',dataset=self.ds,seed=self.seed,tag=tag,train_singles=sum(c!='control' and '+' not in c for c in train),train_doubles=sum('+' in c for c in train),queries=len(targets))
        t=time.monotonic()
        obj=core.fit_predict(fm,fv,fc,self.genes,train,val,targets,identity_mode='multihead_v14',gene2go=self.go)
        prediction=np.stack([obj.predictions[c] for c in targets]).astype(np.float32)
        if not np.all(np.isfinite(prediction)):raise ValueError('nonfinite predictions')
        p.parent.mkdir(parents=True,exist_ok=True)
        valid=np.stack([self.mask(c) for c in targets]) if projected else np.ones_like(prediction,dtype=bool)
        np.savez_compressed(p,prediction=prediction,conditions=np.array(targets),valid_coordinates=valid)
        md={'tag':tag,'input_digest':digest,'fit_seconds':time.monotonic()-t,'alpha':obj.alpha,'noise':obj.noise_ratio,'gamma':obj.gamma,'allowed_conditions':keys,'train':train,'validation':val,'targets':targets,'valid_coordinates':'non_target_only' if projected else 'all','maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        dump(meta,md);emit(event='fit_done',dataset=self.ds,seed=self.seed,tag=tag,seconds=round(md['fit_seconds'],3))
        return prediction.astype(float),obj
    def mask(self,c):
        m=np.ones(len(self.genes),dtype=bool)
        for g in c.split('+'):m[self.gi[g]]=False
        return m
    def subsets(self):
        doubles=[c for c in self.train if '+' in c]
        for scheme in ['pair','endpoint']:
            ids=[pair(c) for c in doubles] if scheme=='pair' else sorted({g for c in self.train if c!='control' for g in c.split('+')})
            ordered=sorted(set(ids),key=lambda x:h(PID,self.ds,self.seed,scheme,x));assignment={v:i%5 for i,v in enumerate(ordered)}
            for fold in range(5):
                held={x for x,v in assignment.items() if v==fold}
                if scheme=='pair':
                    queries=[c for c in doubles if pair(c) in held];train=[c for c in self.train if c not in queries];val=list(self.val)
                else:
                    touches=lambda c:c!='control' and bool(set(c.split('+'))&held)
                    queries=[c for c in doubles if touches(c)];train=[c for c in self.train if not touches(c)];val=[c for c in self.val if not touches(c)]
                yield scheme,fold,train,val,queries
    def ensemble(self,train,val,targets,scheme,fold):
        predictions=[]
        for member in range(8):
            rng=np.random.default_rng(int(h(PID,self.ds,self.seed,scheme,fold,member)[:8],16))
            selected={'control'}
            for cs in [[c for c in train if c!='control' and '+' not in c],[c for c in train if '+' in c]]:
                if cs:selected.update(rng.choice(cs,size=max(1,int(.8*len(cs))),replace=False).tolist())
            tr=[c for c in train if c in selected]
            pred,_=self.fit(tr,val,targets,f'{scheme}_{fold}_bag{member}')
            predictions.append(pred)
        return np.var(np.stack(predictions),axis=0,ddof=1)
    def feature_rows(self,targets,pred,var,train,tag,labels=False):
        known={c for c in train if c!='control' and '+' not in c};deg=Counter(g for c in train if '+' in c for g in c.split('+'));rows=[]
        for i,c in enumerate(targets):
            m=self.mask(c);es=float(np.mean((pred[i,m]-self.ctrl[m])**2));vs=float(np.mean(var[i,m]));ep=c.split('+')
            row={'dataset':self.ds,'seed':self.seed,'pair':pair(c),'condition':c,'scenario':tag,'effect_score':float(np.log1p(es)),'variance_score':float(np.log1p(vs)),'missing_endpoints':sum(g not in known for g in ep),'inverse_degree':float(np.mean([1/(1+deg[g]) for g in ep])),'seen':sum(g in known for g in ep)}
            if labels:
                row['loss_mse']=float(np.mean((pred[i,m]-self.means[c][m])**2));row['true_effect']=float(np.mean((self.means[c][m]-self.ctrl[m])**2))
            rows.append(row)
        return rows
    def warm_replay(self):
        # No external replay pickle is adopted as an identity-cache entry.
        return
    def replay(self):
        self.warm_replay();pred,obj=self.fit(self.train,self.val,self.targets,'full_replay',want_object=True,projected=False)
        diff=float(np.max(np.abs(pred-self.pred)));match=bool(np.allclose(pred,self.pred,atol=1e-6,rtol=1e-6));exact=bool(np.array_equal(pred,self.pred))
        dump(self.out/'replay.json',{'max_abs_error':diff,'within_tolerance':match,'float32_exact':exact})
        if not match:raise ValueError('final mean replay mismatch')
        from baselines import build_baselines
        b=build_baselines({c:self.means[c] for c in set(['control']+self.train)},{c:self.means[c] for c in set(['control']+self.val)},self.genes,self.train,self.val,self.targets,core_fit=obj,factorized_prediction=self.fact,gene2go=self.go)
        np.savez_compressed(self.out/'baseline_predictions.npz',conditions=np.array(self.targets),**b.predictions)
        dump(self.out/'baseline_metadata.json',dict(b.metadata,saturation_alpha=b.saturation_alpha,validation_scores=b.validation_scores))
    def run(self,pilot=False):
        self.warm_replay();allrows=[];t=time.monotonic()
        for scheme,fold,train,val,queries in self.subsets():
            tag=f'{scheme}_{fold}'
            if not queries:continue
            if pilot and (scheme!='pair' or fold!=0):continue
            pred,_=self.fit(train,val,queries,tag+'_base');var=self.ensemble(train,val,queries,scheme,fold)
            rows=self.feature_rows(queries,pred,var,train,tag,labels=True);csvwrite(self.out/f'oof_{tag}.csv',rows);allrows.extend(rows)
        if pilot:
            dump(self.out/'pilot.json',{'seconds':time.monotonic()-t,'oof_queries':len(allrows),'fit_calls':9,'test_truth_scored':False});emit(event='pilot_done',seconds=time.monotonic()-t);return
        counts=Counter(r['pair'] for r in allrows)
        for row in allrows:row['weight']=1/counts[row['pair']]
        csvwrite(self.out/'training_oof_records.csv',allrows)
        weights=np.array([r['weight'] for r in allrows]);weights/=weights.sum()
        losses=np.array([r['loss_mse'] for r in allrows]);m=max(weighted_q(losses,weights,.5),1e-12)
        y=np.log1p(losses/m);names=['effect_score','variance_score','missing_endpoints','inverse_degree']
        x=np.array([[r[n] for n in names] for r in allrows]);models={}
        for name,cols in [('B',2),('S',4)]:
            xx=x[:,:cols];mu=weights@xx;sd=np.sqrt(weights@((xx-mu)**2));sd[sd==0]=1.;xx=(xx-mu)/sd;ym=float(weights@y)
            beta=np.linalg.solve(xx.T@(weights[:,None]*xx)+np.eye(cols),xx.T@(weights*(y-ym)))
            models[name]={'features':names[:cols],'mean':mu.tolist(),'sd':sd.tolist(),'coef':beta.tolist(),'intercept':ym,'m':m}
        dump(self.out/'risk_models.json',models)
        self.replay()
        var=self.ensemble(self.train,self.val,self.targets,'full',0)
        rows=self.feature_rows(self.targets,self.pred,var,self.train,'formal',labels=False)
        for row in rows:
            for name,model in models.items():
                xx=np.array([row[n] for n in model['features']]);row['risk_'+name]=float(((xx-model['mean'])/model['sd'])@model['coef']+model['intercept'])
        csvwrite(self.out/'frozen_query_outputs.csv',rows)
        eps=max(weighted_q([r['true_effect'] for r in allrows],weights,.1),1e-12)
        cuts=[weighted_q(x[:,0],weights,q) for q in [.25,.5,.75]]
        dump(self.out/'evaluation_metadata.json',{'dataset':self.ds,'seed':self.seed,'target_genes_removed':True,'loss_truth':'E18_full_cell_condition_means','oof_effect_quartile_edges':cuts,'oof_effect_quartile_axis':'effect_score','effect_quartile_transform':'log1p','relative_loss_epsilon':eps,'frozen_query_sha256':sha(self.out/'frozen_query_outputs.csv'),'risk_model_sha256':sha(self.out/'risk_models.json'),'total_oof_records':len(allrows),'unique_training_pairs':len(counts),'formal_targets':len(rows),'elapsed_seconds':time.monotonic()-t})
        emit(event='split_complete',dataset=self.ds,seed=self.seed,seconds=round(time.monotonic()-t,2),queries=len(rows),test_truth_scored=False)

def evaluate_labels(ds,seed):
    out=HERE/'results'/f'{ds}_seed{seed}'
    md=json.loads((out/'evaluation_metadata.json').read_text());assert sha(out/'frozen_query_outputs.csv')==md['frozen_query_sha256']
    with (out/'frozen_query_outputs.csv').open() as f:rows=list(csv.DictReader(f))
    with np.load(out/'allowed_input.npz',allow_pickle=False) as z:
        genes=z['genes'].astype(str);ctrl=z['means'][list(z['conditions'].astype(str)).index('control')].astype(float);pred=z['frozen_prediction'].astype(float)
        training_names=[c for c in z['train'].astype(str).tolist() if c!='control'];ii={c:i for i,c in enumerate(z['conditions'].astype(str))};training_mean=np.mean(z['means'][[ii[c] for c in training_names]].astype(float),axis=0)
    with np.load(out/'evaluation_only.npz',allow_pickle=False) as z:truth=z['truth'].astype(float);w=z['weights'].astype(float);conditions=z['conditions'].astype(str)
    with np.load(out/'baseline_predictions.npz',allow_pickle=False) as z:baselines={k:z[k].astype(float) for k in z.files if k!='conditions'}
    gi={g:i for i,g in enumerate(genes)};baseline_rows=[]
    def losses(p,t,ww,mask):
        err=(p[mask]-t[mask])**2;tw=ww[mask];sw=tw.sum();mse=float(err.mean());wmse=float((tw@err)/sw) if sw>0 else float('nan')
        de=t[mask]-ctrl[mask];pe=p[mask]-ctrl[mask];effect=float(np.mean(de**2));den=float(np.linalg.norm(de)*np.linalg.norm(pe));cos=1-float(np.clip(float(de@pe)/den,-1.,1.)) if den>0 else float('nan')
        inds=np.flatnonzero(mask);top=sorted(inds,key=lambda j:(-ww[j],str(genes[j])))[:100]
        d=t[mask]-training_mean[mask];dp=p[mask]-training_mean[mask]
        denom=float(tw@((d-np.average(d,weights=tw))**2)) if sw>0 else 0
        r2=1-float(tw@((d-dp)**2))/denom if denom>0 else float('nan')
        return dict(loss_mse=mse,loss_wmse=wmse,loss_relative=mse/max(effect,md['relative_loss_epsilon']),loss_cosine=cos,true_effect=effect,loss_top100=float(np.mean((p[top]-t[top])**2)),weighted_delta_r2=r2)
    for i,row in enumerate(rows):
        assert row['condition']==conditions[i]
        mask=np.ones(len(genes),dtype=bool)
        for g in row['condition'].split('+'):mask[gi[g]]=False
        row.update(losses(pred[i],truth[i],w[i],mask))
        full=losses(pred[i],truth[i],w[i],np.ones(len(genes),dtype=bool))
        row.update({k+'_including_targets':v for k,v in full.items()})
        for name,mat in {'WitnessCell':pred,**baselines}.items():
            br={k:row[k] for k in ['dataset','seed','pair','condition','seen']};br['method']=name;br.update(losses(mat[i],truth[i],w[i],mask));baseline_rows.append(br)
    csvwrite(out/'scored_queries.csv',rows);csvwrite(out/'baseline_losses.csv',baseline_rows)
    emit(event='labels_evaluated',dataset=ds,seed=seed,queries=len(rows))

def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','pilot','formal','replay','labels']);p.add_argument('--dataset',default='Norman');p.add_argument('--seed',type=int,default=1);a=p.parse_args()
    if a.mode=='prepare':prepare()
    elif a.mode=='labels':evaluate_labels(a.dataset,a.seed)
    else:
        s=Split(a.dataset,a.seed)
        if a.mode=='replay':s.replay()
        else:s.run(pilot=a.mode=='pilot')
if __name__=='__main__':main()
