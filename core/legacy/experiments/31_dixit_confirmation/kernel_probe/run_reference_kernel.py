#!/usr/bin/env python3
"""Fixed-reference kernel constraint, with full-loss and train-P-loss KRR."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
import argparse,hashlib,importlib.util,itertools,json,sys,time
from pathlib import Path
import numpy as np

def mod(name,path):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sys.modules[name]=m;sp.loader.exec_module(m);return m

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--asset-root',type=Path);ap.add_argument('--compact-root',type=Path)
    ap.add_argument('--compact-loader',type=Path,default=Path(__file__).resolve().parent.parent/'compact_loader.py')
    ap.add_argument('--model-code',type=Path);ap.add_argument('--bridge-code',type=Path,required=True)
    ap.add_argument('--geometry-code',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.compact_root is None and a.asset_root is None:ap.error('Supply --compact-root or --asset-root')
    if a.model_code is None and a.asset_root is None:ap.error('Supply --model-code when using compact-only input')
    b=mod('refkernel_bridge',a.bridge_code);geo=mod('refkernel_geometry',a.geometry_code)
    r=a.asset_root;gm=b.module('global_pair_model',a.model_code or r/'experiment_assets/WitnessCell_Partner_Component_Run_20261018/model/global_pair_model.py')
    if a.compact_root is None:
        old=b.module('refkernel_old',r/'source_v19/WitnessCell_Paper_v19_Full/source_data/RNA_controls/run_controls.py');new=b.module('refkernel_new',r/'work_B/run_frozen_rna_confirmation.py')
        task=b.load_task('Dixit',r,gm,new,old);program_path=r/'work_B/PREDEFINED_G2M_PROGRAM.json'
    else:
        loader=mod('refkernel_compact',a.compact_loader);task=loader.load_task(a.compact_root,gm,b)
        program_path=a.compact_root/'PREDEFINED_G2M_PROGRAM.json'
    out=a.out;out.mkdir(parents=True,exist_ok=True)
    if (out/'STARTED.json').exists():raise FileExistsError(out)
    tick=time.monotonic();here=Path(__file__).resolve().parent
    program=json.loads(program_path.read_text());pi=np.array([int(np.flatnonzero(task.symbols==s)[0]) for s in program['genes'] if np.count_nonzero(task.symbols==s)==1])
    assert len(pi)==52
    refs=['+'.join(p) for p in itertools.combinations(sorted(task.nodes),2)];refidx={p:i for i,p in enumerate(refs)}
    ids=np.array([refidx[p] for p in task.pairs]);lut={p:i for i,p in enumerate(task.pairs)};n=len(task.pairs)
    assert len(refs)==45 and sorted(refs)==sorted(task.pairs)
    xr=gm.incidence(refs,task.nodes);x=task.d.x;P=b.projection(xr)
    edges=np.array([[task.nodes.index(g) for g in p.split('+')] for p in refs]);e=np.stack([task.values[g]-task.values['control'] for g in task.nodes])
    K=gm.symmetric_kernel(edges,e@e.T/e.shape[1]);H=P@K@P;H=(H+H.T)/2
    b.dump(out/'STARTED.json',{'utc':b.stamp(),'plan_sha256':b.sha(here/'PLAN.txt'),'script_sha256':b.sha(__file__),
        'primary':'ref_full','secondary':'ref_trainP','selection':'cycle3 inner FULL program SSE','program_symbols':task.symbols[pi].tolist(),
        'reference_pairs':refs,'reference_rule':'all distinct unordered allowed measured-single node pairs','reference_rank':int(np.linalg.matrix_rank(xr)),
        'kernel_sha256':hashlib.sha256(K.tobytes()).hexdigest(),'reference_kernel_sha256':hashlib.sha256(H.tobytes()).hexdigest(),
        'reference_kernel_orthogonality_max_abs':float(np.max(np.abs(xr.T@H))),'status':'post-observation retrospective development'})
    preds={};foldids=np.full(n,-1);frozen=[];audits=[]
    def fill(name,q,v):
        if name not in preds:preds[name]=np.full(n,np.nan)
        preds[name][q]=np.asarray(v).ravel()
    def fit(fit,identity):
        allowed=['control',*task.nodes,*[task.pairs[i] for i in fit]];matrix=np.stack([task.values[v] for v in allowed]);var=np.var(matrix,axis=0,ddof=1)
        floor=float(np.quantile(var[var>0],.01));sigma=np.sqrt(np.maximum(var,floor))
        def h(v):return np.mean(np.asarray(v)[...,pi]/sigma[pi],axis=-1,keepdims=True)
        add=h(task.d.additive);target=h(task.d.y[fit]);alpha=min(gm.ALPHAS,key=lambda z:(float(np.sum((gm.sat(add[fit],z)-target)**2)),z))
        residual=target-gm.sat(add[fit],alpha);t=ids[fit];ht=H[np.ix_(t,t)];pt=b.projection(x[fit]);G={};aa=[]
        for kind in ('ref_full','ref_trainP'):
            ptuse=np.eye(len(fit)) if kind=='ref_full' else pt
            A=ptuse@ht@ptuse;A=(A+A.T)/2;scale=float(np.trace(A)/len(fit))
            if scale<=1e-14:raise ValueError('Degenerate reference kernel')
            normalized=A/scale;ev,V=np.linalg.eigh(normalized);cross=H[:,t]@ptuse/scale
            G[kind]={};solves=[]
            for lam in b.LAMBDAS:
                coef=((V/(ev+lam))@V.T)@residual;G[kind][lam]=cross@coef
                normalres=normalized@coef+lam*coef-residual
                solves.append({'lambda':lam,'linear_system_residual_norm':float(np.linalg.norm(normalres)),
                    'global_endpoint_orthogonality_max_abs':float(np.max(np.abs(xr.T@G[kind][lam])))})
            aa.append({'kind':kind,'scale':scale,'min_eigenvalue':float(ev.min()),'solves':solves})
        audit={'identity':identity,'fit':fit.tolist(),'normalization_conditions':allowed,'variance_floor':floor,'alpha':alpha,'fits':aa}
        audits.append(audit)
        return {'sigma':sigma,'alpha':alpha,'G':G,'add':add,'h':h,'audit':audit}
    def request(state,query,kind,lam):
        # Request composition is deliberately absent from all fitted quantities.
        return gm.sat(state['add'][query],state['alpha'])+state['G'][kind][lam][ids[query]]
    for fold in task.folds:
        f=fold['fold'];tr=np.array([lut[p] for p in fold['train_pairs']]);q=np.array([lut[p] for p in fold['test_pairs']])
        scores={kind:{lam:0. for lam in b.LAMBDAS} for kind in ('ref_full','ref_trainP')};scoresP={kind:{lam:0. for lam in b.LAMBDAS} for kind in scores};zero=0.;inner=[]
        for j,held in enumerate(geo.graph_folds(tr,task.pairs,f,x,'cycle3')):
            ti=np.array(sorted(set(tr)-set(held)));state=fit(ti,f'outer{f}:inner{j}');truth=state['h'](task.d.y[held]);base=gm.sat(state['add'][held],state['alpha']);zero+=float(np.sum((base-truth)**2));pv=b.projection(x[held])
            for kind in scores:
                for lam in b.LAMBDAS:
                    er=request(state,held,kind,lam)-truth;scores[kind][lam]+=float(np.sum(er**2));scoresP[kind][lam]+=float(np.sum((pv@er)**2))
            inner.append({'fit':ti.tolist(),'validation':held.tolist(),'Pdim':len(held)-int(np.linalg.matrix_rank(x[held])),'alpha':state['alpha']})
        choices={kind:min(scores[kind],key=lambda z:(scores[kind][z],-z)) for kind in scores};state=fit(tr,f'outer{f}:final')
        base=gm.sat(state['add'][q],state['alpha']);fill('scalar_saturation',q,base);fill('measured_additive',q,state['add'][q]);foldids[q]=f
        for kind in scores:
            for lam in b.LAMBDAS:fill(f'{kind}_fixed_{lam:g}',q,request(state,q,kind,lam))
            lam=choices[kind];fill(kind,q,request(state,q,kind,lam));fill(kind+'_gated',q,request(state,q,kind,lam) if scores[kind][lam]<=zero else base)
        frozen.append({'fold':f,'tr':tr,'q':q,'state':state,'choices':choices})
        b.dump(out/f'FOLD_{f}.json',{'fold':f,'choices':choices,'inner_FULL':scores,'inner_P':scoresP,'zero_FULL':zero,'inner_folds':inner,'alpha':state['alpha'],
            'inner_noharm_gates':{kind:scores[kind][choices[kind]]<=zero for kind in scores}})
        np.savez_compressed(out/f'MODEL_{f}.npz',reference_pairs=np.array(refs),query=q,train=tr,sigma=state['sigma'],alpha=state['alpha'],additive_program=state['add'],
            **{f'{kind}_lambda_{lam:g}':g for kind,gg in state['G'].items() for lam,g in gg.items()})
        b.emit('reference_kernel_fold_frozen',outer=f,seconds=time.monotonic()-tick,outer_scored=False)
    assert all(np.isfinite(v).all() for v in preds.values())
    np.savez_compressed(out/'PREDICTIONS.npz',foldids=foldids,**preds)
    b.dump(out/'ALL_PREDICTIONS_FROZEN.json',{'utc':b.stamp(),'sha256':b.sha(out/'PREDICTIONS.npz'),'all_predictions_frozen_before_outer_evaluation':True})
    truth=np.full(n,np.nan)
    for item in frozen:truth[item['q']]=item['state']['h'](task.d.y[item['q']]).ravel()
    rows=[];byfold=[]
    for name,pred in preds.items():
        er=pred-truth;full=float(np.sum(er**2));psse=0.
        for f in range(3):
            q=np.flatnonzero(foldids==f);pp=float(np.sum((b.projection(x[q])@er[q])**2));psse+=pp;byfold.append({'method':name,'fold':f,'full_sse':float(np.sum(er[q]**2)),'P_sse':pp})
        rows.append({'method':name,'full_sse':full,'P_sse':psse})
    lookup={z['method']:z for z in rows}
    for row in rows:
        for metric in ('full_sse','P_sse'):
            for base in ('scalar_saturation','measured_additive'):row[metric+'_gain_vs_'+base]=1-row[metric]/lookup[base][metric]
    b.csv_write(out/'RESULTS.csv',rows);b.csv_write(out/'FOLD_RESULTS.csv',byfold);b.dump(out/'SOLVER_AUDIT.json',audits)
    np.savez_compressed(out/'EVALUATION_TRUTH.npz',truth=truth,foldids=foldids,pairs=np.array(task.pairs),program_indices=pi)
    invariant=[]
    for item in frozen:
        q=item['q'];state=item['state']
        for kind,lam in item['choices'].items():
            whole=request(state,q,kind,lam);single=np.vstack([request(state,np.array([v]),kind,lam) for v in q]);reverse=request(state,q[::-1],kind,lam)[::-1]
            deletes=[float(np.max(np.abs(request(state,np.delete(q,j),kind,lam)-np.delete(whole,j,axis=0)))) for j in range(len(q))]
            allreq=request(state,np.arange(n),kind,lam)[q];duplicate=request(state,np.repeat(q,2),kind,lam)[::2]
            correction=state['G'][kind][lam][ids[q]]
            invariant.append({'outer':item['fold'],'method':kind,'single_vs_whole_max_abs':float(np.max(np.abs(single-whole))),
                'reverse_vs_whole_max_abs':float(np.max(np.abs(reverse-whole))),'delete_one_max_abs':max(deletes),
                'add_all_queries_max_abs':float(np.max(np.abs(allreq-whole))),'duplicate_queries_max_abs':float(np.max(np.abs(duplicate-whole))),
                'singleton_nonzero_corrections':int(np.count_nonzero(np.abs(correction)>1e-12)),'query_correction_rms':float(np.sqrt(np.mean(correction**2)))})
    b.dump(out/'REQUEST_INVARIANCE.json',invariant)
    item=frozen[0];tr=item['tr'];q=item['q'];before=fit(tr,'poison_before');rng=np.random.default_rng(2026102017);save_y=task.d.y[q].copy();save_v={task.pairs[i]:task.values[task.pairs[i]].copy() for i in q}
    task.d.y[q]=rng.normal(0,1e7,task.d.y[q].shape)
    for i in q:task.values[task.pairs[i]]=rng.normal(0,1e7,len(task.genes))
    after=fit(tr,'poison_after');task.d.y[q]=save_y;task.values.update(save_v)
    dif=max(float(np.max(np.abs(before['G'][kind][lam]-after['G'][kind][lam]))) for kind in before['G'] for lam in b.LAMBDAS)
    b.dump(out/'LABEL_POISON.json',{'status':'PASS' if dif==0 else 'FAIL','correction_max_abs':dif,'alpha_equal':before['alpha']==after['alpha'],
        'sigma_max_abs':float(np.max(np.abs(before['sigma']-after['sigma']))),'outer_query_means_poisoned':True,'test_halves_never_used':True})
    b.dump(out/'COMPLETE.json',{'utc':b.stamp(),'seconds':time.monotonic()-tick,'primary':'ref_full','status':'complete exploratory fixed-reference model',
        'universe_size':len(refs),'known_measured_single_endpoints':len(task.nodes),'single_query_supported':True,'candidate_batch_independent':True})
    b.emit('reference_kernel_complete',seconds=time.monotonic()-tick)
if __name__=='__main__':main()
