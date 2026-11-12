from pathlib import Path
import sys,importlib.util,json,io,tarfile,hashlib,shutil,csv
import numpy as np
import pandas as pd
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--source-v24',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
R=args.source_v24.resolve()
O=args.out.resolve();O.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('restore',R/'tools/restore_data_v24.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
idx=json.loads((R/'reproduction_store/INDEX.json').read_text());store=m.Store(R,idx); origins={}
def read(rel):
    b=(R/rel).read_bytes() if (R/rel).exists() else b''.join(store.parts(idx['files'][rel]));origins[rel]={'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()};return b
def put(rel,b):
    p=O/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
def copy(rel):put(rel,read(rel))
def js(p,v):put(p,(json.dumps(v,indent=2,allow_nan=True)+'\n').encode())
def npz(rel):return np.load(io.BytesIO(read(rel)),allow_pickle=False)
def loadmod(name,path):
    sp=importlib.util.spec_from_file_location(name,path);x=importlib.util.module_from_spec(sp);sp.loader.exec_module(x);return x
E20='experiments/20_osf_preregistered_replication';E21='experiments/21_frozen_selective_prediction';E22='experiments/22_metric_calibration_stress';E23='experiments/23_joint_prediction_risk_step1';U=E23+'/witnesscell_step1_run'
allpaths=set(idx['files'])|{str(p.relative_to(R)) for p in (R/'experiments').rglob('*') if p.is_file()}
# retain source protocols and small formal outputs, omit obsolete smoke/pilot/cache trees
for rel in sorted(allpaths):
    if '__pycache__' in rel:continue
    keep=False
    if rel.startswith(E20+'/') and not rel.endswith(('.tar.gz','.tgz','.zip')):keep=True
    if rel.startswith(E20+'/') and rel.endswith('.zip'):keep=True
    if rel.startswith(E21+'/'):keep=True
    if rel.startswith(E22+'/'):
        keep=(('/src/' in rel or '/amendments/' in rel or '/audit/' in rel or '/results/formal_e22/' in rel or '/assets/formal_contracts/' in rel) and not rel.endswith('.npz')) or (len(Path(rel).relative_to(E22).parts)==1)
    if rel.startswith(E23+'/witnesscell_step1_protocol_v1/'):keep=True
    if rel.startswith(U+'/'):
        s=rel[len(U)+1:];parts=Path(s).parts
        if len(parts)==1 and (s.endswith('.py') or s in ['execution_lock.json','all_risks_frozen.json','all_selections_frozen.json','decision.json','README.md']):keep=True
        if len(parts)==3 and parts[0]=='results' and parts[2] in ['frozen_query_outputs.csv','scored_queries.csv','frozen_selections.csv','risk_models.json','training_oof_records.csv','evaluation_metadata.json','split_manifest.json','replay.json','baseline_losses.csv','baseline_metadata.json','full_replay.json']:keep=True
        if len(parts)==2 and parts[0]=='evaluation' and s.endswith(('.json','.csv')):keep=True
    if keep:copy(rel)
print('copied',len(origins),'small files',flush=True)
# E20 exact vector identity receipt, remove redundant tar
tr=[r for r in allpaths if r.startswith(E20+'/') and r.endswith(('.tar.gz','.tgz'))];print(tr,flush=True)
receipts=[]
for tarpath in tr:
    with tarfile.open(fileobj=io.BytesIO(read(tarpath)),mode='r:gz') as tf:
        for mem in tf.getmembers():
            if mem.name.endswith('deploy_predictions.npz'):
                comps=mem.name.split('/');i=comps.index('gate19_regenerated');ds,seed=comps[i+1:i+3]
                raw=tf.extractfile(mem).read();a=np.load(io.BytesIO(raw),allow_pickle=False)
                canonical=f'experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo/{ds}/{seed}/deploy_predictions.npz';b=npz(canonical)
                for k in a.files:
                    assert np.array_equal(a[k],b[k],equal_nan=True) if a[k].dtype.kind in 'fc' else np.array_equal(a[k],b[k])
                    receipts.append({'dataset':ds,'seed':seed,'array':k,'shape':list(a[k].shape),'dtype':str(a[k].dtype),'sha256':hashlib.sha256(a[k].tobytes()).hexdigest(),'canonical_source':canonical,'exact_equal':True})
            if mem.name.endswith('paired_by_condition.csv'):put(E20+'/compact/paired_by_condition.csv',tf.extractfile(mem).read())
js(E20+'/compact/identity_receipt.json',{'arrays':receipts,'count':len(receipts),'scope':'Exact array identity checked against v24 canonical E19 predictions during compaction; this digest receipt is not a raw-array re-execution.'})
# E22 formal P: reduce exact expression vectors to loss/correlation moments.
sys.path.insert(0,str(R/E22/'src'));import phase_p_prediction_stress as P
from comparator_core import regenerate_mean_baselines
from metric_core import WeightVector,nir
rows=[];checks=[];nirchecks=[]
for ds in P.DATASETS:
    z=npz(f'{E22}/assets/formal_contracts/{ds}.phase_m_contract.npz');ledger=json.loads(read(f'{E22}/assets/formal_contracts/{ds}.split_ledger.json'))
    cn=z['conditions'].astype(str);genes=z['genes'].astype(str);ii={c:i for i,c in enumerate(cn)};full=z['full_means'].astype(float);truths=z['first_means'].astype(float);control=z['control'].astype(float)
    for seed in (1,2,3):
        split=ledger['split_contract'][str(seed)];targets=split['test'];train=split['train'];mu=np.mean([full[ii[v]] for v in train if v!='control'],axis=0)
        a=npz(f'experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo/{ds}/seed{seed}/deploy_predictions.npz');b=npz(f'experiments/18_dual_head_evidence_gate/predictions/formal_combo/{ds}/seed{seed}/deploy_predictions.npz')
        assert a['conditions'].astype(str).tolist()==targets and np.array_equal(a['genes'].astype(str),genes)
        tm,bc=regenerate_mean_baselines(targets,train,cn,full,z['full_counts'].astype(np.int64),control)
        mats={'WitnessCell_v14':a['prediction'].astype(float),'WitnessCell_v13':b['prediction'].astype(float),'trainMean':tm,'baseControl':bc}
        inds=[i for i,c in enumerate(targets) if z['weight_evaluable'][ii[c]] and np.sum(z['weights'][ii[c]])>0]
        for i in inds:
            c=targets[i];t=truths[ii[c]];w=z['weights'][ii[c]].astype(float);d=t-mu;den=float(np.sum(w*(d-np.average(d,weights=w))**2));uden=float(np.sum((d-d.mean())**2));dt=t-control;dt=dt-dt.mean()
            for method,mat in mats.items():
                pred=mat[i];err=(pred-t)**2;dp=pred-control;dp=dp-dp.mean()
                row={'dataset':ds,'seed':seed,'condition':c,'subgroup':str(a['subgroups'][i]),'method':method,'n':len(t),'sse':float(err.sum()),'weight_sum':float(w.sum()),'weighted_sse':float(np.sum(w*err)),'wmse_native':float(np.sum((w/w.sum())*err)),'weighted_r2_sse':float(np.sum(w*((t-mu)-(pred-mu))**2)),'weighted_sst':den,'r2_sse':float(np.sum(((t-mu)-(pred-mu))**2)),'sst':uden,'cross':float(dt@dp),'truth_norm2':float(dt@dt),'pred_norm2':float(dp@dp)};rows.append(row)
                score=P.safe_score(pred,t,mu,control,WeightVector(genes,w));checks.append(dict(dataset=ds,seed=seed,condition=c,method=method,**score))
        nt=np.stack([truths[ii[targets[i]]] for i in inds])
        for method,mat in mats.items():
            ns=nir(mat[inds],nt,[targets[i] for i in inds],[ds]*len(inds));nirchecks.extend(dict(dataset=ds,seed=seed,condition=c,method=method,nir=v) for c,v in ns.items())
    print('E22 loss moments',ds,flush=True)
frame=pd.DataFrame(rows);frame.to_csv(O/E22/'compact_moments.csv',index=False)
original=pd.read_csv(io.BytesIO(read(E22+'/results/formal_e22/P/per_operation_metrics.csv')));generated=pd.DataFrame(checks)
ids=['dataset','seed','condition','method'];a=original.set_index(ids).sort_index();b=generated.set_index(ids).sort_index();cols=['wmse','weighted_r2_deltapert','mse','r2_deltapert','pearson_deltactrl'];diff=np.nanmax(np.abs(a[cols].values-b[cols].values));assert diff<1e-12,diff
n0=pd.read_csv(io.BytesIO(read(E22+'/results/formal_e22/P/per_operation_nir.csv'))).set_index(ids).sort_index();n1=pd.DataFrame(nirchecks).set_index(ids).sort_index();nd=np.max(np.abs(n0.nir-n1.nir));assert nd<1e-12
print('E22 raw score maxdiff',diff,'NIR',nd,flush=True)
# Gate21 weighted finite-row sufficient statistics
z=npz(E22+'/assets/gate21_gene_contract.npz');print('gatefields',z.files,flush=True)
grows=[]
for i in range(len(z['query_ids'])):
    w=z['canonical_weights'][i].astype(float);t=z['truth'][i].astype(float)
    for arm in ['estimated','geometry']:
        p=z[arm+'_prediction'][i].astype(float);grows.append({'query_id':str(z['query_ids'][i]),'arm':arm,'weight_sum':float(w.sum()),'weighted_sse':float(np.sum(w*(p-t)**2)),'wmse_native':float(np.sum((w/w.sum())*(p-t)**2)) if w.sum()>0 else float('nan'),'scoring_evaluable':bool(z['scoring_evaluable'][i])})
pd.DataFrame(grows).to_csv(O/E22/'compact_gate21_moments.csv',index=False)
# E23 loss moments include all five predictor arms, non-target and final predictor full-coordinate diagnostics.
moments=[];maxdiff=0.;frozenchecks=[]
for ds in ['Norman','Wessels']:
  for seed in [1,2,3]:
    base=f'{U}/results/{ds}_seed{seed}';z=npz(base+'/allowed_input.npz');ev=npz(base+'/evaluation_only.npz');bz=npz(base+'/baseline_predictions.npz');md=json.loads(read(base+'/evaluation_metadata.json'))
    genes=z['genes'].astype(str);ci={c:i for i,c in enumerate(z['conditions'].astype(str))};ctrl=z['means'][ci['control']].astype(float);mu=np.mean(z['means'][[ci[c] for c in z['train'].astype(str) if c!='control']].astype(float),axis=0);gi={g:i for i,g in enumerate(genes)}
    frozen=pd.read_csv(io.BytesIO(read(base+'/frozen_query_outputs.csv')));scored=pd.read_csv(io.BytesIO(read(base+'/scored_queries.csv')));baseline=pd.read_csv(io.BytesIO(read(base+'/baseline_losses.csv')))
    mats={'WitnessCell':z['frozen_prediction'].astype(float),**{k:bz[k].astype(float) for k in bz.files if k!='conditions'}}
    official=npz(f'experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo/{ds}/seed{seed}/deploy_predictions.npz');oi={c:i for i,c in enumerate(official['conditions'].astype(str))};pos=[oi[c] for c in ev['conditions'].astype(str)]
    assert np.array_equal(z['frozen_prediction'],official['prediction'][pos]);frozenchecks.append({'dataset':ds,'seed':seed,'rows':len(frozen),'official_prediction_exact_equal':True})
    for i,c in enumerate(ev['conditions'].astype(str)):
        assert frozen.iloc[i].condition==c;mask=np.ones(len(genes),bool)
        for g in c.split('+'):mask[gi[g]]=False
        for method,mat in mats.items():
          for scope in (['non_target','including_targets'] if method=='WitnessCell' else ['non_target']):
            mm=mask if scope=='non_target' else np.ones(len(genes),bool);p=mat[i];t=ev['truth'][i].astype(float);w=ev['weights'][i].astype(float);ww=w[mm];err=(p[mm]-t[mm])**2;de=t[mm]-ctrl[mm];pe=p[mm]-ctrl[mm];d=t[mm]-mu[mm];dp=p[mm]-mu[mm];top=sorted(np.flatnonzero(mm),key=lambda j:(-w[j],str(genes[j])))[:100]
            rr={'dataset':ds,'seed':seed,'pair':frozen.iloc[i]['pair'],'condition':c,'method':method,'scope':scope,'n':int(mm.sum()),'sse':float(err.sum()),'weight_sum':float(ww.sum()),'weighted_sse':float(ww@err),'effect_sse':float(de@de),'effect_native':float(np.mean(de**2)),'cos_dot':float(de@pe),'cos_den':float(np.linalg.norm(de)*np.linalg.norm(pe)),'top_n':len(top),'top_sse':float(np.sum((p[top]-t[top])**2)),'weighted_r2_sse':float(ww@((d-dp)**2)),'weighted_sst':float(ww@((d-np.average(d,weights=ww))**2)) if ww.sum()>0 else 0.,'relative_epsilon':md['relative_loss_epsilon']};moments.append(rr)
            vals={'loss_mse':rr['sse']/rr['n'],'loss_wmse':rr['weighted_sse']/rr['weight_sum'] if rr['weight_sum'] else float('nan'),'true_effect':rr['effect_native'],'loss_relative':(rr['sse']/rr['n'])/max(rr['effect_native'],rr['relative_epsilon']),'loss_cosine':1-float(np.clip(rr['cos_dot']/rr['cos_den'],-1,1)) if rr['cos_den'] else float('nan'),'loss_top100':rr['top_sse']/rr['top_n'],'weighted_delta_r2':1-rr['weighted_r2_sse']/rr['weighted_sst'] if rr['weighted_sst'] else float('nan')}
            ref=scored.iloc[i] if method=='WitnessCell' else baseline[(baseline.condition==c)&(baseline.method==method)].iloc[0]
            for k,v in vals.items():
                kk=k+('_including_targets' if scope=='including_targets' else '');q=float(ref[kk]);delta=abs(v-q) if np.isfinite(v) and np.isfinite(q) else 0.;maxdiff=max(maxdiff,delta);assert delta<1e-11,(ds,seed,c,method,k,delta)
    print('E23 moments',ds,seed,'maxdiff',maxdiff,flush=True)
pd.DataFrame(moments).to_csv(O/E23/'compact_loss_moments.csv',index=False)
js('BUILD_VERIFICATION.json',{'E20_exact_prediction_array_identities':len(receipts),'E22_original_vs_rescored_max_absolute_delta':float(diff),'E22_NIR_max_absolute_delta':float(nd),'E22_prediction_rows':len(frame),'E23_original_vs_rescored_max_absolute_delta':float(maxdiff),'E23_loss_moment_rows':len(moments),'E23_official_prediction_identity':frozenchecks,'scope':'Computed from decoded original v24 expression/prediction/weight arrays, then verified against frozen per-query scores; compact replay recomputes scores from retained moments and all resampling from rows.'})
js('SOURCE_INPUTS_SHA256.json',origins)
print('BUILD COMPLETE',sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),flush=True)
