#!/usr/bin/env python3
"""Reproduce frozen E20--E23 metrics and uncertainty from compact intermediates.
No predictor training, model selection or raw-cell processing occurs here.
"""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json,os,sys
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent
E20=R/'experiments/20_osf_preregistered_replication';E21=R/'experiments/21_frozen_selective_prediction';E22=R/'experiments/22_metric_calibration_stress';E23=R/'experiments/23_joint_prediction_risk_step1';U=E23/'witnesscell_step1_run'
report={};maxerr=0.
resampling_receipts=json.loads((R/'ORIGINAL_RESAMPLING_ARRAYS_SHA256.json').read_text())
def check_array(rel,key,array):
 ref=resampling_receipts[rel]['arrays'][key];assert list(array.shape)==ref['shape'];assert str(array.dtype)==ref['dtype']
 if hashlib.sha256(array.tobytes()).hexdigest()!=ref['sha256']:
  original=R/rel;assert original.exists(),(rel,key,'resampling byte mismatch without reference array')
  assert hashlib.sha256(original.read_bytes()).hexdigest()==resampling_receipts[rel]['npz_sha256']
  with np.load(original,allow_pickle=False) as z:
   assert np.array_equal(np.isnan(array),np.isnan(z[key]));delta=float(np.nanmax(np.abs(array-z[key])));assert delta<1e-12,(rel,key,delta)
   global maxerr;maxerr=max(maxerr,delta)
def load(p):return json.loads(Path(p).read_text())
def imp(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def eq(a,b,path=''):
 global maxerr
 if isinstance(a,dict):
  assert set(a)==set(b),(path,set(a)^set(b))
  for k in a:eq(a[k],b[k],path+'/'+k)
 elif isinstance(a,(list,tuple)):
  assert len(a)==len(b),path
  for i,(x,y) in enumerate(zip(a,b)):eq(x,y,path+'/'+str(i))
 elif isinstance(a,(int,float,np.number)) and not isinstance(a,(bool,np.bool_)):
  if a is None or b is None:assert a is b,path
  elif not np.isfinite(float(a)) or not np.isfinite(float(b)):assert (not np.isfinite(float(a))) and (not np.isfinite(float(b))),path
  else:
   d=abs(float(a)-float(b));maxerr=max(maxerr,d);assert d<=1e-10,(path,a,b,d)
 else:assert a==b,(path,a,b)
def numcheck(a,b,cols,keys):
 a=a.set_index(keys).sort_index();b=b.set_index(keys).sort_index();assert a.index.equals(b.index)
 for c in cols:
  x=a[c].to_numpy(float);y=b[c].to_numpy(float);assert np.array_equal(np.isnan(x),np.isnan(y)),c
  good=np.isfinite(x)&np.isfinite(y);d=float(np.max(np.abs(x[good]-y[good]))) if good.any() else 0.;assert d<1e-10,(c,d)
  global maxerr;maxerr=max(maxerr,d)
def E20_replay():
 f=pd.read_csv(E20/'compact/paired_by_condition.csv');ref=pd.read_csv(E20/'evidence/paired_summary.csv');rows=[]
 for row in ref.to_dict('records'):
  scope=row['scope'];x=f if scope=='all' else f[(f.dataset==scope)|(f.subgroup==scope)]
  out={'scope':scope,'condition_seed_units':len(x)}
  for k in ['all_gene_mse_v13','all_gene_mse_v14','top100_mse_v13','top100_mse_v14','top100_pcc_v13','top100_pcc_v14']:out[k]=float(x[k].mean())
  for prefix in ['all_gene','top100']:out['relative_'+prefix+'_mse_gain']=1-out[prefix+'_mse_v14']/out[prefix+'_mse_v13']
  out['top100_pcc_delta']=out['top100_pcc_v14']-out['top100_pcc_v13'];eq(out,row,'E20/'+scope);rows.append(out)
 receipt=load(E20/'compact/identity_receipt.json');assert receipt['count']==156 and all(x['exact_equal'] for x in receipt['arrays'])
 return {'aggregate_rows':rows,'source_identity_receipt_arrays':156,'status':'PASS_COMPACT_METRIC_REPLAY'}
def selective(frame,loss,threshold,ref):
 import selective_core as s
 primary=s.evaluate_threshold(frame,'estimated_witness_risk',loss,threshold);eq(primary,{k:ref['primary'][k] for k in primary},'selective/primary')
 boot=s.cluster_bootstrap_ratio(frame,'estimated_witness_risk',loss,threshold,20000,20260811);perm=s.within_seed_permutation(frame,'estimated_witness_risk',loss,threshold,20000,20260812)
 rel='experiments/21_frozen_selective_prediction/results/formal_reveal/inference_arrays.npz' if loss=='estimated_realized_mse' else 'experiments/22_metric_calibration_stress/results/formal_e22/D/gate21/gate21_inference_arrays.npz'
 check_array(rel,'bootstrap_accepted_over_all_ratio' if loss=='estimated_realized_mse' else 'witness_bootstrap_ratio',boot)
 check_array(rel,'random_selection_accepted_mse' if loss=='estimated_realized_mse' else 'witness_random_selection_wmse',perm)
 ci=np.quantile(boot,[.025,.975]).tolist();eq(ci,ref['bootstrap_accepted_over_all_ci95'],'selective/CI')
 p=float((1+np.sum(perm<=primary['accepted_mse']))/(len(perm)+1));eq(p,ref['within_seed_random_selection_p_one_sided'],'selective/permutation')
 return {'primary':primary,'bootstrap_ci95':ci,'permutation_p':p,'replicates':20000,'bootstrap_seed':20260811,'permutation_seed':20260812}
def E21_replay():
 import selective_core as s
 assert s.sha256(E21/'FROZEN_PROTOCOL.json')==(E21/'FROZEN_PROTOCOL.sha256').read_text().split()[0]
 risk=pd.read_csv(E21/'assets/formal_risk_only.csv');out=pd.read_csv(E21/'sealed/formal_outcomes.csv');f=risk.merge(out,on='query_id',validate='one_to_one');t=pd.read_csv(E21/'results/frozen_calibration/witness_thresholds.csv');threshold=float(t.loc[np.isclose(t.target_coverage,.5),'risk_threshold'].iloc[0]);ref=load(E21/'results/formal_reveal/FORMAL_VERDICT.json')
 result=selective(f,'estimated_realized_mse',threshold,ref);result['scientific_verdict']=ref['status'];return result
def E22_replay():
 import phase_p_prediction_stress as p
 f=pd.read_csv(E22/'compact_moments.csv',float_precision='round_trip');r=f[['dataset','seed','condition','subgroup','method']].copy()
 r['wmse']=f.wmse_native;r['mse']=f.sse/f.n;r['weighted_r2_deltapert']=1-f.weighted_r2_sse/f.weighted_sst;r['r2_deltapert']=1-f.r2_sse/f.sst;r['pearson_deltactrl']=f['cross']/np.sqrt(f.truth_norm2*f.pred_norm2)
 original=pd.read_csv(E22/'results/formal_e22/P/per_operation_metrics.csv');numcheck(r,original,['wmse','mse','weighted_r2_deltapert','r2_deltapert','pearson_deltactrl'],['dataset','seed','condition','method'])
 ref=load(E22/'results/formal_e22/P/prediction_verdict.json');result={}
 # Use archival decimal rows after verifying reconstruction, preserving original CSV parsing and bootstrap arithmetic.
 for base in ['WitnessCell_v13','trainMean','baseControl']:
  point=p.comparison_point(original,base);point.pop('splits');boot,br=p.cluster_bootstrap(original,base);point['wmse_ratio_ci95']=np.quantile(boot,[.025,.975]).tolist();point['weighted_r2_difference_ci95']=np.quantile(br,[.025,.975]).tolist()
  check_array('experiments/22_metric_calibration_stress/results/formal_e22/P/cluster_bootstrap_arrays.npz',base+'__wmse_ratio',boot);check_array('experiments/22_metric_calibration_stress/results/formal_e22/P/cluster_bootstrap_arrays.npz',base+'__weighted_r2_difference',br)
  for k,v in point.items():eq(v,ref['comparisons'][base][k],'E22/P/'+base+'/'+k)
  result[base]=point
 combo=original[original.condition.str.contains('+',regex=False)];pt=p.comparison_point(combo,'trainMean');pt.pop('splits');boot,br=p.cluster_bootstrap(combo,'trainMean');pt['wmse_ratio_ci95']=np.quantile(boot,[.025,.975]).tolist();pt['weighted_r2_difference_ci95']=np.quantile(br,[.025,.975]).tolist();result['paper_v22_combinations_only']=pt
 assert abs(pt['dataset_equal_wmse_ratio']-.7513124080560089)<1e-12
 # Phase M continuous DRF and NIR summaries, retaining negative/undefined outputs.
 dr=[];nr=[];elig=[]
 for ds in p.DATASETS:
  v=pd.read_csv(E22/f'results/formal_e22/M/{ds}.metric_validity.csv');n=pd.read_csv(E22/f'results/formal_e22/M/{ds}.nir_controls.csv')
  for col in v.columns:
   if not col.startswith('drf__'):continue
   _,neg,pos,metric=col.split('__',3);x=v.loc[v.evaluable.astype(bool),col].to_numpy(float);x=x[np.isfinite(x)];dr.append(dict(dataset=ds,metric=metric,negative_control=neg,positive_control=pos,valid_n=len(x),mean_drf=float(np.mean(x)) if len(x) else np.nan,median_drf=float(np.median(x)) if len(x) else np.nan,fraction_drf_le_zero=float(np.mean(x<=0)) if len(x) else np.nan))
  for (seed,control),local in n.groupby(['seed','control']):nr.append(dict(dataset=ds,seed=seed,control=control,conditions=len(local),mean_nir=local.nir.mean(),median_nir=local.nir.median()))
 numcheck(pd.DataFrame(dr),pd.read_csv(E22/'results/formal_e22/M/drf_summary.csv'),['valid_n','mean_drf','median_drf','fraction_drf_le_zero'],['dataset','metric','negative_control','positive_control'])
 numcheck(pd.DataFrame(nr),pd.read_csv(E22/'results/formal_e22/M/nir_control_summary.csv'),['conditions','mean_nir','median_nir'],['dataset','seed','control'])
 g=pd.read_csv(E22/'compact_gate21_moments.csv',float_precision='round_trip');gf=pd.read_csv(E22/'results/formal_e22/D/gate21/gate21_wmse_rows.csv')
 for arm in ['estimated','geometry']:
  a=g[(g.arm==arm)&g.scoring_evaluable].set_index('query_id');b=gf[gf.scoring_evaluable].set_index('query_id');assert a.index.equals(b.index);eq(a.wmse_native.tolist(),b[arm+'_wmse'].tolist(),'E22/D/moments')
 # A1 index reset is necessary before seed-local permutation; A2 preserves complete-case risk alignment.
 finite=gf.loc[gf.scoring_evaluable].reset_index(drop=True);dref=load(E22/'results/formal_e22/D/gate21/decision_verdict.json');result['gate21_complete_case']=selective(finite,'estimated_wmse',.0923227147328771,dref);result['gate21_full_status']=dref['GATE21_WMSE'];result['gate21_complete_case_status']=dref['GATE21_WMSE_COMPLETE_CASE'];result['bootstrap_replicates']=20000;result['bootstrap_seed']=20260822
 return result
def E23_replay():
 m=pd.read_csv(E23/'compact_loss_moments.csv',float_precision='round_trip');m['loss_mse']=m.sse/m.n;m['loss_wmse']=m.weighted_sse/m.weight_sum;m['true_effect']=m.effect_native;m['loss_relative']=m.loss_mse/np.maximum(m.effect_native,m.relative_epsilon);m['loss_cosine']=1-np.clip(m.cos_dot/m.cos_den,-1,1);m['loss_top100']=m.top_sse/m.top_n;m['weighted_delta_r2']=1-m.weighted_r2_sse/m.weighted_sst
 metrics=['loss_mse','loss_wmse','true_effect','loss_relative','loss_cosine','loss_top100','weighted_delta_r2']
 for (ds,seed),part in m.groupby(['dataset','seed']):
  d=U/f'results/{ds}_seed{seed}';baseline=pd.read_csv(d/'baseline_losses.csv');numcheck(part[part.scope=='non_target'],baseline,metrics,['dataset','seed','condition','method'])
  f=pd.read_csv(d/'frozen_query_outputs.csv');s=pd.read_csv(d/'scored_queries.csv');assert (f.condition==s.condition).all()
  for field in ['risk_S','risk_B','effect_score','variance_score']:eq(f[field].tolist(),s[field].tolist(),'E23/frozen/'+field)
  full=part[(part.scope=='including_targets')&(part.method=='WitnessCell')].copy();full=full.rename(columns={x:x+'_including_targets' for x in metrics});numcheck(full,s,[x+'_including_targets' for x in metrics],['condition'])
 ev=imp('step1evaluate',U/'evaluate.py');base=imp('step1baselines',U/'summarize_baselines.py');batches,_=ev.load_batches(U/'results');bb,_=base.load_batches(U/'results');ref=load(U/'evaluation/evaluation.json');bref=load(U/'evaluation/baseline_summary.json');result={};boot_hashes={}
 for ds in ['Norman','Wessels']:
  r,arrays=ev.analyze_dataset(ds,batches,5000,1729);eq(ev.json_safe(r),ref['datasets'][ds],'E23/'+ds)
  for k,a in arrays.items():
   check_array('experiments/23_joint_prediction_risk_step1/witnesscell_step1_run/evaluation/bootstrap_ratios.npz',k,a);boot_hashes[k]=hashlib.sha256(a.tobytes()).hexdigest()
  b,_,_=base.analyze_dataset(ds,bb,5000,1729);eq(base.json_safe(b),bref['datasets'][ds],'E23/baselines/'+ds)
  result[ds]={'risk':r,'baselines':b}
 # Reconstruct all frozen acceptance flags, never rerank with outcomes.
 flags=[]
 for b in batches.values():
  for c in ev.COVERAGES:
   for scorer in ev.SCORERS:
    a=ev.acceptance(b,scorer,c);flags.extend((b['dataset'],b['seed'],pair,c,scorer,bool(flag)) for pair,flag in zip(b['pairs'],a))
 archived=pd.read_csv(U/'evaluation/selection_manifest.csv');saved={(r.dataset,int(r.seed),r.pair,float(r.coverage_requested),r.scorer,bool(r.accepted)) for r in archived.itertuples()};assert set(flags)==saved
 frozen=set()
 for path in (U/'results').glob('*/frozen_selections.csv'):
  for r in pd.read_csv(path).itertuples():frozen.add((r.dataset,int(r.seed),r.pair,float(r.coverage),r.method,bool(r.accepted)))
 assert set(flags)==frozen
 result['selection_flags_recomputed']=len(flags);assert len(flags)==6224
 result['bootstrap_arrays_sha256']=boot_hashes;result['risk_gate']=ref['gate'];result['prediction_competitiveness_gate']=bref['prediction_competitiveness_gate'];return ev.json_safe(result)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
 sys.path[:0]=[str(E21),str(E22/'src')]
 for name,fn in [('E20',E20_replay),('E21',E21_replay),('E22',E22_replay),('E23',E23_replay)]:
  print('Replaying',name,flush=True);report[name]=fn()
 report['max_absolute_reproduction_delta']=maxerr;report['status']='PASS_COMPACT_REPLAY';report['scope']='Frozen prediction-loss moments, row metrics, selection and original bootstrap protocols. No raw-cell processing or model retraining.'
 (args.out/'REPLAY_REPORT.json').write_text(json.dumps(report,indent=2,allow_nan=True)+'\n');print(json.dumps({'status':report['status'],'max_absolute_reproduction_delta':maxerr,'report':str(args.out/'REPLAY_REPORT.json')},indent=2))
if __name__=='__main__':main()
