#!/usr/bin/env python3
"""Produce exact compact score inputs and array certificates from the full v24 store."""
from pathlib import Path
import argparse, gzip, hashlib, importlib.util, io, json, shutil
import numpy as np

def main():
 p=argparse.ArgumentParser();p.add_argument('--source-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.source_root.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
 idx=json.loads((r/'reproduction_store/INDEX.json').read_text());spec=importlib.util.spec_from_file_location('legacy_store',r/'tools/restore_data_v24.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);store=mod.Store(r,idx)
 sources=[]
 def read(rel):
  path=r/rel
  b=path.read_bytes() if path.exists() else b''.join(store.parts(idx['files'][rel]))
  if rel in idx['files']:
   x=idx['files'][rel];assert len(b)==x['bytes'] and hashlib.sha256(b).hexdigest()==x['sha256']
  return b
 def put(rel,dest):
  b=read(rel);q=out/dest;q.parent.mkdir(parents=True,exist_ok=True)
  q.write_bytes(gzip.compress(b,compresslevel=9,mtime=0) if dest.endswith('.gz') else b)
  sources.append({'source_path':rel,'source_bytes':len(b),'source_sha256':hashlib.sha256(b).hexdigest(),'compact_path':dest,'transformation':'gzip_lossless' if dest.endswith('.gz') else 'unchanged'})
 E={17:'17_endpoint_sota_reproduction',18:'18_dual_head_evidence_gate',19:'19_v14_incremental_amplitude_gate',24:'24_v14_official_six_metric_release'}
 ablations={'endpoint':(17,'formal_paired_three_seed'),'dual_head':(18,'joint_k10_formal_12_score'),'released_vs_prior':(19,'24_release_reproduction_v14_vs_v13'),'released_vs_backbone':(19,'22_formal_all12_vs_v1')}
 for label,(num,stem) in ablations.items():
  for name in ['paired_by_condition.csv','paired_summary.csv','verdict.json']:
   put(f'experiments/{E[num]}/results/{stem}/{name}',f'data/ablations/{label}/{name}'+('.gz' if name.endswith('.csv') else ''))
 score_roots={'endpoint':f'experiments/{E[17]}/results/formal_official','dual_head':f'experiments/{E[18]}/results/joint_k10_official','released':f'experiments/{E[24]}/payload/results/26_v14_official_six_metric'}
 for label,stem in score_roots.items():
  for name in ['published_plus_witness_per_operation.csv','formal_top100_leaderboard.csv','formal_equal_weight_leaderboard.csv','formal_both_panels_leaderboard.csv','dataset_equal_weight_top100_sensitivity.csv','dataset_panel_leaderboard.csv','official_top100_per_metric_leaderboard.csv','paired_cluster_bootstrap.json','paired_cluster_bootstrap.csv','formal_verdict.json']:
   put(stem+'/aggregate/'+name,f'data/official/{label}/{name}'+('.gz' if name.endswith('.csv') else ''))
 for dataset in ['Norman','Replogle_exp6','Schmidt','Wessels']:
  for seed in [1,2,3]:
   put(score_roots['released']+f'/{dataset}/raw_seed{seed}.csv',f'data/official/released/raw/{dataset}_seed{seed}.csv.gz')
 put(score_roots['released']+'/assembly_manifest.json','data/official/released/assembly_manifest.json')
 put(score_roots['released']+'/exact_fallback_audit.json','data/official/released/fallback_source_audit.json')
 # Keep the released model's complete local import closure exactly, not all development runners.
 for name in ['identity_head.py','official_split.py','dual_head.py','go_identity_head.py','witnesscell_v1_core.py','witnesscell_v4_core.py','multihead_v14.py','probe_multhead_training.py','run_v14_from_cache.py']:
  put(f'experiments/{E[19]}/src/{name}','source/v14/'+name)
 put(f'experiments/{E[19]}/ALGORITHM_CONTRACT.md','source/ALGORITHM_CONTRACT.md')
 put('experiments/15_scperturbench_sota/module/aggregate_official_full.py','source/aggregate_official_full.py')
 for num,label in [(17,'v12_endpoint'),(18,'v13_dual')]:
  for f in sorted((r/'experiments'/E[num]/'src').glob('*.py')):
   put(str(f.relative_to(r)),f'source/{label}/{f.name}')
 put(f'experiments/{E[19]}/src/compare_v14_v13.py','source/v14/compare_v14_v13.py')
 # Directly inspect every archived deploy array. Compact certificates explicitly
 # record inspections; their source matrices are intentionally not rebundled.
 array_records=[];fallback=[]
 for dataset in ['Norman','Replogle_exp6','Schmidt','Wessels']:
  for seed in [1,2,3]:
   arrays={};rawhash={};meta={}
   for num in [17,18,19]:
    rel=f'experiments/{E[num]}/predictions/formal_combo/{dataset}/seed{seed}/deploy_predictions.npz';b=read(rel);rawhash[num]=hashlib.sha256(b).hexdigest();z=np.load(io.BytesIO(b),allow_pickle=False);arrays[num]={k:z[k] for k in z.files};meta[num]=json.loads(read(rel.replace('deploy_predictions.npz','manifest.json')))
    assert not any(k in z.files for k in ['truth','truth_variance']);assert all(z[k].dtype.kind!='O' for k in z.files)
    array_records.append({'experiment':num,'dataset':dataset,'seed':seed,'source_path':rel,'source_sha256':rawhash[num],'source_bytes':len(b),'arrays':[{'key':k,'shape':list(z[k].shape),'dtype':str(z[k].dtype),'array_sha256':hashlib.sha256(z[k].tobytes()).hexdigest()} for k in z.files],'no_truth_or_object_arrays':True,'identity_head':meta[num]['identity_head']})
   current,prior=arrays[19],arrays[18];inactive=not meta[19]['identity_head']['correction_active'];delta=current['prediction'].astype(float)-prior['prediction'].astype(float)
   if inactive:assert np.array_equal(current['prediction'],prior['prediction'])
   fallback.append({'dataset':dataset,'seed':seed,'active':not inactive,'max_abs_prediction_delta':float(np.max(np.abs(delta))),'prediction_delta_sse':float(np.sum(delta*delta)),'conditions_equal':bool(np.array_equal(current['conditions'],prior['conditions'])),'inactive_exact_fallback':bool(not inactive or np.array_equal(current['prediction'],prior['prediction'])),'v14_sha256':rawhash[19],'prior_sha256':rawhash[18]})
 (out/'data/ARRAY_CERTIFICATES.json').write_text(json.dumps({'scope':'Producer-time full-array inspection; default compact replay does not re-load original gene-level matrices','records':array_records,'v14_vs_prior':fallback},indent=2)+'\n')
 (out/'SOURCE_MAP.json').write_text(json.dumps({'source_release':'grn_repo_v24','source_zip_sha256':'dae2d0f616bdbd025a0452c63e317848d504cf5832c2fd12f0e471f0c9895330','files':sources,'large_array_certificates':'data/ARRAY_CERTIFICATES.json','preservation':'All compact CSV payloads losslessly decode to original source bytes; no quantization or rounded replacement values'},indent=2)+'\n')
 print(json.dumps({'status':'PASS','source_files':len(sources),'deploy_arrays_inspected':len(array_records),'bytes':sum(p.stat().st_size for p in out.rglob('*') if p.is_file())},indent=2))

if __name__=='__main__':main()
