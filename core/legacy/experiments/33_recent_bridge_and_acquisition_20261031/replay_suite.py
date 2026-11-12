"""Reproduce the nine packaged rounds in a new disposable directory.

Inputs are restored ORIGINAL files, not the compressed object store itself.
Scientific sources and frozen historical outcomes are never edited in place.
Default scope is frozen-weight prediction/score replay; smoke/admission use
their original deterministic end-to-end replay. Inner-cache regeneration is
an optional separate first step, reported explicitly.
"""
from pathlib import Path
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import argparse,json,shutil,subprocess,sys
from compact_bridge_replay import numeric_equal,compare_csv

SPEC=[
 ('01','WitnessCell_NM_NCS_Bridge_20261031','bridge_experiment_20261031'),
 ('02','bridge_retrospective_20261031','bridge_retrospective_20261031'),
 ('03','WitnessCell_Targeted_Bridge_Calibration_20261031','target_bridge_calibration_20261031'),
 ('04','WitnessCell_Formal_Hidden_Mechanism_Bridge_20261031','formal_bridge_20261031'),
 ('05','WitnessCell_OneStep_Minimal_Smoke_20261031','one_step_smoke_20261031'),
 ('06','WitnessCell_Training_Admission_Calibration_20261031','admission_calibration_20261031'),
 ('07','direction_alignment_20261031','direction_alignment_20261031'),
 ('08','target_aligned_run_20261031','target_aligned_run_20261031'),
 ('09','target_aligned_outer_20261031','target_aligned_outer_20261031')]


def main():
 p=argparse.ArgumentParser();p.add_argument('--restored',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
 p.add_argument('--stages',default='01,02,03,04,05,06,07,08,09');p.add_argument('--rebuild-inner-cache',type=Path)
 a=p.parse_args();src=a.restored.resolve();out=a.out.resolve();selected=set(a.stages.split(','))
 assert not out.exists(),'Output must be new';out.mkdir(parents=True);(out/'logs').mkdir()
 dirs={};commands=[];checks=[]
 for id_,name,alias in SPEC:
  direct=src/name
  matches=[direct] if direct.is_dir() else [q for q in src.rglob(name) if q.is_dir()]
  assert len(matches)==1,(name,'missing or ambiguous',matches)
  dst=out/alias;shutil.copytree(matches[0],dst);dirs[id_]=dst
 def run(root,script,*args):
  log=out/'logs'/f'{len(commands):02d}_{Path(script).stem}.log'
  with log.open('w') as f:r=subprocess.run([sys.executable,str(root/script),*map(str,args)],cwd=root,stdout=f,stderr=subprocess.STDOUT)
  commands.append({'script':str(root/script),'args':list(map(str,args)),'returncode':r.returncode,'log':str(log.relative_to(out))})
  print(json.dumps({'finished':Path(script).name,'returncode':r.returncode}),flush=True)
  assert r.returncode==0,log.read_text()[-8000:]
 def compare_scores(before,after):
  matched=[]
  for old in sorted(before.iterdir()):
   new=after/old.name
   if old.name=='CACHE_EXECUTION.json':continue
   if old.suffix=='.csv':assert new.exists();compare_csv(old,new);matched.append(old.name)
   elif old.suffix=='.json':assert new.exists();numeric_equal(json.loads(old.read_text()),json.loads(new.read_text()),old.name);matched.append(old.name)
  return matched
 if a.rebuild_inner_cache:
  script=a.rebuild_inner_cache.resolve()
  run(script.parent,script.name,'--root',dirs['01'],'--restore')
 if selected&{'01','02'}:
  helper=Path(__file__).resolve().parent/'compact_bridge_replay.py'
  run(helper.parent,helper.name,'--bridge',dirs['01'],'--retrospective',dirs['02'],'--out',out/'bridge_replay')
  checks.append({'stages':['01','02'],'result':'PASS','evidence':'bridge_replay/COMPACT_BRIDGE_REPLAY_RESULT.json'})
 if '03' in selected:
  original=dirs['03'];old=original/'calibration/scores';root=out/'stage03_rescore'
  shutil.copytree(original,root,ignore=shutil.ignore_patterns('scores','scores_expected'))
  assert not (root/'calibration/scores').exists()
  run(root,'calibration/independent_score.py')
  checks.append({'stage':'03','result':'PASS','compared':compare_scores(old,root/'calibration/scores')})
 if '04' in selected:
  original=dirs['04'];old=original/'scores';root=out/'stage04_rescore'
  shutil.copytree(original,root,ignore=shutil.ignore_patterns('scores','scores_expected'))
  assert not (root/'scores').exists()
  run(root,'score_cached.py');run(root,'verify_frozen_statistics.py')
  checks.append({'stage':'04','result':'PASS','compared':compare_scores(old,root/'scores')})
 if '05' in selected:
  root=dirs['05'];run(root,'independent_audit.py');run(root,'reproduce.py','--out',out/'fresh_smoke')
  checks.append({'stage':'05','result':'PASS','scope':'original complete deterministic replay + independent audit'})
 if '06' in selected:
  root=dirs['06'];run(root,'audit.py');run(root,'reproduce.py','--out',out/'fresh_admission')
  checks.append({'stage':'06','result':'PASS','scope':'original complete deterministic replay + independent audit'})
 if '07' in selected:
  root=dirs['07'];names=('alignment_diagnostic.json','summary.json','reference_probe.json')
  old={n:json.loads((root/n).read_text()) for n in names}
  run(root,'diagnose_alignment.py');run(root,'probe_reference.py')
  for n,v in old.items():numeric_equal(v,json.loads((root/n).read_text()),n)
  checks.append({'stage':'07','result':'PASS','compared':list(names),'ignored_metadata':['opened_npz','npz_opened']})
 if '08' in selected:
  root=dirs['08'];run(root,'independentcheck.py');run(root,'reproduce.py','--out',out/'fresh_aligned_calibration')
  checks.append({'stage':'08','result':'PASS','scope':'fixed local fit/selection/score replay; original failed gate preserved'})
 if '09' in selected:
  root=dirs['09'];run(root,'independent_outer_audit.py');run(root,'reproduce_outer.py',out/'fresh_outer')
  checks.append({'stage':'09','result':'PASS','scope':'30 prediction arrays exact and five JSON exact'})
 result={'status':'PASS','requested_stages':sorted(selected),'checks':checks,'commands':commands,
         'scope':'reproduction of historical frozen outcomes; no new tuning, labels, targets or scientific result',
         'inner_cache_rebuild_requested':bool(a.rebuild_inner_cache)}
 (out/'REPLAY_SUITE_RESULT.json').write_text(json.dumps(result,indent=2))
 print(json.dumps({'status':'PASS','checks':checks},indent=2))


if __name__=='__main__':main()
