#!/usr/bin/env python3
"""Current paper v36 entry point. Default tasks perform no biological model training."""
from pathlib import Path
import argparse,subprocess,sys,json,os,time,signal
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--task',default='all',help='all, or comma-separated task IDs from REPRODUCTION_INDEX.json');ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out.resolve()
if out.exists() or out.is_relative_to(R):ap.error('Use a fresh output directory outside this release.')
index=json.loads((R/'REPRODUCTION_INDEX.json').read_text());tasks=index['tasks'];selected=list(tasks) if a.task=='all' else a.task.split(',')
if any(k not in tasks for k in selected):ap.error('Unknown task. Valid: '+', '.join(tasks))
out.mkdir(parents=True);env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONOPTIMIZE='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
report=dict(release='28.0.0',paper='v36',status='RUNNING',tasks=[],base_model_training_enabled=False)
for task in selected:
 spec=tasks[task];cmd=[sys.executable,'-B',str(R/spec['script'])]+[x.replace('{out}',str(out/task)).replace('{root}',str(R)) for x in spec['args']]
 print('Running '+task,flush=True);start=time.monotonic()
 with (out/(task+'.log')).open('w') as log:
  process=subprocess.Popen(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  try:code=process.wait(timeout=spec.get('timeout_seconds',900));status='PASS' if code==0 else 'FAIL'
  except subprocess.TimeoutExpired:
   os.killpg(process.pid,signal.SIGTERM)
   try:process.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
   code=process.returncode;status='TIMEOUT'
 row=dict(task=task,status=status,returncode=code,seconds=time.monotonic()-start,scope=spec['scope'],command=cmd);report['tasks'].append(row);(out/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n');print(task+': '+status,flush=True)
report['status']='PASS' if all(x['status']=='PASS' for x in report['tasks']) else 'FAIL';(out/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if report['status']!='PASS':sys.exit(1)
