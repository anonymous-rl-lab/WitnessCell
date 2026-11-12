#!/usr/bin/env python3
"""Unified offline v41 replay with per-task checkpoints and explicit refit routes."""
from pathlib import Path
import argparse,json,subprocess,sys,os,time,concurrent.futures,hashlib
R=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);ap.add_argument('--tasks',nargs='+',default=['all']);ap.add_argument('--jobs',type=int,default=2);ap.add_argument('--resume',action='store_true');ap.add_argument('--compatible',action='store_true');ap.add_argument('--list',action='store_true');a=ap.parse_args()
meta=json.loads((R/'REPRODUCTION_INDEX.json').read_text());tasks=meta['tasks']
if a.list:
 for key,task in tasks.items():print(key,'[default]' if task.get('default',True) else '[optional]',task['scope'])
 raise SystemExit(0)
if not a.out:ap.error('--out is required unless --list is used')
out=a.out.resolve()
if (out==R or R in out.parents or out in R.parents) or (out.exists() and not a.resume):ap.error('Choose fresh output outside the release, or --resume an existing run.')
if not 1<=a.jobs<=8:ap.error('--jobs must be 1 through 8')
names=[k for k,v in tasks.items() if v.get('default',True)] if a.tasks==['all'] else a.tasks
if len(names)!=len(set(names)) or not set(names)<=set(tasks):ap.error('Unknown or repeated task; use --list')
out.mkdir(parents=True,exist_ok=True);env=os.environ.copy()
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
env['PYTHONDONTWRITEBYTECODE']='1';env['MPLBACKEND']='Agg'
release_digest=hashlib.sha256((R/'MANIFEST.json').read_bytes()).hexdigest()
subprocess.run([sys.executable,'-B',str(R/'tools/verify.py'),'--out',str(out/'RELEASE_INTEGRITY.json')],check=True,env=env)
cmd=[sys.executable,'-B',str(R/'tools/prepare.py'),'--out',str(out/'preparation')]
if not any(tasks[name].get('requires_han_cache',False) for name in names):cmd+=['--aliases-only']
if a.compatible:cmd+=['--compatible']
subprocess.run(cmd,check=True,env=env)
def dump(path,value):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)
def run(name):
 task=tasks[name];dest=out/name;stat=out/(name+'.status.json');fingerprint=hashlib.sha256((release_digest+json.dumps(task,sort_keys=True)).encode()).hexdigest()
 if a.resume and stat.exists():
  old=json.loads(stat.read_text())
  if old.get('status')=='PASS' and old.get('task_fingerprint')==fingerprint:print(name,'SKIP: already passed',flush=True);return old
 if dest.exists():dest.rename(out/(name+'.previous_'+str(time.time_ns())))
 cmd=[sys.executable,'-B',str(R/task['script']),*[x.replace('{out}',str(dest)).replace('{root}',str(R)) for x in task['args']]]
 report={'task':name,'status':'RUNNING','task_fingerprint':fingerprint,'scope':task['scope'],'command':cmd};dump(stat,report);start=time.monotonic()
 try:
  with (out/(name+'.log')).open('w') as log:p=subprocess.run(cmd,cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=task.get('timeout_seconds',1800))
  report.update(status='PASS' if p.returncode==0 else 'FAIL',returncode=p.returncode)
 except Exception as e:report.update(status='FAIL',error=str(e))
 report['seconds']=time.monotonic()-start;dump(stat,report);print(name,report['status'],round(report['seconds'],2),flush=True);return report
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:results=list(pool.map(run,names))
report={'release':meta['release'],'paper':meta['paper'],'status':'PASS' if all(x['status']=='PASS' for x in results) else 'FAIL','predictor_refit_tasks':[name for name in names if tasks[name].get('refits_predictors',False)],'precision':'compatible cache restoration + unchanged scoring tolerances' if a.compatible else 'exact cache restoration + task-specific numeric tolerances','tasks':results}
dump(out/'RUN_REPORT.json',report)
raise SystemExit(0 if report['status']=='PASS' else 1)
