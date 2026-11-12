"""Bounded two-worker supervisor; no test scoring until all risk outputs freeze."""
import os,sys,time,json,hashlib,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
MAX_WALL=1200
WORKERS=2
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n')
def emit(**x):print(json.dumps(x),flush=True)
assert json.loads((HERE/'projection_subsample_equivalence.json').read_text())['non_target_float32_exact']
assert json.loads((HERE/'results/Norman_seed1/projection_equivalence.json').read_text())['non_target_float32_exact']
assert json.loads((HERE/'results/Norman_seed1/replay.json').read_text())['float32_exact']
for ds in ['Norman','Wessels']:assert (HERE/f'results/{ds}_seed1/pilot.json').exists()
sources={p.name:sha(p) for p in sorted(HERE.glob('*.py'))}
write('execution_lock.json',{'protocol':'WitnessCell-step1-v1-20261017','frozen_at_unix':time.time(),'sources':sources,'workers':WORKERS,'cpu_threads_per_worker':1,'wall_limit_seconds':MAX_WALL,'formal_seeds':[1,2,3],'datasets':['Norman','Wessels'],'epochs':'not applicable to deterministic closed-form model','scouter':'pending; no available vector/training input','estimated_formal_minutes':[5,10],'labels_opened_before_freeze':False})
pending=[(d,s) for d in ['Norman','Wessels'] for s in [1,2,3]];active={};done=[];started=time.monotonic();last=0
try:
 while pending or active:
  if time.monotonic()-started>MAX_WALL:raise TimeoutError('Formal wall budget exceeded')
  while pending and len(active)<WORKERS:
   ds,seed=pending.pop(0);key=f'{ds}_seed{seed}';log=(HERE/f'formal_{key}.log').open('w')
   proc=subprocess.Popen([sys.executable,'-u',str(HERE/'run_step1.py'),'formal','--dataset',ds,'--seed',str(seed)],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
   active[key]=(proc,log,time.monotonic());emit(event='started',job=key,pid=proc.pid)
  for key,(proc,log,t) in list(active.items()):
   rc=proc.poll()
   if rc is not None:
    log.close();del active[key]
    if rc:raise RuntimeError(f'{key} failed with exit {rc}; see formal_{key}.log')
    done.append(key);emit(event='finished',job=key,seconds=round(time.monotonic()-t,2))
  if time.monotonic()-last>=10:
   state={'status':'running','elapsed_seconds':round(time.monotonic()-started,1),'completed':done,'active':list(active),'pending':pending};write('formal_status.json',state);emit(**state);last=time.monotonic()
  time.sleep(1)
 for name,digest in sources.items():assert sha(HERE/name)==digest,('source_changed',name)
 outputs={d.name:sha(d/'frozen_query_outputs.csv') for d in sorted((HERE/'results').iterdir()) if d.is_dir()}
 write('all_risks_frozen.json',{'frozen_at_unix':time.time(),'outputs':outputs,'elapsed_seconds':time.monotonic()-started,'test_truth_scored':False})
 write('formal_status.json',{'status':'all_six_risk_outputs_frozen','elapsed_seconds':time.monotonic()-started,'completed':done,'next':'separate label evaluation and independent aggregation'})
 emit(event='all_six_frozen',seconds=round(time.monotonic()-started,2))
except BaseException as exc:
 for proc,log,t in active.values():
  if proc.poll() is None:proc.terminate()
  log.close()
 write('formal_status.json',{'status':'failed_or_bounded_stop','error':str(exc),'completed':done,'elapsed_seconds':time.monotonic()-started})
 raise
