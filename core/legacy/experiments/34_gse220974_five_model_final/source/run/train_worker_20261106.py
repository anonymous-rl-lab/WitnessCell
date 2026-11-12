"""Operational wrapper: thread policy, job lock, progress and safe checkpoint exit.

Imports unchanged training code. Model, optimizer, batches, RNG, budgets and
selection are inherited. A requested pause exits only after an atomic partial
checkpoint; it does not interrupt a completed epoch before its early-stop check.
"""
from pathlib import Path
import argparse,fcntl,hashlib,json,os,signal,sys,time
R=Path(__file__).resolve().parent

def atomic(path,value):
 tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(value,indent=2));tmp.replace(path)

def main():
 p=argparse.ArgumentParser();p.add_argument('--kind',choices=['gears','scouter'],required=True)
 p.add_argument('--fold',type=int,required=True);p.add_argument('--seed',type=int,required=True)
 p.add_argument('--threads',type=int,required=True);a=p.parse_args()
 out=R/'results/detected49'/f'{a.kind}_seed{a.seed}'/f'fold_{a.fold}';out.mkdir(parents=True,exist_ok=True)
 lock=(out/'WORKER.lock').open('a')
 try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 except BlockingIOError:print('Job already locked; no duplicate worker.',flush=True);return 73
 stop=out/'STOP_WORKER';pause=[False];steps=[0];pending={};started=time.time()
 def request_pause(*_):pause[0]=True
 signal.signal(signal.SIGTERM,request_pause);signal.signal(signal.SIGINT,request_pause)
 import train_cell_models_runtime as t
 original_threads=t.torch.set_num_threads
 t.torch.set_num_threads=lambda ignored:original_threads(a.threads)
 original_threads(a.threads)
 original_step=t.torch.optim.Adam.step
 original_save=t.torch.save
 original_replace=Path.replace
 metadata=dict(kind=a.kind,fold=a.fold,seed=a.seed,threads=a.threads,pid=os.getpid(),started_unix=started,
  wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),runtime_fingerprint=t.runtime_fingerprint(),training_rules_changed=False)
 atomic(out/'WORKER_RUNTIME.json',metadata)
 def heartbeat(event,**kwargs):
  atomic(out/'WORKER_HEARTBEAT.json',dict(event=event,steps_in_process=steps[0],unix=time.time(),**kwargs))
 def step(self,*args,**kwargs):
  result=original_step(self,*args,**kwargs);steps[0]+=1
  if steps[0]%10==0:heartbeat('TRAINING')
  return result
 def save(obj,path,*args,**kwargs):
  if isinstance(obj,dict) and 'partial' in obj and str(path).endswith('resume.pt.tmp'):
   partial=obj.get('partial')
   pending[str(Path(path))]=dict(partial=partial is not None,epoch_completed=obj['epoch'],
    next_start=None if partial is None else partial['next_start'],batches_done=None if partial is None else len(partial['losses']))
  return original_save(obj,path,*args,**kwargs)
 def replace(self,target):
  result=original_replace(self,target)
  info=pending.pop(str(self),None)
  if info is not None:
   heartbeat('CHECKPOINT',checkpoint=str(target),**info)
   if info['partial'] and (pause[0] or stop.exists()):
    heartbeat('PAUSED_AT_CHECKPOINT',checkpoint=str(target),**info)
    print(json.dumps(dict(event='PAUSED_AT_CHECKPOINT',checkpoint=str(target),**info)),flush=True)
    raise SystemExit(75)
  return result
 t.torch.optim.Adam.step=step;t.torch.save=save;Path.replace=replace
 heartbeat('STARTING')
 try:t.run(a.fold,a.kind,a.seed)
 finally:
  # Do not let final metadata overwrite PAUSED_AT_CHECKPOINT.
  t.torch.optim.Adam.step=original_step;t.torch.save=original_save;Path.replace=original_replace
 if (out/'FREEZE.json').exists():heartbeat('FROZEN')
 return 0
if __name__=='__main__':sys.exit(main())
