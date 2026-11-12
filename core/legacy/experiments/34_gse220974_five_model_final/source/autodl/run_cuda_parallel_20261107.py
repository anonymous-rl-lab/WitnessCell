"""Two GPU workers, adopting existing workers and awaiting CPU migration safely."""
from pathlib import Path
import argparse,fcntl,json,os,subprocess,sys,time
R=Path(__file__).resolve().parent
import run_cuda_queue_20261107 as q
def read(p,default=None):return json.loads(p.read_text()) if p.exists() else default
def folder(job):return R/'results/detected49'/f'gears_seed{job[1]}'/f'fold_{job[0]}'
def locked(job):
    p=folder(job)/'WORKER.lock';p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('a') as f:
        try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);return False
        except BlockingIOError:return True
def ready(job):return job not in [(3,0),(3,1)] or (R/'CPU_MIGRATION_READY.json').exists()
def main():
    lock=(R/'GPU_QUEUE.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert read(R/'CUDA_PROBE.json')['status']=='PASS'
    assert read(R/'CUDA_REFERENCE_REPLAY.json')['status']=='PASS'
    jobs=[(3,2),(3,0),(3,1),(4,0),(4,1),(4,2)]
    active={};done=read(R/'GPU_COMPLETED.json',[]);failures=[]
    complete={(x['fold'],x['seed']) for x in done}
    pending=[j for j in jobs if j not in complete]
    logs=R/'gpu_logs';logs.mkdir(exist_ok=True)
    for j in list(pending):
        if locked(j):active[j]={'process':None,'adopted':True};pending.remove(j)
    while pending or active:
        for j,a in list(active.items()):
            process=a['process'];alive=locked(j) if process is None else process.poll() is None
            if alive:continue
            del active[j]
            if (folder(j)/'FREEZE.json').exists():
                try:
                    row=q.archive(*j);done.append(row);complete.add(j);q.write(R/'GPU_COMPLETED.json',done)
                except Exception as e:failures.append(dict(fold=j[0],seed=j[1],error=repr(e)))
            else:failures.append(dict(fold=j[0],seed=j[1],returncode=None if process is None else process.returncode))
        candidates=[j for j in pending if ready(j)]
        while candidates and len(active)<2:
            j=candidates.pop(0);pending.remove(j)
            if (folder(j)/'FREEZE.json').exists():
                row=q.archive(*j);done.append(row);complete.add(j);q.write(R/'GPU_COMPLETED.json',done);continue
            if locked(j):active[j]={'process':None,'adopted':True};continue
            env=os.environ.copy();env.update(OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',MALLOC_ARENA_MAX='2')
            with (logs/f'gears_seed{j[1]}_fold{j[0]}.log').open('a') as log:
                proc=subprocess.Popen([sys.executable,'-u',str(R/'gpu_adapter_20261107.py'),'--kind','gears','--fold',str(j[0]),'--seed',str(j[1]),'--threads','2'],cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT)
            active[j]={'process':proc,'adopted':False}
        q.write(R/'GPU_STATUS.json',dict(status='RUNNING_WITH_BLOCKED_JOBS' if failures else 'RUNNING',worker_limit=2,
            active=[dict(fold=j[0],seed=j[1],adopted=a['adopted'],pid=None if a['process'] is None else a['process'].pid,heartbeat=read(folder(j)/'WORKER_HEARTBEAT.json',{})) for j,a in active.items()],
            queued=[dict(fold=j[0],seed=j[1],ready=ready(j)) for j in pending],complete=done,failures=failures,unix=time.time()))
        time.sleep(10)
    q.write(R/'GPU_STATUS.json',dict(status='COMPLETE' if not failures else 'BLOCKED',complete=done,failures=failures,unix=time.time()))
if __name__=='__main__':main()
