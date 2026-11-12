"""Single-GPU queue; launch only explicit jobs after the numerical probe passes."""
from pathlib import Path
import argparse,fcntl,hashlib,json,os,subprocess,sys,time,zipfile
R=Path(__file__).resolve().parent
def write(p,x):
    q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(x,indent=2));q.replace(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def archive(fold,seed):
    job=R/'results/detected49'/f'gears_seed{seed}'/f'fold_{fold}'
    subprocess.run([sys.executable,str(R/'archive_completed_jobs.py'),'--verify-one',str(job)],cwd=R,check=True)
    target=R/'gpu_completed';target.mkdir(exist_ok=True)
    dest=target/f'gears_seed{seed}_fold{fold}.zip';manifest={}
    files=[p for p in job.rglob('*') if p.is_file() and (p.suffix=='.json' or p.name=='predictions.npz' or p==job/'refit/last.pt')]
    files += [R/n for n in ['gpu_adapter_20261107.py','train_cell_models.py','train_cell_models_runtime.py','lazy_cell_data.py','fast_gears_ops.py','PROTOCOL.md','CUDA_PROBE.json']]
    with zipfile.ZipFile(dest.with_suffix('.tmp'),'w',zipfile.ZIP_DEFLATED,compresslevel=1) as z:
        for p in files:
            name=str(p.relative_to(R));z.write(p,name);manifest[name]={'sha256':sha(p),'bytes':p.stat().st_size}
        z.writestr('MANIFEST.json',json.dumps(manifest,indent=2))
    dest.with_suffix('.tmp').replace(dest)
    with zipfile.ZipFile(dest) as z:
        for n,m in manifest.items():assert hashlib.sha256(z.read(n)).hexdigest()==m['sha256']
    return dict(fold=fold,seed=seed,archive=str(dest),sha256=sha(dest),bytes=dest.stat().st_size,unix=time.time())
def main():
    p=argparse.ArgumentParser();p.add_argument('--jobs',required=True,help='fold:seed,fold:seed');a=p.parse_args()
    lock=(R/'GPU_QUEUE.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert json.loads((R/'CUDA_PROBE.json').read_text())['status']=='PASS'
    jobs=[tuple(map(int,x.split(':'))) for x in a.jobs.split(',')]
    done=[];logdir=R/'gpu_logs';logdir.mkdir(exist_ok=True)
    for i,(fold,seed) in enumerate(jobs):
        out=R/'results/detected49'/f'gears_seed{seed}'/f'fold_{fold}'
        write(R/'GPU_STATUS.json',dict(status='RUNNING',active=dict(fold=fold,seed=seed),complete=done,remaining=len(jobs)-i,unix=time.time()))
        if not (out/'FREEZE.json').exists():
            env=os.environ.copy();env.update(OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',MALLOC_ARENA_MAX='2')
            with (logdir/f'gears_seed{seed}_fold{fold}.log').open('a') as log:
                proc=subprocess.Popen([sys.executable,'-u',str(R/'gpu_adapter_20261107.py'),'--kind','gears','--fold',str(fold),'--seed',str(seed),'--threads','2'],cwd=R,env=env,stdout=log,stderr=subprocess.STDOUT)
                while proc.poll() is None:
                    write(R/'GPU_STATUS.json',dict(status='RUNNING',active=dict(fold=fold,seed=seed,pid=proc.pid),complete=done,remaining=len(jobs)-i,unix=time.time()))
                    time.sleep(10)
                if proc.returncode:
                    write(R/'GPU_STATUS.json',dict(status='BLOCKED',returncode=proc.returncode,active=dict(fold=fold,seed=seed),complete=done,unix=time.time()));return proc.returncode
        done.append(archive(fold,seed));write(R/'GPU_COMPLETED.json',done)
    write(R/'GPU_STATUS.json',dict(status='COMPLETE',complete=done,unix=time.time()));return 0
if __name__=='__main__':raise SystemExit(main())
