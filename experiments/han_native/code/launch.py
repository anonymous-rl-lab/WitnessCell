from common import *
import argparse,subprocess,sys,time,concurrent.futures

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=6);ap.add_argument('--tasks',default='all');a=ap.parse_args();tasks=json.loads((ROOT/'TASKS.json').read_text());ids=list(range(len(tasks))) if a.tasks=='all' else [int(x) for x in a.tasks.split(',')]
    (ROOT/'logs').mkdir(exist_ok=True)
    def run(i):
        with (ROOT/'logs'/f'task{i:02d}.log').open('a') as f:
            f.write('\n'+now()+' launch\n');f.flush();p=subprocess.run([sys.executable,str(ROOT/'code/fit_task.py'),'--task',str(i)],stdout=f,stderr=subprocess.STDOUT)
        return dict(task=i,exit_code=p.returncode,utc=now())
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as ex:
        futures={ex.submit(run,i):i for i in ids}
        for fut in concurrent.futures.as_completed(futures):
            r=fut.result();results.append(r);dump(ROOT/'BATCH_PROGRESS.json',dict(completed=len(results),requested=len(ids),workers=a.workers,results=results));print(json.dumps(r),flush=True)
    dump(ROOT/'BATCH_FINISHED.json',dict(utc=now(),results=results))
if __name__=='__main__':main()
