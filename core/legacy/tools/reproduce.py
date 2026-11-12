#!/usr/bin/env python3
"""Bounded offline replay of current-paper evidence; no downloads or tuning."""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--task',default='all',help='all or an identifier in REPRODUCTION_INDEX.json')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--seconds-per-task',type=float,default=None,
                   help='Override task timeout; default 300 seconds, or the task-specific limit')
    a=p.parse_args();out=a.out.resolve()
    if out.exists() or out.is_relative_to(ROOT):p.error('Output must be new and outside release')
    if a.seconds_per_task is not None and (not math.isfinite(a.seconds_per_task) or a.seconds_per_task<=0):p.error('Time limit must be finite and positive')
    tasks=json.loads((ROOT/'REPRODUCTION_INDEX.json').read_text())['tasks']
    if a.task!='all' and a.task not in tasks:p.error('Unknown task; choose '+', '.join(tasks))
    selected=list(tasks) if a.task=='all' else [a.task];out.mkdir(parents=True)
    env=os.environ.copy();env.pop('PYTHONPATH',None)
    env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONOPTIMIZE='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',MPLBACKEND='Agg')
    report={'status':'RUNNING','release':(ROOT/'VERSION').read_text().strip(),'paper':'22','tasks':[]}
    for key in selected:
        t=tasks[key];args=[s.replace('{root}',str(ROOT)).replace('{out}',str(out/key)) for s in t['args']]
        cmd=[sys.executable,'-B',str(ROOT/t['script']),*args];row={'task':key,'scope':t['scope'],'command':cmd,'log':key+'.log'}
        limit=float(a.seconds_per_task if a.seconds_per_task is not None else t.get('timeout_seconds',300))
        if not math.isfinite(limit) or limit<=0:raise ValueError('Invalid per-task timeout '+key)
        row['timeout_seconds']=limit
        print('Running '+key,flush=True);start=time.monotonic()
        with (out/row['log']).open('w') as log:
            proc=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=(os.name=='posix'))
            try:row['returncode']=proc.wait(timeout=limit);row['status']='PASS' if row['returncode']==0 else 'FAIL'
            except subprocess.TimeoutExpired:
                if os.name=='posix':os.killpg(proc.pid,signal.SIGTERM)
                else:proc.terminate()
                try:proc.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    if os.name=='posix':os.killpg(proc.pid,signal.SIGKILL)
                    else:proc.kill()
                    proc.wait()
                row.update(status='TIMEOUT',returncode=proc.returncode)
        row['seconds']=time.monotonic()-start;report['tasks'].append(row)
        (out/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n');print(key+': '+row['status'],flush=True)
    report['passed']=sum(x['status']=='PASS' for x in report['tasks']);report['status']='PASS' if report['passed']==len(selected) else 'FAIL'
    (out/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':report['status'],'passed':report['passed'],'total':len(selected),'out':str(out)},indent=2))
    if report['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
