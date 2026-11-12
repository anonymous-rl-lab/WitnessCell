#!/usr/bin/env python3
"""Offline full follow-up rerun using lossless condition-mean inputs only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse,csv,hashlib,json,subprocess,sys,time,datetime
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')
def compare_table(left,right):
    with left.open() as f:a=list(csv.DictReader(f))
    with right.open() as f:b=list(csv.DictReader(f))
    assert len(a)==len(b),(left,'row count')
    largest=0.
    for x,y in zip(a,b):
        assert x.keys()==y.keys(),(left,'columns')
        for key,value in x.items():
            try:v,w=float(value),float(y[key])
            except (ValueError,TypeError):assert value==y[key],(left,key);continue
            diff=abs(v-w);largest=max(largest,diff)
            assert diff<=2e-11*max(1.,abs(v)),(left,key,value,y[key])
    return {'rows':len(a),'max_abs_difference':largest}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    out=a.output.resolve()
    if out.exists() and any(out.iterdir()):raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True,exist_ok=True);logs=out/'logs';logs.mkdir()
    for item in json.loads((ROOT/'SOURCE_MANIFEST.json').read_text()):
        assert sha(ROOT/item['path'])==item['sha256'],('input/source changed',item['path'])
    source=ROOT/'source_helpers';common=['--compact-root',str(ROOT/'compact_assets'),'--bridge-code',str(source/'run_bridge.py'),'--geometry-code',str(source/'run_geometry_probe.py')]
    modelsource=source/'experiment_assets/WitnessCell_Partner_Component_Run_20261018/model/global_pair_model.py'
    steps=[
        ('postfit_reference','model/run_reference.py',['--asset-root',str(source),*common,'--out',str(out/'model')]),
        ('reference_kernel','kernel_probe/run_reference_kernel.py',['--model-code',str(modelsource),*common,'--out',str(out/'kernel')]),
        ('algebra','theory_audit/check_reference_algebra.py',['--out',str(out/'algebra.json')]),
        ('independent_real_audit','theory_audit/audit_real_reference.py',['--model-dir',str(out/'model'),'--kernel-dir',str(out/'kernel'),'--out',str(out/'audit')]),
        ('export_single_query','model/export_deploy.py',['--source-run',str(out/'model'),'--out',str(out/'deploy')])]
    bootstrap='''import os,sys,runpy\nfrom pathlib import Path\nallowed=[Path(sys.argv.pop(1)).resolve(),Path(sys.argv.pop(1)).resolve()]\ndef guard(event,args):\n if event=="open" and isinstance(args[0],(str,bytes,os.PathLike)):\n  p=Path(os.fsdecode(args[0])).resolve()\n  if str(p).startswith("/workspace/scratch/") and not any(p==r or r in p.parents for r in allowed):raise RuntimeError("Unexpected historical workspace read: "+str(p))\nsys.addaudithook(guard)\nscript=sys.argv.pop(1);sys.argv[0]=script;sys.path.insert(0,str(Path(script).parent));runpy.run_path(script,run_name="__main__")'''
    results=[];start=time.monotonic()
    for name,script,args in steps:
        tick=time.monotonic();print(json.dumps({'start':name}),flush=True)
        with (logs/(name+'.log')).open('w') as handle:
            p=subprocess.run([sys.executable,'-c',bootstrap,str(ROOT),str(out),str(ROOT/script),*args],cwd=ROOT,stdout=handle,stderr=subprocess.STDOUT,timeout=180)
        if p.returncode:raise RuntimeError(f'{name} failed; inspect {logs/(name+".log")}')
        results.append({'stage':name,'seconds':time.monotonic()-tick});print(json.dumps(results[-1]),flush=True)
    expected=json.loads((ROOT/'EXPECTED_REPRODUCTION.json').read_text());verified=[]
    for record in expected['arrays']:
        with np.load(out/record['output_file'],allow_pickle=False) as z:
            value=z[record['key']]
            assert list(value.shape)==record['shape'],record
            assert str(value.dtype)==record['dtype'],record
            actual=hashlib.sha256(value.tobytes()).hexdigest()
            assert actual==record['sha256'],('prediction/state array differs',record['output_file'],record['key'])
        verified.append((record['output_file'],record['key']))
    tables={}
    for entry in expected['tables']:
        tables[entry['output_file']]=compare_table(ROOT/entry['reference_file'],out/entry['output_file'])
    assert json.loads((out/'model/COMPLETE.json').read_text())['status']=='PASS'
    assert json.loads((out/'audit/REAL_AUDIT_RECEIPT.json').read_text())['status']=='PASS'
    receipt={'status':'PASS','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'seconds':time.monotonic()-start,
             'stages':results,'exact_named_arrays':len(verified),'tables':tables,
             'network_required':False,'historical_workspace_reads_blocked':True,
             'large_fit_states_regenerated_from_exact_condition_means':True}
    dump(out/'REPRODUCTION_COMPLETE.json',receipt);print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()
