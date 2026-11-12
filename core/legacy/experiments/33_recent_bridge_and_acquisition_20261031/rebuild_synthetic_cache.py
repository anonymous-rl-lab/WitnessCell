"""Regenerate deterministic synthetic NPZ caches with ORIGINAL scientific code.

Only NPZ cache files are restored. Original JSON, timing, protocol and hash
manifests remain unchanged. Every restored file must match the original SHA256.
No new tuning, changed seeds or evaluation-label fitting is performed.

python rebuild_synthetic_cache.py --restored /path/to/restored --restore
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
import argparse,hashlib,json,shutil,subprocess,sys,tempfile,time
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def locate(root,name,alias):
    for n in (name,alias):
        p=root/n
        if p.is_dir():return p
    raise FileNotFoundError(name)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--restored',type=Path,required=True)
    ap.add_argument('--restore',action='store_true',required=True)
    ap.add_argument('--report',type=Path)
    ap.add_argument('--only',choices=('all','formal','causal'),default='all')
    args=ap.parse_args();root=args.restored.resolve();started=time.monotonic();results=[]
    outreport=args.report or root/'SYNTHETIC_CACHE_REBUILD_RESULT.json'
    with tempfile.TemporaryDirectory(prefix='witnesscell_synthetic_rebuild_') as td:
        temp=Path(td)
        def run(path,argv):
            log=temp/(path.parent.name+'_'+path.stem+'.log')
            with log.open('w') as f:
                r=subprocess.run([sys.executable,'-B',str(path),*map(str,argv)],cwd=path.parent,
                                 stdout=f,stderr=subprocess.STDOUT)
            assert r.returncode==0,log.read_text()[-12000:]
        def restore_files(src,generated,expected,label):
            count=0;already=0
            for name,h in expected.items():
                if not name.endswith('.npz'):continue
                p=generated/name;assert p.is_file(),str(p)
                assert sha(p)==h,(label,name,'regeneration SHA mismatch',sha(p),h)
                target=src/name
                if target.exists():assert sha(target)==h,(label,name,'existing SHA mismatch');already+=1
                else:target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
                count+=1
            return {'stage':label,'npz_sha256_exact':count,'already_present_exact':already,
                    'original_json_and_manifest_unchanged':True}
        if args.only in ('all','formal'):
            original=locate(root,'WitnessCell_Formal_Hidden_Mechanism_Bridge_20261031','formal_bridge_20261031')
            manifest=read(original/'frozen/FREEZE_MANIFEST.json')
            for name,h in manifest['dependencies'].items():assert sha(original/name)==h,name
            work=temp/'formal';work.mkdir()
            for name in set(manifest['dependencies'])|{'runner.py','model.py','CONFIG.json'}:
                p=work/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(original/name,p)
            run(work/'runner.py',[])
            results.append(restore_files(original/'frozen',work/'frozen',manifest['files'],'formal'))
            assert results[-1]['npz_sha256_exact']==1
            print(json.dumps(results[-1]),flush=True)
        if args.only in ('all','causal'):
            bridge=locate(root,'WitnessCell_NM_NCS_Bridge_20261031','bridge_experiment_20261031')
            original=bridge/'causal_run';manifest=read(original/'formal/FREEZE_MANIFEST.json')
            for name,h in manifest['dependencies'].items():assert sha(original/name)==h,name
            work=temp/'causal_generated'
            run(original/'formal_run.py',['--output',work])
            results.append(restore_files(original/'formal',work,manifest['files'],'causal'))
            assert results[-1]['npz_sha256_exact']==288
            print(json.dumps(results[-1]),flush=True)
    result={'status':'PASS','original_npz_sha256_exact':sum(r['npz_sha256_exact'] for r in results),
            'results':results,'seconds':time.monotonic()-started,'python':sys.version,
            'scientific_sources_modified':False,'evaluation_labels_used_in_fitting':False}
    outreport.parent.mkdir(parents=True,exist_ok=True);outreport.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
