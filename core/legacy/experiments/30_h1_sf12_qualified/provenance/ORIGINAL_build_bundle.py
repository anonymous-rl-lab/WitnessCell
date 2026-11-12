from pathlib import Path
import shutil,json,hashlib,zipfile,subprocess,sys
R=Path(__file__).resolve().parent
name='WitnessCell_Qualified_Rectangles_Fixed_Test_20261026'
dest=R/'deliverables'/name
dest.mkdir(parents=True,exist_ok=True)
files=[R/n for n in ['REPORT_CN.txt','README.txt','run_five.py','score_five.py','replay_compact.py','build_bundle.py']]
files+=list((R/'source').glob('*.py'))+list((R/'source').glob('*.json'))
files+=list((R/'audit').glob('*.py'))+list((R/'audit').glob('*.json'))
files+=[R/'run'/n for n in ['PROTOCOL.json','SMOKE.json','ALL_PREDICTIONS_FROZEN.json','SCORES.json','EVALUATION_STATS.npz','COMPACT_VECTORS.npz']]
for i in range(1,6):files +=[R/'run'/f'Q{i}'/n for n in ['TRACE.json','TRAIN_STATS.npz']]
for p in files:
    q=dest/p.relative_to(R);q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
proc=subprocess.run([sys.executable,'replay_compact.py'],cwd=dest,text=True,capture_output=True)
assert proc.returncode==0,proc.stderr
(dest/'COMPACT_REPLAY.json').write_text(proc.stdout)
manifest={str(p.relative_to(dest)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(dest.rglob('*')) if p.is_file() and p.name!='PACKAGE_MANIFEST.json'}
(dest/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
zp=dest.with_suffix('.zip')
with zipfile.ZipFile(zp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for p in sorted(dest.rglob('*')):
        if p.is_file():z.write(p,str(p.relative_to(dest.parent)))
with zipfile.ZipFile(zp) as z:
    assert z.testzip() is None
    for n,m in manifest.items():assert hashlib.sha256(z.read(name+'/'+n)).hexdigest()==m['sha256']
    extracted=R/'deliverables/clean_verify';z.extractall(extracted)
proc=subprocess.run([sys.executable,'replay_compact.py'],cwd=extracted/name,text=True,capture_output=True)
assert proc.returncode==0,proc.stderr
print(json.dumps(dict(path=str(zp),bytes=zp.stat().st_size,sha256=hashlib.sha256(zp.read_bytes()).hexdigest(),files=len(manifest)+1,clean_replay=json.loads(proc.stdout)),indent=2))
