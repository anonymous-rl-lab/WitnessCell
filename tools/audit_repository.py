"""Read-only repository size, duplicate, hygiene and entry-point audit."""
from pathlib import Path
from collections import defaultdict
import argparse,hashlib,json,re
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
if a.out.resolve().is_relative_to(R):ap.error('Reports belong outside the repository')
files=[p for p in R.rglob('*') if p.is_file() and '.git' not in p.relative_to(R).parts]
groups=defaultdict(list);secrets=[];cache=[];active_paths=[]
pattern=re.compile(r'AKIA[0-9A-Z]{16}|gh[pousr]_[A-Za-z0-9]{36}|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----')
for p in files:
 rel=str(p.relative_to(R));b=p.read_bytes();groups[(len(b),hashlib.sha256(b).hexdigest())].append(rel)
 if p.suffix=='.pyc' or '__pycache__' in p.parts:cache.append(rel)
 if p.suffix in {'.py','.md','.txt','.json','.yml','.yaml','.sh','.toml','.cfg'} and len(b)<10_000_000:
  text=b.decode('utf-8',errors='replace')
  for n,line in enumerate(text.splitlines(),1):
   if pattern.search(line):secrets.append(dict(path=rel,line=n))
  # This scanner contains the literal detection signatures, not runtime locations.
  if rel!='tools/audit_repository.py' and (rel.startswith('tools/') or rel=='reproduce.py') and re.search(r'/workspace/|/mnt/data/|/home/agent/',text):active_paths.append(rel)
dups=[dict(bytes=size,sha256=h,paths=names) for (size,h),names in groups.items() if len(names)>1 and size>0]
large=[dict(path=str(p.relative_to(R)),bytes=p.stat().st_size) for p in files if p.stat().st_size>50*1024**2]
meta=json.loads((R/'REPRODUCTION_INDEX.json').read_text());missing=[k for k,v in meta['tasks'].items() if not (R/v['script']).is_file()]
report=dict(status='PASS' if not secrets and not cache and not active_paths and not large and not missing else 'FAIL',files=len(files),bytes=sum(p.stat().st_size for p in files),largest_file=max((dict(path=str(p.relative_to(R)),bytes=p.stat().st_size) for p in files),key=lambda x:x['bytes']),over_50MiB=large,duplicate_groups=dups,duplicate_bytes=sum(x['bytes']*(len(x['paths'])-1) for x in dups),interpreter_caches=cache,credential_pattern_findings=secrets,absolute_workspace_paths_in_active_wrappers=active_paths,missing_task_entries=missing,scope='Static hygiene and physical duplication scan; no claim of exhaustive secret detection or proof that every historical source module is needed.')
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='duplicate_groups'},indent=2));raise SystemExit(0 if report['status']=='PASS' else 1)
