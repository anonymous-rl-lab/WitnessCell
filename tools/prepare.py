"""Restore exact duplicate paths and deterministic Han caches. No fitting or labels."""
from pathlib import Path
import json,hashlib,shutil,subprocess,sys,argparse,gzip
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def restore_packed():
 path=R/'provenance/PACKED_ASSETS.json'
 rows=json.loads(path.read_text()) if path.exists() else []
 for row in rows:
  src,dst=R/row['packed'],R/row['path']
  for p in (src,dst):
   if not p.resolve().is_relative_to(R) or p.is_symlink():raise ValueError('Unsafe packed path')
  if sha(src)!=row['packed_sha256']:raise ValueError('Packed checksum mismatch')
  if dst.exists():
   if sha(dst)!=row['sha256']:raise ValueError('Existing unpacked payload differs')
  else:
   dst.parent.mkdir(parents=True,exist_ok=True)
   tmp=dst.with_name(dst.name+'.restoring')
   with gzip.open(src,'rb') as inp,tmp.open('wb') as output:shutil.copyfileobj(inp,output)
   if tmp.stat().st_size!=row['bytes'] or sha(tmp)!=row['sha256']:
    tmp.unlink();raise ValueError('Restored payload differs')
   tmp.replace(dst)
 return len(rows)
def restore_aliases():
 rows=json.loads((R/'provenance/ALIASES.json').read_text()) if (R/'provenance/ALIASES.json').exists() else []
 for row in rows:
  src,dst=R/row['canonical'],R/row['path']
  for p in (src,dst):
   if not p.resolve().is_relative_to(R) or p.is_symlink():raise ValueError('Unsafe alias path')
  if sha(src)!=row['sha256']:raise ValueError('Canonical checksum mismatch: '+row['canonical'])
  if dst.exists():
   if sha(dst)!=row['sha256']:raise ValueError('Existing alias differs: '+row['path'])
  else:dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
 return len(rows)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--compatible',action='store_true');ap.add_argument('--aliases-only',action='store_true');a=ap.parse_args();out=a.out.resolve()
 if out.is_relative_to(R):ap.error('Reports belong outside the release')
 out.mkdir(parents=True,exist_ok=True);np=restore_packed();n=restore_aliases()
 if not a.aliases_only:
  cmd=[sys.executable,'-B',str(R/'tools/restore_han_cache.py'),'--root',str(R/'experiments/han_native'),'--manifest',str(R/'provenance/HAN_DERIVED_CACHES.json'),'--out',str(out/'HAN_CACHE_RESTORE.json')]
  if a.compatible:cmd+=['--compatible']
  subprocess.run(cmd,check=True)
 (out/'PREPARE.json').write_text(json.dumps({'status':'PASS','aliases_restored_or_verified':n,'packed_assets_restored_or_verified':np,'han_cache_rebuilt':not a.aliases_only},indent=2)+'\n')
 print('Prepared',n,'duplicate logical paths')
if __name__=='__main__':main()
