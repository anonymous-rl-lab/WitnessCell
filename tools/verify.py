"""Verify physical assets, exact restored aliases and current companion provenance."""
from pathlib import Path
import argparse,csv,gzip,hashlib,json
R=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def safe(path):
 p=R/path
 if not p.resolve().is_relative_to(R) or p.is_symlink():raise ValueError('Unsafe release path: '+path)
 return p
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);ap.add_argument('--restored',action='store_true');ap.add_argument('--compatible',action='store_true');a=ap.parse_args()
 if a.out and a.out.resolve().is_relative_to(R):ap.error('Write reports outside the repository')
 rows=json.loads((R/'MANIFEST.json').read_text())['files']
 for row in rows:
  p=safe(row['path']);assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],row['path']
 packed=json.loads((R/'provenance/PACKED_ASSETS.json').read_text()) if (R/'provenance/PACKED_ASSETS.json').exists() else []
 packmap={x['path']:x for x in packed}
 for row in packed:
  p=safe(row['packed']);assert sha(p)==row['packed_sha256']
  if a.restored:assert safe(row['path']).exists() and sha(safe(row['path']))==row['sha256']
 aliases=json.loads((R/'provenance/ALIASES.json').read_text());derived=json.loads((R/'provenance/HAN_DERIVED_CACHES.json').read_text())
 for row in aliases:
  p=safe(row['path']);c=safe(row['canonical'])
  if c.exists():assert sha(c)==row['sha256']
  else:assert row['canonical'] in packmap and packmap[row['canonical']]['sha256']==row['sha256']
  if a.restored:assert p.exists() and sha(p)==row['sha256'],row['path']
 if a.restored:
  for row in derived:
   p=safe('experiments/han_native/'+row['path']);assert p.exists()
   if not a.compatible:assert sha(p)==row['sha256'],row['path']
 count=0
 if a.restored:
  with (R/'provenance/V41_INTEGRATION_MAP.csv').open() as f:
   for row in csv.DictReader(f):
    p=safe(row['integrated_path']);assert p.stat().st_size==int(row['bytes']) and sha(p)==row['sha256'],row['source_path'];count+=1
 report=dict(status='PASS',release='31.0.0',paper='v41',language_revision='R1',physical_files_verified=len(rows),exact_aliases=len(aliases),restored_aliases_verified=len(aliases) if a.restored else 0,v41_source_files_verified=count,packed_assets=len(packed),derived_caches=len(derived),cache_byte_identity_checked=a.restored and not a.compatible,scope='File integrity; numerical acceptance is reported by reproduce.py tasks.')
 if a.out:a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
