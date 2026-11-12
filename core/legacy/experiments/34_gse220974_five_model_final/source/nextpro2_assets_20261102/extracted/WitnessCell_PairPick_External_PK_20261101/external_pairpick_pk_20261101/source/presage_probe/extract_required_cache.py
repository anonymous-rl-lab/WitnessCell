from pathlib import Path
import tarfile,json,hashlib,shutil,time
ROOT=Path(__file__).resolve().parent
status=json.loads((ROOT/'DOWNLOAD_STATUS.json').read_text())
assert status['status']=='PASS'
lines=(ROOT/'upstream/sample_files/prior_files/sample.knowledge_experimental.txt').read_text().splitlines()
required={x.strip().replace('../','',1) for x in lines if x.strip()}
assert len(required)==40
found=[];t=time.monotonic()
with tarfile.open(ROOT/'cache.tar.gz','r|gz') as z:
 for item in z:
  if item.name not in required:continue
  assert item.isfile() and not item.name.startswith('/') and '..' not in Path(item.name).parts
  dest=ROOT/item.name;dest.parent.mkdir(parents=True,exist_ok=True)
  h=hashlib.sha256();n=0
  with z.extractfile(item) as src,dest.open('wb') as out:
   while True:
    block=src.read(1024*1024)
    if not block:break
    out.write(block);h.update(block);n+=len(block)
  assert n==item.size
  found.append(dict(path=item.name,bytes=n,sha256=h.hexdigest()))
  print(json.dumps(dict(extracted=len(found),total=len(required),path=item.name,seconds=time.monotonic()-t)),flush=True)
assert {x['path'] for x in found}==required
(ROOT/'REQUIRED_CACHE_MANIFEST.json').write_text(json.dumps(dict(status='PASS',source_md5=status['md5'],files=found,seconds=time.monotonic()-t),indent=2)+'\n')
print('PASS extracted required 40 source files only',flush=True)
