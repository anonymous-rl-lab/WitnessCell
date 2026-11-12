from pathlib import Path
import concurrent.futures,hashlib,json,time,urllib.request,os
ROOT=Path(__file__).resolve().parent
URL="https://zenodo.org/records/15587986/files/cache.tar.gz?download=1"
SIZE=3392114905
CHUNK=32*1024*1024
PARTS=ROOT/"cache_parts";PARTS.mkdir(exist_ok=True)
def get(i):
 start=i*CHUNK;end=min(SIZE,start+CHUNK)-1;expected=end-start+1;p=PARTS/f"part_{i:04d}"
 if p.exists() and p.stat().st_size==expected:return i,expected,0.
 for attempt in range(4):
  t=time.monotonic();tmp=p.with_suffix(".tmp")
  try:
   req=urllib.request.Request(URL,headers={"Range":f"bytes={start}-{end}"})
   with urllib.request.urlopen(req,timeout=60) as r,tmp.open("wb") as f:
    if r.status!=206 or r.headers.get("Content-Range")!=f"bytes {start}-{end}/{SIZE}":raise ValueError("range mismatch")
    while True:
     block=r.read(1024*1024)
     if not block:break
     f.write(block)
   if tmp.stat().st_size!=expected:raise ValueError("size mismatch")
   tmp.replace(p);return i,expected,time.monotonic()-t
  except Exception as e:
   print(json.dumps({"part":i,"attempt":attempt,"error":str(e)}),flush=True)
   if attempt==3:raise
   time.sleep(min(2**attempt,8))
N=(SIZE+CHUNK-1)//CHUNK;t=time.monotonic();completed=0;total=0;records=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
 fs=[ex.submit(get,i) for i in range(N)]
 for f in concurrent.futures.as_completed(fs):
  i,n,dt=f.result();completed+=1;total+=n;records.append({"part":i,"bytes":n,"seconds":dt})
  info={"status":"downloading","completed_parts":completed,"total_parts":N,"bytes":total,"total_bytes":SIZE,"seconds":time.monotonic()-t}
  (ROOT/"DOWNLOAD_STATUS.json").write_text(json.dumps(info,indent=2)+"\n")
  print(json.dumps(info),flush=True)
out=ROOT/"cache.tar.gz";h=hashlib.md5()
with out.open("wb") as f:
 for i in range(N):
  with (PARTS/f"part_{i:04d}").open("rb") as r:
   while True:
    b=r.read(1024*1024)
    if not b:break
    f.write(b);h.update(b)
if out.stat().st_size!=SIZE or h.hexdigest()!="269273abc4f1d0da500898e35c514889":raise ValueError("archive checksum mismatch")
info={"status":"PASS","bytes":SIZE,"md5":h.hexdigest(),"seconds":time.monotonic()-t,"chunks":records}
(ROOT/"DOWNLOAD_STATUS.json").write_text(json.dumps(info,indent=2)+"\n")
print(json.dumps({k:v for k,v in info.items() if k!="chunks"}),flush=True)
