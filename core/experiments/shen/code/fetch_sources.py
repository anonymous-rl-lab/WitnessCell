#!/usr/bin/env python3
"""Fetch the original 2017 source workbook; no dependence on CellSimBench files."""
import argparse,hashlib,json,sys,time,zipfile
from pathlib import Path
from urllib.request import Request,urlopen
from datetime import datetime,timezone
FILES={
'Shen2017_TableS2.xlsx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.4225/MediaObjects/41592_2017_BFnmeth4225_MOESM42_ESM.xlsx',
'Shen2017_TableS1.xlsx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.4225/MediaObjects/41592_2017_BFnmeth4225_MOESM41_ESM.xlsx'}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--timeout',type=int,default=60);a=p.parse_args()
 a.out.mkdir(parents=True,exist_ok=True);receipt=[]
 for name,url in FILES.items():
  dst=a.out/name;item=dict(name=name,url=url,utc=datetime.now(timezone.utc).isoformat())
  try:
   if not dst.exists():
    tmp=dst.with_name(dst.name+'.part')
    with urlopen(Request(url,headers={'User-Agent':'WitnessCell-Shen-Replication/1.0'}),timeout=a.timeout) as src,tmp.open('wb') as f:
     while True:
      x=src.read(1<<20)
      if not x:break
      f.write(x)
    with zipfile.ZipFile(tmp) as z:
     if z.testzip() is not None:raise ValueError('XLSX CRC failure')
     if 'xl/workbook.xml' not in z.namelist():raise ValueError('Not an XLSX workbook')
    tmp.replace(dst)
   with zipfile.ZipFile(dst) as z:
    if z.testzip() is not None:raise ValueError('Existing workbook corrupt')
   item.update(status='DOWNLOADED',size=dst.stat().st_size,sha256=hashlib.sha256(dst.read_bytes()).hexdigest())
  except Exception as e:item.update(status='NOT_DOWNLOADED',error=repr(e))
  receipt.append(item);print(name,item['status'],flush=True)
 (a.out/'download_receipts.json').write_text(json.dumps(receipt,indent=2)+'\n')
 if any(r['status']!='DOWNLOADED' for r in receipt):sys.exit(2)
if __name__=='__main__':main()
