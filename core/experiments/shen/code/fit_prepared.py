#!/usr/bin/env python3
"""Fit unchanged v31 models, then freeze ALL predictions before any formal score.

No evaluation truth files are read. Source qualification and public preregistration
remain investigator actions; this runner does not certify either of them.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='1'
from pathlib import Path
import argparse,json,hashlib,datetime,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'vendor'))
import witnesscell_growth as wc

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--stage',choices=['smoke','formal'],required=True)
    p.add_argument('--max-seconds-per-cell',type=float,default=600)
    a=p.parse_args()
    if a.out.exists(): p.error('Fresh output path required; never overwrite frozen predictions')
    receipt=json.loads((a.prepared/'PREPARATION_RECEIPT.json').read_text())
    cells=['HeLa'] if a.stage=='smoke' else ['HeLa','A549','293T']
    relative='smoke_input' if a.stage=='smoke' else 'training_inputs'
    a.out.mkdir(parents=True)
    entries={}
    for cell in cells:
        source=a.prepared/relative/f'{cell}.json'
        inp=json.loads(source.read_text())
        meta=wc.predict(inp,a.out/cell,max_seconds=a.max_seconds_per_cell)
        target=a.out/cell/'results/high/predictions.npz'
        entries[cell]={'input_sha256':sha(source),'predictions_sha256':sha(target),
                       'metadata_sha256':sha(a.out/cell/'fit_metadata.json'),
                       'prediction_file':str(target.relative_to(a.out)),
                       'query_outcomes_read':False,'selections':meta['selections']}
    freeze={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':a.stage,
            'protocol_sha256':sha(ROOT/'protocol/PROTOCOL.json'),
            'prepared_receipt_sha256':sha(a.prepared/'PREPARATION_RECEIPT.json'),
            'source_declared':receipt['provenance']['source_name'],
            'all_predictions_complete':True,'query_outcomes_read':False,
            'models':entries,'evaluation_performed':False}
    (a.out/'PREDICTIONS_FROZEN.json').write_text(json.dumps(freeze,indent=2)+'\n')
    print(json.dumps({'status':'PREDICTIONS_FROZEN_NOT_SCORED','stage':a.stage,'contexts':cells}))
if __name__=='__main__':main()
