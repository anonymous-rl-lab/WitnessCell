#!/usr/bin/env python3
"""Score a fully frozen run; refuse missing/changed predictions or misaligned truth.

This program evaluates points, repeat-specific results and fixed-projector
endpoint-reweighting sensitivity. It does not provide population-level inference.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,json,hashlib,datetime
import numpy as np
from shen_adapter import score
from endpoint_sensitivity import endpoint_sensitivity
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--fits',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists():p.error('Fresh score directory required')
    f=json.loads((a.fits/'PREDICTIONS_FROZEN.json').read_text());stage=f['stage']
    cells=['HeLa'] if stage=='smoke' else ['HeLa','A549','293T']
    if not f['all_predictions_complete'] or set(f['models'])!=set(cells):p.error('Incomplete frozen prediction set')
    if sha(a.prepared/'PREPARATION_RECEIPT.json')!=f['prepared_receipt_sha256']:p.error('Changed prepared inputs')
    if sha(ROOT/'protocol/PROTOCOL.json')!=f['protocol_sha256']:p.error('Protocol differs from fitted run')
    # Verify ALL files before opening any evaluation outcomes.
    inputs={}
    for cell in cells:
        source=a.prepared/('smoke_input' if stage=='smoke' else 'training_inputs')/f'{cell}.json'
        item=f['models'][cell]
        if sha(source)!=item['input_sha256']:p.error('Fitted input mismatch')
        if sha(a.fits/item['prediction_file'])!=item['predictions_sha256']:p.error('Prediction hash mismatch')
        if sha(a.fits/cell/'fit_metadata.json')!=item['metadata_sha256']:p.error('Metadata hash mismatch')
        inputs[cell]=json.loads(source.read_text())
    out={};sensitivity_inputs={}
    for cell in cells:
        inp=inputs[cell]
        truthpath=a.prepared/('smoke_input/HeLa_truth.npz' if stage=='smoke' else f'sealed_evaluation/{cell}.npz')
        if stage=='formal':
            r=json.loads((a.prepared/'PREPARATION_RECEIPT.json').read_text())
            if sha(truthpath)!=r['summary'][cell]['sealed_truth_sha256']:p.error('Truth hash changed')
        with np.load(truthpath,allow_pickle=False) as truth, np.load(a.fits/f['models'][cell]['prediction_file'],allow_pickle=False) as z:
            if truth['genes'].tolist()!=inp['genes'] or truth['pairs'].tolist()!=inp['queries']:p.error('Truth identities misaligned')
            qids=['|'.join(sorted(q)) for q in inp['queries']]
            if z['global_ids'].tolist()!=qids:p.error('Prediction identities misaligned')
            y=truth['y'];base=z['M0'];lookup={g:i for i,g in enumerate(inp['genes'])};s=np.array(inp['singles'])
            add=np.array([s[lookup[g]]+s[lookup[h]] for g,h in inp['queries']])
            sensitivity_inputs[cell]={'genes':inp['genes'],'pairs':inp['queries'],'y':y.copy(),'base':base.copy(),'pred':z['rank1_selected'].copy()}
            out[cell]={name:score(inp['genes'],inp['queries'],y,base,pred)
                       for name,pred in [('M0',base),('rank1_primary',z['rank1_selected']),
                                         ('spectral_secondary',z['spectral_selected']),('single_additive',add)]}
    gains=[out[c]['rank1_primary']['partner_gain'] for c in cells]
    report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':stage,
            'source_declared':f['source_declared'],'models':out,
            'mean_within_context_partner_gain':None if any(x is None for x in gains) else float(np.mean(gains)),
            'definition':'Equal context average of relative gains; not pooled cross-context SSE',
            'conditional_sensitivity':endpoint_sensitivity(sensitivity_inputs),
            'uncertainty_computed':'conditional sensitivity only','formal_score':stage=='formal',
            'no_sign_count_success_gate':True}
    a.out.mkdir(parents=True);(a.out/'SCORES.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'SCORED','stage':stage,'source':f['source_declared']}))
if __name__=='__main__':main()
