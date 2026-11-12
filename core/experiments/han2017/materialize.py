"""Restore deterministic Han duplicate arrays in a fresh directory; no fitting."""
from pathlib import Path
import argparse,hashlib,json,shutil,sys
import numpy as np
R=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--mode',choices=['exact','compatible'],default='exact');a=ap.parse_args();out=a.out.resolve()
 if out.exists() or out.is_relative_to(R):ap.error('Output must be new and outside the compact experiment.')
 ledger=json.loads((R/'COMPACTION_LEDGER.json').read_text());recipes=json.loads((R/'RECIPES.json').read_text())['derived_files']
 for x in ledger['retained_files']:assert sha(R/'snapshot'/x['path'])==x['sha256'],x['path']
 shutil.copytree(R/'snapshot',out);sys.path.insert(0,str(out/'code'));import run,common
 status=[]
 for panel in run.PANELS:
  data=common.load(out/f'data/{panel}_frozen_predictions.npz');truth=common.load(out/f'data/{panel}_outcomes.npz')['outcomes']
  for seed in range(20261108,20261113):
   for fold in range(3):
    d=out/'data/splits'/panel/str(seed)/str(fold);ident=json.loads((d/'identities.json').read_text());q=np.asarray(ident['query']);tr=np.asarray(ident['train'])
    np.savez_compressed(d/'calibration.npz',ids=tr,y=truth[tr]);np.savez_compressed(d/'query_outcomes.npz',ids=q,y=truth[q])
    d=out/'results'/panel/str(seed)/str(fold);state=json.loads((d/'PARAMETERS.json').read_text());_,arms=run.arms(data,q,state,panel)
    save=dict(ids=q,target_control=arms['target_control'])
    if 'anchored_target_control' in arms:save['anchored_target_control']=arms['anchored_target_control']
    np.savez_compressed(d/'QUERY_PREDICTIONS.npz',**save)
 for x in recipes:
  p=out/x['path'];same=sha(p)==x['sha256'];status.append(dict(path=x['path'],exact_file_bytes=same))
  if a.mode=='exact':assert same,('A numerical or NPZ-container environment difference prevents byte-exact reconstruction',x['path'])
 report=dict(status='PASS',mode=a.mode,regenerated_files=len(status),exact_file_matches=sum(x['exact_file_bytes'] for x in status),all_derived_files_byte_exact=all(x['exact_file_bytes'] for x in status),base_model_fits=0,files=status,scope='Stored scientific inputs, coefficients and reference scores remain byte-identical. Exact mode additionally requires original bytes for all 180 deterministic duplicates. Compatible mode delegates fresh inference and scoring acceptance to the stagewise replay comparator; this materialization alone does not establish numerical equivalence.')
 (out/'COMPACT_RESTORE_RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='files'}),flush=True)
if __name__=='__main__':main()
