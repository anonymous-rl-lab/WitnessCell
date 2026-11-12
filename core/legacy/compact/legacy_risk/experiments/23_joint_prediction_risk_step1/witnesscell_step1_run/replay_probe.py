import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
import argparse,hashlib,json,pickle,sys,time,resource
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--out',required=True);p.add_argument('--dataset',default='Norman');p.add_argument('--seed',type=int,default=1)
a=p.parse_args();repo=Path(a.repo);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(repo/'experiments/19_v14_incremental_amplitude_gate/src'))
import witnesscell_v4_core as core
t=time.monotonic()
ledger=json.loads((repo/f'experiments/22_metric_calibration_stress/assets/formal_contracts/{a.dataset}.split_ledger.json').read_text())['split_contract'][str(a.seed)]
with np.load(repo/f'experiments/18_dual_head_evidence_gate/cache/{a.dataset}.condition_moments.npz',allow_pickle=False) as z:
    names=z['conditions'].astype(str).tolist();idx={c:i for i,c in enumerate(names)};genes=z['genes'].astype(str).tolist()
    allowed=sorted(set(['control']+ledger['train']+ledger['validation']))
    # Only these slices are supplied to fit; no test truth is inspected/scored.
    means={c:z['means'][idx[c]].astype(float) for c in allowed}
    variances={c:z['variances'][idx[c]].astype(float) for c in allowed}
    counts={c:int(z['counts'][idx[c]]) for c in allowed}
with (repo/'experiments/15_scperturbench_sota/module/data/gears_assets/gene2go_all.pkl').open('rb') as f:go=pickle.load(f)
assert not set(ledger['test'])&set(means)
print(json.dumps({'event':'fit_start','dataset':a.dataset,'seed':a.seed,'training_singles':sum(c!='control' and '+' not in c for c in ledger['train']),'training_doubles':sum('+' in c for c in ledger['train']),'time':time.time()}),flush=True)
fit_start=time.monotonic()
fit=core.fit_predict(means,variances,counts,genes,ledger['train'],ledger['validation'],ledger['test'],identity_mode='multihead_v14',gene2go=go)
fit_sec=time.monotonic()-fit_start
with np.load(repo/f'experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo/{a.dataset}/seed{a.seed}/deploy_predictions.npz',allow_pickle=False) as z:
    targets=z['conditions'].astype(str).tolist();actual=np.stack([fit.predictions[c] for c in targets]).astype(np.float32);expected=z['prediction']
    assert np.array_equal(z['genes'].astype(str),genes)
    diff=float(np.max(np.abs(actual-expected)));match=bool(np.allclose(actual,expected,rtol=1e-6,atol=1e-6));exact=bool(np.array_equal(actual,expected))
result={'dataset':a.dataset,'seed':a.seed,'fit_seconds':fit_sec,'total_seconds':time.monotonic()-t,'max_abs_error':diff,'within_tolerance':match,'float32_exact':exact,'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'alpha':fit.alpha,'noise_ratio':fit.noise_ratio,'gamma':fit.gamma,'test_truth_evaluated':False}
(out/'replay_result.json').write_text(json.dumps(result,indent=2)+'\n')
if match:
    with (out/'replay_fit.pkl').open('wb') as f:pickle.dump(fit,f)
print(json.dumps(result),flush=True)
if not match:raise SystemExit('REPLAY MISMATCH')
