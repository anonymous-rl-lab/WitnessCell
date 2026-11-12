"""Replay all frozen neural predictions using the untouched upstream models."""
from train_cell_models import *
from lazy_cell_data import LazyData
import gc

def main():
 torch.set_num_threads(2);rows=[]
 root=ROOT/'results/detected49'
 for kind in ('gears','scouter'):
  for seed in range(3):
   for f in range(5):
    out=root/f'{kind}_seed{seed}'/f'fold_{f}'
    assert (out/'FREEZE.json').exists(),f'Incomplete model: {out}'
    freeze=json.loads((out/'FREEZE.json').read_text());assert sha(out/'predictions.npz')==freeze['predictions_sha256']
    d=LazyData(f,kind);model=d.newmodel()
    # No runtime forward substitutions: native upstream class and parameters.
    assert type(model) is (GM.GEARS_Model if kind=='gears' else SM.ScouterModel)
    model.load_state_dict(torch.load(out/'refit/last.pt',weights_only=True));model.eval()
    replay=predict(model,d,d.q);stored=load(out/'predictions.npz')
    assert np.array_equal(stored['query'],d.q);assert np.array_equal(stored['genes'],d.m['genes'])
    diff=replay-stored['pred_effect'];maximum=float(np.max(abs(diff)));rms=float(np.sqrt(np.mean(diff**2)))
    np.testing.assert_allclose(replay,stored['pred_effect'],rtol=1e-4,atol=5e-5)
    row=dict(model=kind,seed=seed,fold=f,status='PASS',maximum_absolute_difference=maximum,rms_difference=rms,native_source=True,weight_sha256=sha(out/'refit/last.pt'))
    rows.append(row);print(json.dumps(row),flush=True)
    del d,model,replay,stored;gc.collect()
 result=dict(status='PASS',models_replayed=30,original_upstream_forward=True,tolerance=dict(rtol=1e-4,atol=5e-5),rows=rows)
 writej(root/'evaluation_main_20261105/NATIVE_WEIGHT_REPLAY.json',result)
if __name__=='__main__':main()
