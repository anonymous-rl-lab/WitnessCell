"""One actual train-only subsample: full v14 vs exact non-target projection."""
import json,time
import numpy as np
import run_step1 as r
s=r.Split('Norman',1)
scheme,fold,train,val,queries=next(s.subsets())
rng=np.random.default_rng(int(r.h(r.PID,s.ds,s.seed,scheme,fold,0)[:8],16))
selected={'control'}
for cs in [[c for c in train if c!='control' and '+' not in c],[c for c in train if '+' in c]]:
    if cs:selected.update(rng.choice(cs,size=max(1,int(.8*len(cs))),replace=False).tolist())
tr=[c for c in train if c in selected]
full,fo=s.fit(tr,val,queries,'verification_bag_full',want_object=True,projected=False)
projection,po=s.fit(tr,val,queries,'verification_bag_projection',want_object=True,projected=True)
mask=np.stack([s.mask(c) for c in queries])
result={'dataset':s.ds,'seed':s.seed,'scenario':'pair_fold0_bag0','queries':len(queries),'training_singles':sum(c!='control' and '+' not in c for c in tr),'full_gate_active':bool(fo.identity_head.correction_active),'non_target_max_abs_error':float(np.max(np.abs(full[mask]-projection[mask]))),'non_target_float32_exact':bool(np.array_equal(full[mask],projection[mask])),'parameters_equal':bool((fo.alpha,fo.noise_ratio,fo.gamma)==(po.alpha,po.noise_ratio,po.gamma)),'test_labels_read':False}
assert result['non_target_float32_exact'] and result['parameters_equal']
r.dump(r.HERE/'projection_subsample_equivalence.json',result)
print(json.dumps(result),flush=True)
