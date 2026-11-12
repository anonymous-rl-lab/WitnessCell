import json,numpy as np,torch
import train_cell_models_runtime as t
from gpu_adapter_20261107 import on_cuda
torch.set_num_threads(2)
torch.backends.cuda.matmul.allow_tf32=False
torch.backends.cudnn.allow_tf32=False
d=t.Data(2,'gears');m=d.newmodel()
m.load_state_dict(torch.load(t.ROOT/'validation_reference/refit/last.pt',map_location='cpu',weights_only=True))
saved=t.load(t.ROOT/'validation_reference/predictions.npz')
assert np.array_equal(d.q,saved['query'])
cpu=t.predict(m,d,d.q);gpu=t.predict(on_cuda(m),d,d.q)
np.testing.assert_allclose(cpu,saved['pred_effect'],rtol=1e-4,atol=5e-5)
np.testing.assert_allclose(gpu,saved['pred_effect'],rtol=1e-4,atol=5e-5)
proof=dict(status='PASS',cpu_reference_max=float(abs(cpu-saved['pred_effect']).max()),cuda_reference_max=float(abs(gpu-saved['pred_effect']).max()),outer_truth_used=False,original_weight_sha256=t.sha(t.ROOT/'validation_reference/refit/last.pt'))
t.writej(t.ROOT/'CUDA_REFERENCE_REPLAY.json',proof);print(json.dumps(proof),flush=True)
