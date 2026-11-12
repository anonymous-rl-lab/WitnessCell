"""Held-in numerical and speed gate. No outer truth or outer metric is used."""
import copy, json, time, statistics
import numpy as np
import torch
import train_cell_models_runtime as t
from gpu_adapter_20261107 import on_cuda, environment

torch.set_num_threads(2)
torch.backends.cuda.matmul.allow_tf32=False
torch.backends.cudnn.allow_tf32=False
torch.backends.cudnn.benchmark=False
assert torch.cuda.is_available()
d=t.Data(3,'gears')
Y,labels,ctrl,pi,nz=d.train_arrays(d.fit,2)
torch.manual_seed(2)
cpu=d.newmodel();cpu.load_state_dict(torch.load(t.ROOT/'initial_states/seed2.pt',weights_only=True))
gpu=on_cuda(copy.deepcopy(cpu))
order=np.random.default_rng(2).permutation(len(Y))
timings={'cpu':[],'cuda':[]};losses={'cpu':[],'cuda':[]}
opts={key:torch.optim.Adam(model.parameters(),lr=.001,weight_decay=.0005)
      for key,model in [('cpu',cpu),('cuda',gpu)]}
single={}
for i in range(12):
    ix=order[32*i:32*(i+1)];x=torch.from_numpy(ctrl[ix]);y=torch.from_numpy(Y[ix]);indices=[pi[j] for j in ix]
    predictions={};gradients={}
    for key,model in [('cpu',cpu),('cuda',gpu)]:
        model.train();start=time.perf_counter();opt=opts[key];opt.zero_grad(set_to_none=True)
        pred=model(t.batch(x,indices))
        loss=t.GU.loss_fct(pred,y,np.where(labels[ix]=='control','ctrl',labels[ix]),torch.tensor(d.ctrlmean,dtype=torch.float32),.1,nz)
        loss.backward()
        assert torch.isfinite(loss) and all(torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None)
        predictions[key]=pred.detach().numpy()
        if i==0: gradients[key]=torch.cat([p.grad.detach().cpu().reshape(-1) for p in model.parameters() if p.grad is not None]).numpy()
        torch.nn.utils.clip_grad_value_(model.parameters(),1.);opt.step()
        if key=='cuda': torch.cuda.synchronize()
        timings[key].append(time.perf_counter()-start);losses[key].append(float(loss.detach()))
    if i==0:
        a,b=predictions['cpu'],predictions['cuda'];ga,gb=gradients['cpu'],gradients['cuda']
        single=dict(prediction_max=float(abs(a-b).max()),prediction_relative_l2=float(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-30)),
                    gradient_relative_l2=float(np.linalg.norm(ga-gb)/max(np.linalg.norm(ga),1e-30)))
    if i in (0,5,11): print(json.dumps(dict(step=i+1,seconds={k:timings[k][-1] for k in timings})),flush=True)
inner_cpu=t.predict(cpu,d,d.val,True);inner_gpu=t.predict(gpu,d,d.val,True)
relative=float(np.linalg.norm(inner_cpu-inner_gpu)/max(np.linalg.norm(inner_cpu),1e-30))
medians={k:statistics.median(v[2:]) for k,v in timings.items()}
result=dict(status='PASS' if single['prediction_relative_l2']<1e-4 and single['gradient_relative_l2']<1e-3 and relative<.05 and medians['cuda']<medians['cpu'] else 'FAIL',
    environment=environment(),single_step=single,inner_prediction_relative_l2_after_12_steps=relative,
    median_seconds=medians,speedup=medians['cpu']/medians['cuda'],peak_gpu_memory_bytes=torch.cuda.max_memory_allocated(),
    durations=timings,losses=losses,outer_outcomes_used=False)
t.writej(t.ROOT/'CUDA_PROBE.json',result)
print(json.dumps({k:v for k,v in result.items() if k not in ('durations','losses')}),flush=True)
assert result['status']=='PASS'
