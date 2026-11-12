"""CUDA placement adapter; the four frozen training source files stay unchanged.

The network forward runs on CUDA in FP32, and its output returns to CPU with
autograd intact. The existing CPU loss, batch order, scheduler, selection,
checkpoint cadence and training budgets are therefore reused without rewriting.
No AMP or TF32 is used. Checkpoints are serialized on CPU for portable recovery.
"""
import copy, hashlib, json, os, types
from pathlib import Path
import torch

R = Path(__file__).resolve().parent
ADAPTER_SHA = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

def environment():
    import numpy, torch_geometric
    return dict(adapter_sha256=ADAPTER_SHA, device='cuda:0', gpu=torch.cuda.get_device_name(0),
        torch=torch.__version__, cuda=torch.version.cuda, numpy=numpy.__version__,
        torch_geometric=torch_geometric.__version__, fp32=True, amp=False, tf32=False,
        forward_device='cuda', loss_device='cpu', training_rules_changed=False,
        initialization='frozen CPU initial state for each seed')

def on_cuda(model):
    model.to('cuda:0')
    for name in ('G_coexpress','G_coexpress_weight','G_sim','G_sim_weight'):
        setattr(model,name,getattr(model,name).to('cuda:0'))
    model.args = dict(model.args, device='cuda:0')
    forward = model.forward
    def bridge(self, data):
        obj = types.SimpleNamespace(x=data.x.to('cuda:0'), pert_idx=data.pert_idx,
                                    batch=data.batch.to('cuda:0'))
        return forward(obj).cpu()
    model.forward = types.MethodType(bridge,model)
    return model

def cpu_tree(obj):
    if isinstance(obj,torch.Tensor): return obj.detach().cpu()
    if isinstance(obj,dict):
        result=copy.copy(obj)
        for k,v in obj.items(): result[k]=cpu_tree(v)
        return result
    if isinstance(obj,list): return [cpu_tree(x) for x in obj]
    if isinstance(obj,tuple): return tuple(cpu_tree(x) for x in obj)
    return obj

def install():
    assert torch.cuda.is_available(), 'CUDA is required; do not silently use CPU'
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    torch.backends.cudnn.benchmark=False
    import train_cell_models_runtime as t
    old_new=t.Data.newmodel
    old_save=torch.save
    old_load=torch.load
    old_write=t.writej
    meta=environment()
    def newmodel(self):
        assert self.kind=='gears', 'Only remaining GEARS jobs are approved here'
        model=old_new(self)
        seed=int(torch.initial_seed())
        initial=R/'initial_states'/f'seed{seed}.pt'
        assert initial.exists(), f'Missing frozen initialization for seed {seed}'
        model.load_state_dict(old_load(initial,map_location='cpu',weights_only=True))
        return on_cuda(model)
    def save(obj,path,*args,**kwargs):
        value=cpu_tree(obj)
        if isinstance(value,dict) and 'optimizer' in value and 'runtime_fingerprint' in value:
            value['execution_fingerprint']=meta
            value['cuda_rng']=torch.cuda.get_rng_state_all()
        return old_save(value,path,*args,**kwargs)
    def load(path,*args,**kwargs):
        kwargs.setdefault('map_location','cpu')
        value=old_load(path,*args,**kwargs)
        if isinstance(value,dict) and 'execution_fingerprint' in value:
            assert value['execution_fingerprint']['adapter_sha256']==ADAPTER_SHA, 'CUDA adapter changed during resume'
            if 'cuda_rng' in value: torch.cuda.set_rng_state_all(value['cuda_rng'])
        return value
    def write(path,value):
        if Path(path).name=='FREEZE.json': value=dict(value,execution_fingerprint=meta)
        return old_write(path,value)
    t.Data.newmodel=newmodel;torch.save=save;torch.load=load;t.writej=write
    return t,meta

if __name__=='__main__':
    t,meta=install()
    print(json.dumps(dict(event='CUDA_EXECUTION',**meta)),flush=True)
    import train_worker_20261106 as w
    raise SystemExit(w.main())
