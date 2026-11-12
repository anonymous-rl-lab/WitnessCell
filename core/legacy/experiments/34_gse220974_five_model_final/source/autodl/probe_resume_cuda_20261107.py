import json,torch,numpy as np
from gpu_adapter_20261107 import install
t,meta=install();torch.set_num_threads(2);rows=[]
for seed,phase in [(0,'select_schedule5'),(1,'select_schedule1')]:
    torch.manual_seed(seed);d=t.Data(3,'gears');model=d.newmodel()
    ck=torch.load(t.ROOT/'results/detected49'/('gears_seed%d'%seed)/'fold_3'/phase/'resume.pt',weights_only=False)
    model.load_state_dict(ck['model']);opt=torch.optim.Adam(model.parameters(),lr=.001,weight_decay=.0005);opt.load_state_dict(ck['optimizer'])
    Y,labels,ctrl,pi,nz=d.train_arrays(d.fit,seed);p=ck['partial'];ix=p['order'][p['next_start']:p['next_start']+32]
    model.train();opt.zero_grad(set_to_none=True);pred=model(t.batch(torch.from_numpy(ctrl[ix]),[pi[i] for i in ix]))
    loss=t.GU.loss_fct(pred,torch.from_numpy(Y[ix]),np.where(labels[ix]=='control','ctrl',labels[ix]),torch.tensor(d.ctrlmean,dtype=torch.float32),.1,nz)
    loss.backward();assert torch.isfinite(loss) and all(torch.isfinite(v.grad).all() for v in model.parameters() if v.grad is not None)
    torch.nn.utils.clip_grad_value_(model.parameters(),1.);opt.step()
    rows.append(dict(seed=seed,phase=phase,epoch_completed=ck['epoch'],next_start=p['next_start'],loss=float(loss.detach()),model_optimizer_restored=True))
    del model,opt,d,ck;torch.cuda.empty_cache()
proof=dict(status='PASS',rows=rows,formal_checkpoint_modified=False,outer_outcomes_used=False);t.writej(t.ROOT/'CUDA_RESUME_PROBE.json',proof);print(json.dumps(proof),flush=True)
