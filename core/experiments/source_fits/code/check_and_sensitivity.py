"""Check recovered feature/weight replay and report conditional sensitivities."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import json
import numpy as np
import torch
from fit_models import ROOT,load_input,Projection,dump
import sf12_legacy_shared_core as sf12
from partner_model_v4 import PartnerHeadV4,PartnerConfig
torch.set_num_threads(1)
checks=[];sens=[]
for cell in ['HeLa','A549','293T','Jurkat']:
    d=load_input(ROOT/'inputs'/f'{cell}.json');p=ROOT/'results_phase2'/cell
    z=np.load(p/'PREDICTIONS.npz');knowledge=np.load(ROOT/'inputs'/f'{cell}_GenePT64.npz')['embedding']
    selected=json.loads((p/'SELECTION_BEFORE_OUTER.json').read_text())['selected']
    # Reconstruct the actual frozen V4 network predictions from weights.
    kk=knowledge.astype(np.float32);kk-=kk.mean(0);kk/=np.maximum(np.sqrt(np.mean(kk*kk,axis=1,keepdims=True)),1e-8)
    kt=torch.tensor(kk);s=d['singles'].mean(1);s/=max(np.sqrt(np.mean(s*s)),1e-12);st=torch.tensor(s,dtype=torch.float32)
    q=torch.tensor(d['query'],dtype=torch.long);row=selected['network_selected'];delta=np.zeros(len(q));chunk_error=0.
    if row['key']!='M0':
        for seed in [17,29]:
            state=torch.load(p/'final'/f"{row['key']}_seed{seed}.pt",weights_only=False)
            model=PartnerHeadV4(PartnerConfig(**state['config']));model.load_state_dict(state['state_dict']);model.eval()
            def forward(e):
                a,b=e[:,0],e[:,1]
                with torch.no_grad():return model(kt[a],kt[b],torch.zeros((1,64)),st[a],st[b],torch.zeros(())).numpy()*state['scale']
            raw=forward(q);chunk=np.concatenate([forward(t) for t in q.split(53)])
            chunk_error=max(chunk_error,float(np.max(abs(raw-chunk))))
            shift=state['endpoint_shift'];delta+=(raw-shift[d['query'][:,0]]-shift[d['query'][:,1]])/2
        error=float(np.max(abs(z['network_selected']-(z['M0']+row['beta']*delta[:,None]))));assert error<1e-7
    else:error=0.
    # Original dense SF12 fitter vs algebraically equivalent efficient adapter.
    idx=np.arange(min(150,len(d['edges'])));e=d['edges'][idx];P=Projection(e,len(d['genes']))
    ft,fq=sf12.features(d['singles'].mean(1)[:,None],knowledge,e,d['query'][:10]);r=d['y'][idx].mean(1)[:,None]
    state=sf12.fit_fixed_m0(r,P.X,ft,rho=.1,amplitude=1.)
    a=P.apply(ft[:,0,:]);gram=a.T@a;rhs=a.T@P.apply(r[:,0]);sc=np.sqrt(np.maximum(np.diag(gram),0));sc=np.maximum(sc,max(sc.max()*1e-10,1e-15))
    theta=np.linalg.solve(gram/sc[:,None]/sc[None,:]+.1*np.eye(12),rhs/sc)/sc
    pred=fq[:,0,:]@theta-P.extend(ft[:,0,:]@theta,d['query'][:10]);ref=sf12.predict(state,np.zeros((10,1)),Projection(d['query'][:10],len(d['genes'])).X,fq)[1][:,0]
    eq=float(np.max(abs(pred-ref)));assert eq<1e-10
    checks.append(dict(cell=cell,network_checkpoint_replay_max_abs=error,network_chunk_max_abs=chunk_error,sf12_original_adapter_max_abs=eq))
    z=np.load(ROOT/'evaluation'/f'{cell}_scored_predictions.npz');edges=z['edges'];P=Projection(edges,len(d['genes']));y=z['y'];base=z['M0']
    eb=P.apply(y-base);eo=P.apply(y-z['original_rank1']);en=P.apply(y-z['combined_selected'])
    b=(eb*eb).sum(1);o=(eo*eo).sum(1);a=(en*en).sum(1);rng=np.random.default_rng(20261108);boot=[]
    for _ in range(5000):
        ix=rng.integers(0,len(a),len(a));boot.append([100*(b[ix].sum()-a[ix].sum())/b[ix].sum(),100*(o[ix].sum()-a[ix].sum())/b[ix].sum()])
    jack=[]
    for g in range(len(d['genes'])):
        ix=np.flatnonzero(~np.any(edges==g,axis=1));pp=Projection(edges[ix],len(d['genes']));pb=pp.apply(y[ix]-base[ix]);pn=pp.apply(y[ix]-z['combined_selected'][ix]);jack.append(100*(1-np.sum(pn*pn)/np.sum(pb*pb)))
    sens.append(dict(cell=cell,selected_partner_gain_ci=np.quantile(boot,[.025,.975],axis=0)[:,0],
        improvement_over_original_percentage_points_ci=np.quantile(boot,[.025,.975],axis=0)[:,1],
        leave_one_endpoint_gain_range=[min(jack),max(jack)],leave_one_endpoint_positive=sum(np.array(jack)>0),endpoints=len(jack),
        caveat='descriptive 5000 pair-block resamples of fixed projected errors; shared endpoints, fitting and model search uncertainty not covered; delete-endpoint projections recomputed'))
dump(ROOT/'diagnostics/RECOVERED_MODEL_CHECKS.json',checks);dump(ROOT/'evaluation/SENSITIVITY.json',sens)
print(json.dumps({'checks':checks,'sensitivity':sens},default=lambda x:x.tolist(),indent=2))
