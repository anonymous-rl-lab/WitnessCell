"""Functional adapters of recovered SF12 and actual trainable PartnerHeadV4.

Uses training responses only. RNA source weights are not applied to fitness.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[k]='1'
import argparse, copy, json, time
from pathlib import Path
import numpy as np
import torch
from fit_models import ROOT, Projection, load_input, dump, sha, native, stats, solve_beta
import sf12_legacy_shared_core as sf12
from partner_model_v4 import PartnerConfig, PartnerHeadV4
torch.set_num_threads(1)


def one_fit(d,idx,query,knowledge,plan,outpath):
    n=len(d['genes']);edges=d['edges'][idx];y=d['y'][idx]
    _,predict=native.baseline(d['singles'],edges,y,time.monotonic()+1800)
    base=predict(query);resid=y-predict(edges);P=Projection(edges,n)
    out={};audit={};states={}
    # Original SF12 feature map, original normalized ridge, new scalar readout.
    ft,fq=sf12.features(d['singles'].mean(1)[:,None],knowledge,edges,query)
    ft=ft[:,0,:];fq=fq[:,0,:];pt=P.apply(ft)
    gram=pt.T@pt;rhs=pt.T@P.apply(resid.mean(1))
    sc=np.sqrt(np.maximum(np.diag(gram),0));sc=np.maximum(sc,max(sc.max()*1e-10,1e-15))
    normal=gram/sc[:,None]/sc[None,:]
    for rho in plan['sf12_rho']:
        theta=np.linalg.solve(normal+rho*np.eye(12),rhs/sc)/sc
        delta=fq@theta-P.extend(ft@theta,query)
        key=f'sf12_r{rho:g}';out[key]=np.repeat(delta[:,None],2,axis=1)
        states[key+'__theta']=theta;states[key+'__shift']=P.coeff(ft@theta)
    np.savez_compressed(outpath/'sf12_parameters.npz',**states)
    # Reuse the actual V4 head class. Fitness has no RNA readout identity, so q=0.
    kk=knowledge.astype(np.float32);kk-=kk.mean(0)
    kk/=np.maximum(np.sqrt(np.mean(kk*kk,axis=1,keepdims=True)),1e-8)
    kt=torch.tensor(kk);s=d['singles'].mean(1);s=s/max(np.sqrt(np.mean(s*s)),1e-12)
    st=torch.tensor(s,dtype=torch.float32);all_edges=np.vstack([edges,query])
    ei=torch.tensor(all_edges,dtype=torch.long);m=len(edges)
    xx=torch.tensor(P.X,dtype=torch.float32);left=torch.tensor(P.inv@P.X.T,dtype=torch.float32)
    target=P.apply(resid.mean(1));scale=max(np.sqrt(np.mean(target*target)),1e-10)
    yt=torch.tensor(target/scale,dtype=torch.float32)
    sums={};histories=[]
    for wd in plan['weight_decay']:
        for seed in plan['network_seeds']:
            torch.manual_seed(seed)
            model=PartnerHeadV4(PartnerConfig(use_extra=False))
            optimizer=torch.optim.AdamW(model.parameters(),lr=plan['network_lr'],weight_decay=wd)
            def forward(e):
                a,b=e[:,0],e[:,1]
                return model(kt[a],kt[b],torch.zeros((1,64)),st[a],st[b],torch.zeros(()))
            history=[]
            for step in range(1,max(plan['checkpoints'])+1):
                model.train();raw=forward(ei[:m]);pred=raw-xx@(left@raw)
                loss=torch.mean((pred-yt)**2)
                optimizer.zero_grad();loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
                if step in plan['checkpoints']:
                    model.eval()
                    with torch.no_grad():
                        raw=forward(ei).numpy()*scale
                        delta=raw[m:]-P.extend(raw[:m],query)
                        rev=forward(ei[m:][:,[1,0]]).numpy()*scale
                    assert np.max(abs(raw[m:]-rev))<1e-6*max(1,scale)
                    key=f'v4_w{wd:g}_s{step}'
                    sums.setdefault(key,np.zeros(len(query)));sums[key]+=delta/len(plan['network_seeds'])
                    torch.save(dict(state_dict=model.state_dict(),config=vars(model.config),scale=scale,
                        endpoint_shift=P.coeff(raw[:m]),step=step,seed=seed,weight_decay=wd),
                        outpath/f'{key}_seed{seed}.pt')
                    history.append(dict(step=step,loss=float(loss.detach()),symmetry_max=float(np.max(abs(raw[m:]-rev)))))
            histories.append(dict(seed=seed,weight_decay=wd,history=history))
    for key,v in sums.items():out[key]=np.repeat(v[:,None],2,axis=1)
    audit={'histories':histories,'output_scale':scale,'network_parameters':sum(x.numel() for x in model.parameters()),
           'new_functional_adapter':True,'source_RNA_weights_used':False,'sf12_feature_rank':int(np.linalg.matrix_rank(normal))}
    return base,out,audit


def run(cell):
    plan=json.loads((ROOT/'PLAN_PHASE2.json').read_text());out=ROOT/'results_phase2'/cell
    out.mkdir(parents=True,exist_ok=True)
    if (out/'FREEZE.json').exists():print(cell,'already frozen');return
    d=load_input(ROOT/'inputs'/f'{cell}.json');knowledge=np.load(ROOT/'inputs'/f'{cell}_GenePT64.npz')['embedding']
    acc={};base_sse=0.;base_partner=0.;begin=time.monotonic()
    for f in range(3):
        sub=out/f'inner_{f}';sub.mkdir(exist_ok=True)
        tr=np.flatnonzero(d['fold']!=f);va=np.flatnonzero(d['fold']==f);q=d['edges'][va];y=d['y'][va]
        base,cs,audit=one_fit(d,tr,q,knowledge,plan,sub);P=Projection(q,len(d['genes']))
        base_sse+=np.sum((base-y)**2);base_partner+=np.sum(P.apply(base-y)**2)
        np.savez_compressed(sub/'predictions.npz',edges=q,ids=d['ids'][va],y=y,M0=base,**cs)
        dump(sub/'audit.json',audit)
        for key,c in cs.items():
            acc.setdefault(key,np.zeros(6));acc[key]+=stats(base,y,np.zeros_like(base),c,P)
        print(json.dumps(dict(cell=cell,phase='inner_complete',fold=f,seconds=time.monotonic()-begin)),flush=True)
    rows=[dict(key=k,**solve_beta(v,float(base_sse))) for k,v in acc.items()]
    zero=dict(key='M0',beta=0.,eligible=True,partner_sse=float(base_partner),rna_sse=float(base_sse))
    groups={}
    for group,prefix in [('selected',''),('sf12_selected','sf12'),('network_selected','v4')]:
        groups[group]=min([zero]+[r for r in rows if r['eligible'] and r['key'].startswith(prefix)],key=lambda r:(r['partner_sse'],r['rna_sse'],r['key']))
    dump(out/'SELECTION_BEFORE_OUTER.json',dict(cell=cell,baseline_sse=base_sse,baseline_partner_sse=base_partner,selected=groups,candidates=rows,outer_truth_read=False))
    sub=out/'final';sub.mkdir(exist_ok=True)
    base,cs,audit=one_fit(d,np.arange(len(d['edges'])),d['query'],knowledge,plan,sub)
    preds={'M0':base}
    for label,row in groups.items():preds[label]=base.copy() if row['key']=='M0' else base+row['beta']*cs[row['key']]
    np.savez_compressed(out/'PREDICTIONS.npz',edges=d['query'],ids=d['query_ids'],**preds,
        **{k+'__uncalibrated_delta':v for k,v in cs.items()})
    dump(sub/'audit.json',audit)
    dump(out/'FREEZE.json',dict(cell=cell,seconds=time.monotonic()-begin,code_sha256=sha(__file__),plan_sha256=sha(ROOT/'PLAN_PHASE2.json'),
        input_sha256=sha(ROOT/'inputs'/f'{cell}.json'),knowledge_sha256=sha(ROOT/'inputs'/f'{cell}_GenePT64.npz'),
        files={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()},outer_truth_read=False))
    print(json.dumps(dict(cell=cell,phase='frozen',selected=groups,seconds=time.monotonic()-begin)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cell',required=True);args=p.parse_args();run(args.cell)
