"""Reconstruct locked predictions from saved fitted weights; no retraining."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,argparse,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from common import loadnpz,dump,Project,digest
import baselines

def matrix(path):
    v=loadnpz(path);n=int(v['n'])
    if 'upper' in v:
        g=np.zeros((n,n));ii=np.triu_indices(n);g[ii]=v['upper'];g[(ii[1],ii[0])]=v['upper'];return g
    return float(v['strength'])*np.outer(v['u'],v['u'])
def run(root,out):
    global ROOT
    ROOT=root;tasks=json.loads((root/'TASKS.json').read_text());rows=[];checks=[]
    for t in tasks:
        d=loadnpz(root/'inputs'/f"{t['panel']}_inputs.npz");yy=loadnpz(root/'inputs'/f"{t['panel']}_outcomes.npz")['outcomes'];ids=json.loads((root/'inputs/splits'/t['panel']/str(t['seed'])/str(t['fold'])/'identities.json').read_text());tr,qu=np.array(ids['train']),np.array(ids['query']);lut={g:i for i,g in enumerate(d['genes'])};ed=np.array([[lut[a],lut[b]] for a,b in d['pairs']]);n=len(d['genes']);o=root/'fitted'/t['panel']/t['baseline']/str(t['seed'])/str(t['fold']);sel=json.loads((o/'SELECTION.json').read_text());q=loadnpz(o/'QUERY.npz');inner=loadnpz(o/'INNER.npz');y=yy[qu];nr=y.shape[1];Pq=Project(ed[qu],n)
        b=baselines.predict(json.loads((o/'baseline.json').read_text()),d['singles'],ed[qu]);bgap=float(np.max(abs(b-q['B'])));assert bgap<1e-10
        predictions={'baseline':b,'additive':d['singles'][ed[qu,0]]+d['singles'][ed[qu,1]]}
        for family in ('rank1','spectral'):
            d_inner=inner[family+'_D'];A=H=0.;largest_inner_gap=0.
            for k in range(3):
                va=inner['folds']==k;G=matrix(o/f'{family}_inner{k}.npz');state=json.loads((o/f'inner{k}_baseline.json').read_text());base=baselines.predict(state,d['singles'],ed[tr[va]]);dd=np.repeat(G[ed[tr[va],0],ed[tr[va],1]][:,None],nr,axis=1)
                if sel[family]['selected']=='G0':dd*=0
                # Match the original arithmetic: (baseline + correction) - baseline.
                dd=(base+dd)-base;largest_inner_gap=max(largest_inner_gap,float(np.max(abs(dd-d_inner[va]))));assert np.max(abs(base-inner['B'][va]))<1e-10
                pp=Project(ed[tr[va]],n);r,z=pp.p(yy[tr[va]]-base),pp.p(dd);A+=np.sum(r*z);H+=np.sum(z*z)
            alpha=float(np.clip(A/H,0,1)) if H>1e-14 else 0.;assert abs(alpha-sel[family]['alpha'])<1e-10
            G=matrix(o/f'{family}_final.npz');D=np.repeat(G[ed[qu,0],ed[qu,1]][:,None],nr,axis=1)
            if sel[family]['selected']=='G0':D*=0
            dgap=float(np.max(abs(D-q[family+'_D'])));assert dgap<1e-10
            PD=Pq.p(D);gamma=sel[family]['gamma'];predictions.update({family+'_original':b+D,family+'_off':b+D-PD,family+'_controlled':b+D+(alpha-1)*PD,family+'_whole_same':b+alpha*D,family+'_whole_full':b+gamma*D})
            checks.append(t|dict(family=family,baseline_max_abs=bgap,inner_D_max_abs=largest_inner_gap,final_D_max_abs=dgap,alpha_max_abs=abs(alpha-sel[family]['alpha'])))
        for name,pred in predictions.items():
            e=y-pred;pe=Pq.p(e)
            for j in range(nr):rows.append(t|dict(method=name,repeat=j+1,full_sse=float(np.sum(e[:,j]**2)),partner_sse=float(np.sum(pe[:,j]**2))))
    df=pd.DataFrame(rows);ref=pd.read_csv(root/'results/FOLD_REPEAT_SCORES.csv');keys=['panel','baseline','seed','fold','method','repeat'];m=df.merge(ref,on=keys,suffixes=('_replay','_reference'),validate='one_to_one');assert len(m)==len(ref)
    gaps={k:float(np.max(abs(m[k+'_replay']-m[k+'_reference']))) for k in ('full_sse','partner_sse')};assert max(gaps.values())<1e-8
    out.mkdir(parents=True,exist_ok=True);df.to_csv(out/'REPLAY_SCORES.csv',index=False);pd.DataFrame(checks).to_csv(out/'WEIGHT_REPLAY_CHECKS.csv',index=False);dump(out/'REPLAY_CHECKS.json',dict(status='PASS',batches=len(tasks),rows=len(df),score_max_abs=gaps,max_alpha_abs=max(x['alpha_max_abs'] for x in checks),max_inner_D_abs=max(x['inner_D_max_abs'] for x in checks)))
    print(json.dumps(gaps))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--out',type=Path,default=ROOT/'replay_results');a=ap.parse_args();run(a.root.resolve(),a.out.resolve())
