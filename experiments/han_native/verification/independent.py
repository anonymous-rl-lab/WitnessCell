import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,numpy as np,pandas as pd
from scipy.linalg import lstsq
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'));from common import loadnpz,dump,digest
import baselines

def P(edges,n,a):
    X=np.zeros((len(edges),n));X[np.arange(len(edges)),edges[:,0]]=1;X[np.arange(len(edges)),edges[:,1]]=1
    return a-X@lstsq(X,a,lapack_driver='gelsy',check_finite=False)[0]
def main():
    freeze=json.loads((ROOT/'FREEZE.json').read_text())
    for name,h in freeze['files'].items():assert digest(ROOT/name)==h,('frozen file changed',name)
    rows=[];checks=[];tasks=json.loads((ROOT/'TASKS.json').read_text())
    for t in tasks:
        o=ROOT/'runs'/t['panel']/t['baseline']/str(t['seed'])/str(t['fold']);d=loadnpz(ROOT/'prepared'/t['panel']/str(t['seed'])/str(t['fold'])/'training.npz');q=loadnpz(o/'predictions.npz');y=loadnpz(ROOT/'prepared'/t['panel']/str(t['seed'])/str(t['fold'])/'scoring_only.npz')['y'];s=json.loads((o/'SELECTION.json').read_text());n=len(d['genes']);nr=y.shape[1]
        rebuilt=baselines.predict(json.loads((o/'baseline.json').read_text()),d['singles'],d['query_edges']);bgap=float(np.max(abs(rebuilt-q['B'])));assert bgap==0
        names=['baseline']+[f+'_'+p for f in ('rank1','spectral') for p in ('original','whole_same','off','controlled','whole_full')]
        errors=[y-(q['B'] if k=='baseline' else q[k]) for k in names];stack=np.concatenate(errors,axis=1);projected=P(d['query_edges'],n,stack)
        for k,name in enumerate(names):
            e,ep=errors[k],projected[:,k*nr:(k+1)*nr]
            for j in range(nr):rows.append(t|dict(method=name,repeat=j+1,full_sse=float(e[:,j]@e[:,j]),partner_sse=float(ep[:,j]@ep[:,j])))
        for family in ('rank1','spectral'):
            inner=loadnpz(o/f'{family}_selected_inner.npz');e=d['y']-inner['baseline'];di=inner['D'];a=b=0.
            for k in range(3):
                va=d['inner_fold']==k;v=P(d['edges'][va],n,np.c_[e[va],di[va]]);r,z=v[:,:nr],v[:,nr:];a+=np.sum(r*z);b+=np.sum(z*z)
            alpha=float(np.clip(a/b,0,1)) if b>1e-14 else 0.;agap=abs(alpha-s[family]['alpha']);assert agap<1e-10
            G=loadnpz(o/f'{family}_final.npz')['G'];dd=G[d['query_edges'][:,0],d['query_edges'][:,1]][:,None];dd=np.repeat(dd,nr,axis=1)
            if s[family]['selected']=='G0':dd*=0
            dgap=float(np.max(abs(dd-q[family+'_D'])));assert dgap==0
            projected_direction=P(d['query_edges'],n,dd);c=rebuilt+dd+(alpha-1)*projected_direction;cgap=float(np.max(abs(c-q[family+'_controlled'])));assert cgap<1e-10
            checks.append(t|dict(family=family,baseline_replay_max_abs=bgap,direction_replay_max_abs=dgap,alpha_qr_max_abs=agap,control_replay_max_abs=cgap))
    df=pd.DataFrame(rows);orig=pd.read_csv(ROOT/'results/FOLD_REPEAT_SCORES.csv');keys=['panel','baseline','seed','fold','method','repeat'];m=df.merge(orig,on=keys,suffixes=('_independent','_native'),validate='one_to_one')
    gaps={k:float(np.max(abs(m[k+'_independent']-m[k+'_native']))) for k in ('partner_sse','full_sse')};assert max(gaps.values())<1e-8
    df.to_csv(ROOT/'verification/INDEPENDENT_SCORES.csv',index=False);pd.DataFrame(checks).to_csv(ROOT/'verification/PREDICTION_AND_ALPHA_CHECKS.csv',index=False)
    dump(ROOT/'verification/INDEPENDENT_CHECKS.json',dict(status='PASS',tasks=len(tasks),score_rows=len(df),frozen_files=len(freeze['files']),sse_max_abs=gaps,max_alpha_abs=max(x['alpha_qr_max_abs'] for x in checks),max_control_prediction_abs=max(x['control_replay_max_abs'] for x in checks)))
    print(json.dumps(gaps))
if __name__=='__main__':main()
