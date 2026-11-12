"""Inspect archived Han fitted modes; no new fitting or pathway discovery."""
from refit_global import *
import pandas as pd
def run():
    src=R/'inputs/han_modes';rows=[];loads=[]
    for t in json.loads((src/'INDEX.json').read_text()):
        d=np.load(src/(t['panel']+'_singles.npz'));z=np.load(src/t['file']);n=int(z['n']);G=np.zeros((n,n))
        if 'upper' in z:
            ii=np.triu_indices(n);G[ii]=z['upper'];G[(ii[1],ii[0])]=z['upper']
        else:G=float(z['strength'])*np.outer(z['u'],z['u'])
        vals,U=np.linalg.eigh(G);order=np.argsort(-abs(vals));vals,U=vals[order],U[:,order];sg=d['singles'].mean(axis=1);rank=int(np.sum(abs(vals)>1e-8*max(1.,max(abs(vals)))))
        selected=t['selected']!='G0';m=vals[0]**2/np.sum(vals**2) if np.sum(vals**2)>0 else 0
        rows.append({k:v for k,v in t.items() if k!='file'}|dict(fitted_rank=rank,used_rank=rank if selected else 0,leading_single_abs_r=float(abs(np.corrcoef(U[:,0],sg)[0,1])) if selected else None,leading_squared_eigenvalue_share=m if selected else None))
        if selected:
            for k in range(min(3,rank)):
                u=U[:,k].copy()
                if u[np.argmax(abs(u))]<0:u=-u
                for j in np.argsort(-abs(u))[:10]:loads.append({k0:v for k0,v in t.items() if k0!='file'}|dict(mode=k+1,eigenvalue=vals[k],gene=str(d['genes'][j]),loading=u[j],single_mean=sg[j]))
    pd.DataFrame(rows).to_csv(R/'results/HAN_MODE_DIAGNOSTICS.csv',index=False);pd.DataFrame(loads).to_csv(R/'results/HAN_MODE_LOADINGS.csv',index=False)
    print(pd.DataFrame(rows).groupby(['panel','baseline','family']).agg(selected=('used_rank',lambda x:int((x>0).sum())),rank_min=('used_rank','min'),rank_max=('used_rank','max'),r_min=('leading_single_abs_r','min'),r_max=('leading_single_abs_r','max')).to_string())
if __name__=='__main__':run()
