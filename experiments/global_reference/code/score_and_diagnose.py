"""Score locked matched-reference fits and diagnose original fitted functions."""
from refit_global import *
import pandas as pd

MITO=set('ATP5J2 MRPL3 MRPL13 MRPL18 MRPL27 MRPL32 MRPL33 MRPL43 MRPS14 MRPS27 MRPS30 MRPS35 NDUFA8 PET112 PTCD1 SSBP1'.split())
def corr(x,y):return float(np.corrcoef(x,y)[0,1]) if np.std(x)>0 and np.std(y)>0 else 0.
def ss(y,p,P):return float(np.sum((y-p)**2)),float(np.sum(P.p(y-p)**2))
def mode(G):
    vals,u=np.linalg.eigh(G);j=np.argmax(abs(vals));v=u[:,j]
    if v[np.argmax(abs(v))]<0:v=-v
    return float(vals[j]),v
def score():
    configs=json.loads((R/'CONFIGURATIONS.json').read_text());assert all((R/'fitted'/c['label']/'LOCK.json').exists() for c in configs)
    rows=[];diag=[];loadings=[];modules=[];atts=[];checks=[]
    for conf in configs:
        label=conf['label'];I=R/'inputs'/label;O=R/'fitted'/label
        lock=json.loads((O/'LOCK.json').read_text());assert sha(O/'PREDICTIONS.npz')==lock['predictions_sha256']
        tr=dict(np.load(I/'training.npz'));orig=dict(np.load(I/'original.npz'));new=dict(np.load(O/'PREDICTIONS.npz'));y=np.load(I/'outcomes.npz')['y'];ed=tr['edges'];qu=tr['query_edges'];s=tr['singles'];n=len(s);P=Project(qu,n);Pt=Project(ed,n);sel=json.loads((O/'SELECTION.json').read_text())
        predictions={'native_B':orig['M0'],'native_rank1':orig['rank1_selected'],'native_spectral':orig['spectral_selected'],'global_B':new['B'],'global_rank1':new['B']+new['rank1_D'],'global_spectral':new['B']+new['spectral_D'],'global_rank1_control':new['rank1_controlled'],'global_spectral_control':new['spectral_controlled']}
        for name,pred in predictions.items():
            for rep in [-1,0,1]:
                yy=y if rep==-1 else y[:,rep:rep+1];pp=pred if rep==-1 else pred[:,rep:rep+1];f,p=ss(yy,pp,P)
                ref=predictions['native_B' if name.startswith('native') else 'global_B'];rf,rp=ss(yy,ref if rep==-1 else ref[:,rep:rep+1],P)
                rows.append(dict(configuration=label,group=label.split('_fold')[0],method=name,repeat='pooled' if rep==-1 else rep+1,full_sse=f,partner_sse=p,reference_full_sse=rf,reference_partner_sse=rp,full_gain_pct=100*(1-f/rf),partner_gain_pct=100*(1-p/rp)))
        # Independent QR projection agrees with the Gram-pseudoinverse projector.
        q,_=np.linalg.qr(P.X,mode='reduced');qq=y-orig['rank1_selected'];qr=qq-q@(q.T@qq)
        checks.append(dict(configuration=label,projection_max_abs=float(np.max(abs(qr-P.p(qq)))),native_replay_max_abs=lock['native_baseline_replay_max_abs']))
        # Rebuild original out-of-fold baselines using the archived fitting rule.
        oo=np.zeros_like(tr['y'])
        for k in range(3):
            va=tr['folds']==k;st=fit_base(s,ed[~va],tr['y'][~va],False);oo[va]=predict(st,ed[va])
        D=orig['rank1_selected']-orig['M0'];d=P.p(D);r=P.p(y-orig['M0']);den=np.sum(r*r)
        singlemean=s.mean(axis=1)
        for template in ['mean_singles','repeat_singles']:
            st=np.repeat(singlemean[:,None],2,axis=1) if template=='mean_singles' else s
            ti=st[ed[:,0]]*st[ed[:,1]];tq=P.p(st[qu[:,0]]*st[qu[:,1]]);aa=np.zeros(2);bb=np.zeros(2)
            for k in range(3):
                va=tr['folds']==k;pk=Project(ed[va],n);tt=pk.p(ti[va]);rr=pk.p(tr['y'][va]-oo[va]);aa+=np.sum(tt*rr,axis=0);bb+=np.sum(tt*tt,axis=0)
            coef=np.divide(aa,bb,out=np.zeros(2),where=bb>1e-14);dp=tq*coef
            product_gain=100*(2*np.sum(r*dp)-np.sum(dp*dp))/den
            ct=np.divide(np.sum(tq*d,axis=0),np.sum(tq*tq,axis=0),out=np.zeros(2),where=np.sum(tq*tq,axis=0)>1e-14);dt=tq*ct;dr=d-dt
            gain=lambda z:100*(2*np.sum(r*z)-np.sum(z*z))/den
            g1,g2,g=gain(dt),gain(dr),gain(d);assert abs(g1+g2-g)<1e-9
            atts.append(dict(configuration=label,template=template,product_beta_repeat1=coef[0],product_beta_repeat2=coef[1],product_gain_pct=product_gain,original_rank1_gain_pct=g,global_direction_gain_pp=g1,orthogonal_gain_pp=g2,global_share_of_net_gain=g1/g if abs(g)>1e-12 else None))
        # Compare a fitted matrix mode to the mode of zero-filled training residuals.
        mods=json.loads((I/'M0_models.json').read_text());sc=np.std(s,axis=0);sc[sc==0]=1;btr=np.column_stack([fm.FittedModel.from_dict(m).predict(ed)*sc[j] for j,m in enumerate(mods)])
        pr=Pt.p(tr['y']-btr).mean(axis=1);T=np.zeros((n,n));T[ed[:,0],ed[:,1]]=pr;T[ed[:,1],ed[:,0]]=pr
        _,v=mode(T);origsel=json.loads((I/'selection.json').read_text());G=orig['rank1_G'] if origsel['rank1']['selected']!='G0' else np.zeros((n,n));lam,u=mode(G)
        if label.startswith('K562'):
            mask=np.array([str(g) in MITO for g in tr['genes']]);
            if u[mask].mean()<u[~mask].mean():u=-u
        else:mask=np.zeros(n,dtype=bool)
        diag.append(dict(configuration=label,n_genes=n,training_pairs=len(ed),query_pairs=len(qu),selected=origsel['rank1']['selected'],fitted_mode_single_abs_r=abs(corr(u,singlemean)) if np.linalg.norm(G)>0 else None,residual_mode_single_abs_r=abs(corr(v,singlemean)),fitted_residual_mode_abs_r=abs(corr(u,v)) if np.linalg.norm(G)>0 else None,lambda_value=lam,global_refit_rank1_selected=sel['rank1']['selected'],global_refit_spectral_selected=sel['spectral']['selected']))
        for j,g in enumerate(tr['genes']):loadings.append(dict(configuration=label,gene=str(g),single_mean=singlemean[j],mitochondrial=int(mask[j]),fitted_rank1_loading=u[j] if np.linalg.norm(G)>0 else None,training_residual_loading=v[j]))
        if label.startswith('K562_high'):
            contrib=np.sum(r*r-(r-d)**2,axis=1);newr=P.p(y-new['B']);newd=P.p(new['rank1_D']);contribnew=np.sum(newr*newr-(newr-newd)**2,axis=1)
            count=mask[qu].sum(axis=1)
            for mname,con in [('original',contrib),('global_refit',contribnew)]:
                for cl,cm in [('none',count==0),('one',count==1),('both',count==2),('any',count>=1)]:
                    modules.append(dict(configuration=label,model=mname,pair_class=cl,pairs=int(cm.sum()),total_pairs=len(qu),query_share=float(cm.mean()),gain_sse=float(con[cm].sum()),total_gain_sse=float(con.sum()),net_gain_share=float(con[cm].sum()/con.sum())))
    df=pd.DataFrame(rows);df.to_csv(R/'results/MATCHED_REFERENCE_SCORES.csv',index=False)
    pooled=df[df.repeat=='pooled'].groupby(['group','method'],as_index=False)[['full_sse','partner_sse','reference_full_sse','reference_partner_sse']].sum();pooled['full_gain_pct']=100*(1-pooled.full_sse/pooled.reference_full_sse);pooled['partner_gain_pct']=100*(1-pooled.partner_sse/pooled.reference_partner_sse);pooled.to_csv(R/'results/MATCHED_REFERENCE_POOLED.csv',index=False)
    for name,data in [('MODE_DIAGNOSTICS',diag),('GENE_LOADINGS',loadings),('MODULE_GAIN_ATTRIBUTION',modules),('GLOBAL_ATTRIBUTION',atts),('NUMERICAL_CHECKS',checks)]:pd.DataFrame(data).to_csv(R/'results'/f'{name}.csv',index=False)
    print(pooled[['group','method','full_gain_pct','partner_gain_pct']].to_string(index=False));print(pd.DataFrame(diag).to_string(index=False));print(pd.DataFrame(atts).to_string(index=False))
if __name__=='__main__':score()
