"""Training-only raw-count diagnostics; original evaluation target unchanged."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import json,hashlib
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from fit_models import ROOT,load_input,Projection,dump

def corr(a,b):
    return float(np.corrcoef(a,b)[0,1]) if np.std(a)*np.std(b)>0 else None

formal=ROOT.parent/'shen_takeover/Shen2017_Formal_20261108'
source=ROOT/'references/canonical_counts.csv.gz'
raw=pd.read_csv(source if source.exists() else formal/'measurement/canonical_counts.csv')
summaries=[];temporal=[];saved=[]
for cell in ['HeLa','A549','293T']:
    d=load_input(ROOT/'inputs'/f'{cell}.json');allowed=set(d['ids']);genes=set(d['genes'])
    frame=raw[raw.cell_line==cell].copy()
    frame['pair']=['|'.join(sorted((a,b))) for a,b in zip(frame.gene_a,frame.gene_b)]
    frame=frame[frame['pair'].isin(allowed)|(frame.gene_a=='NTC')|(frame.gene_b=='NTC')].copy()
    # Exclude same-gene targeting, nontraining doubles before any quality summary.
    frame=frame[(frame.gene_a!=frame.gene_b)|(frame.gene_a=='NTC')]
    controls=frame[(frame.gene_a=='NTC')&(frame.gene_b=='NTC')&(frame['count']>=20)].copy()
    controls['log']=np.log2(controls['count']);ref=controls.groupby(['replicate','day'])['log'].mean().to_dict()
    slopes=[];trajectories=[];lowcounts=[]
    for (rep,cid),g in frame.groupby(['replicate','construct_id'],sort=False):
        pair=g['pair'].iloc[0];a=g.gene_a.iloc[0];b=g.gene_b.iloc[0]
        if a=='NTC' and b=='NTC':continue
        good=g[g['count']>=20];points=len(good)
        lowcounts.append(dict(pair=pair,rep=int(rep),points=points,low_fraction=float(np.mean(g['count']<20))))
        if points<3:continue
        t=good.day.to_numpy(float);z=np.array([np.log2(c)-ref[(rep,day)] for c,day in zip(good['count'],good.day)])
        tc=t-t.mean();slope=float(tc@(z-z.mean())/(tc@tc));err=z-z.mean()-slope*tc
        se=float(np.sqrt(np.sum(err*err)/(points-2)/(tc@tc)))
        ss=float(np.sum((z-z.mean())**2));r2=1-float(err@err)/ss if ss>0 else 0.
        slopes.append(dict(cell=cell,pair=pair,cid=cid,rep=int(rep),slope=slope,se=se,r2=r2,points=points))
        zz=dict(zip(t,z))
        if 3. in zz:
            for day in [14.,21.,28.]:
                if day in zz:trajectories.append(dict(pair=pair,cid=cid,rep=int(rep),day=day,slope=(zz[day]-zz[3.])/(day-3.)))
    sf=pd.DataFrame(slopes);pivot=sf.pivot(index=['pair','cid'],columns='rep',values='slope').dropna()
    counts=pivot.groupby(level=0).size();matched=pivot.loc[pivot.index.get_level_values(0).isin(counts[counts>=3].index)]
    means=matched.groupby(level=0).mean();pred=means.reindex(d['ids']).to_numpy()
    assert np.isfinite(pred).all();maxdiff=float(np.max(abs(pred-d['y'])))
    assert maxdiff<1e-12
    pairrows=matched[matched.index.get_level_values(0).isin(allowed)];center=pairrows-pairrows.groupby(level=0).transform('mean')
    selected_meta=sf.merge(pairrows.reset_index()[['pair','cid']],on=['pair','cid'])
    df=pd.DataFrame(lowcounts);df=df[df.pair.isin(allowed)]
    P=Projection(d['edges'],len(d['genes']));py=P.apply(d['y'])
    halves=[];half_ids=[]
    for pid,g in pairrows.groupby(level=0):
        if len(g)<6:continue
        keys=[(hashlib.sha256(('guidehalf|'+str(cid)).encode()).hexdigest(),k) for k,cid in enumerate(g.index.get_level_values(1))]
        ix=[k for _,k in sorted(keys)];a=g.iloc[ix[::2]].to_numpy().mean();b=g.iloc[ix[1::2]].to_numpy().mean()
        halves.append([a,b]);half_ids.append(pid)
    half_r=None;half_dim=None
    if len(halves)>len(d['genes']):
        lookup={g:i for i,g in enumerate(d['genes'])};edges=np.array([[lookup[x] for x in pid.split('|')] for pid in half_ids]);H=Projection(edges,len(lookup));hp=H.apply(np.array(halves));half_r=corr(hp[:,0],hp[:,1]);half_dim=int(len(edges)-np.linalg.matrix_rank(H.X))
    summary=dict(cell=cell,training_pairs=len(d['edges']),matched_constructs_median=float(pairrows.groupby(level=0).size().median()),
        matched_constructs_min=int(pairrows.groupby(level=0).size().min()),matched_constructs_max=int(pairrows.groupby(level=0).size().max()),
        canonical_training_max_abs=maxdiff,full_repeat_r=corr(pred[:,0],pred[:,1]),partner_repeat_r=corr(py[:,0],py[:,1]),
        single_repeat_r=corr(d['singles'][:,0],d['singles'][:,1]),construct_repeat_r=corr(pairrows[1],pairrows[2]),
        within_pair_construct_repeat_r=corr(center[1],center[2]),median_slope_fit_r2=float(selected_meta.r2.median()),
        median_slope_se=float(selected_meta.se.median()),fraction_construct_repeat_three_points=float(np.mean(selected_meta.points==3)),
        fraction_raw_training_counts_below20=float(df.low_fraction.mean()),guide_half_partner_r=half_r,
        guide_half_pairs=len(halves),guide_half_dimensions=half_dim)
    summaries.append(summary)
    td=pd.DataFrame(trajectories)
    for day in [14.,21.,28.]:
        v=td[td.day==day].pivot(index=['pair','cid'],columns='rep',values='slope').dropna()
        cnt=v.groupby(level=0).size();v=v[v.index.get_level_values(0).isin(cnt[cnt>=3].index)];mean=v.groupby(level=0).mean()
        ids=[p for p in d['ids'] if p in mean.index]
        lookup={g:i for i,g in enumerate(d['genes'])};e=np.array([[lookup[x] for x in p.split('|')] for p in ids]);y=mean.loc[ids].to_numpy();pp=Projection(e,len(lookup));py=pp.apply(y)
        temporal.append(dict(cell=cell,day=day,training_pairs=len(ids),full_repeat_r=corr(y[:,0],y[:,1]),partner_repeat_r=corr(py[:,0],py[:,1]),median_construct_count=float(cnt.reindex(ids).median())))
    saved.append(sf)
    print(json.dumps(summary),flush=True)
out=ROOT/'diagnostics';out.mkdir(exist_ok=True)
dump(out/'MEASUREMENT_SUMMARY.json',summaries);dump(out/'TEMPORAL_DIAGNOSTICS.json',temporal)
pd.concat(saved).to_csv(out/'TRAINING_CONSTRUCT_DIAGNOSTICS.csv.gz',index=False)
