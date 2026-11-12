"""Rebuild every eligible Han pair, retaining previous QC; no model scoring."""
from common import *
import argparse,shutil
import pandas as pd

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--han-dir',type=Path,required=True);ap.add_argument('--knowledge',type=Path,required=True);args=ap.parse_args()
    if (ROOT/'data/PREPARED.json').exists():return
    h=args.han_dir;refs=pd.read_csv(h/'source/reference.tsv',sep='\t',header=None,names=['index','gene','sequence'])
    guide=refs.gene.to_numpy();n=len(guide);assert n==700
    counts={k:np.fromfile(h/'counts'/f'{run}.counts.u64',dtype='<u8').reshape(n,n) for k,run in dict(plasmid='SRR5258688',t0='SRR5258687',r1='SRR5258686',r2='SRR5258685').items()}
    valid=(counts['t0']>=50)&(counts['plasmid']>0);safe=guide=='0Safe';ctrl=valid&safe[:,None]&safe[None,:]
    ii,jj=np.where(valid);aa=guide[ii];bb=guide[jj];zs=[];norm=[]
    for key in ['r1','r2']:
        lfc=np.log2(np.maximum(counts[key][ii,jj],1)/counts['plasmid'][ii,jj])
        controls=lfc[safe[ii]&safe[jj]];center=float(np.median(controls));scale=float(np.std(controls,ddof=1))
        zs.append((lfc-center)/scale);norm.append(dict(repeat=key,center=center,scale=scale,n_controls=len(controls)))
    records=pd.DataFrame(dict(a=aa,b=bb,r1=zs[0],r2=zs[1]))
    sm=records[(records.a=='0Safe')^(records.b=='0Safe')].copy();sm['gene']=np.where(sm.a=='0Safe',sm.b,sm.a)
    sing=sm.groupby('gene')[['r1','r2']].median().sort_index()
    dd=records[(records.a!='0Safe')&(records.b!='0Safe')&(records.a!=records.b)].copy()
    dd['gene_a']=np.minimum(dd.a.to_numpy(),dd.b.to_numpy());dd['gene_b']=np.maximum(dd.a.to_numpy(),dd.b.to_numpy());dd['ab']=(dd.a==dd.gene_a).astype(int)
    group=dd.groupby(['gene_a','gene_b'])
    support=group.agg(constructs=('ab','size'),orientation_ab=('ab','sum'),r1=('r1','median'),r2=('r2','median')).reset_index()
    support['orientation_ba']=support.constructs-support.orientation_ab
    support['eligible']=(support.constructs>=3)&(support.orientation_ab>=1)&(support.orientation_ba>=1)&support.gene_a.isin(sing.index)&support.gene_b.isin(sing.index)
    support[['gene_a','gene_b','constructs','orientation_ab','orientation_ba','eligible']].to_csv(ROOT/'data/RAW_PAIR_SUPPORT.csv',index=False)
    good=support[support.eligible].sort_values(['gene_a','gene_b']);pairs=good[['gene_a','gene_b']].to_numpy(dtype=str);y=good[['r1','r2']].to_numpy();genes=sing.index.to_numpy(dtype=str);single=sing.to_numpy()
    raw=dict(genes=genes,singles=single,pairs=pairs);np.savez_compressed(ROOT/'data/full_raw_inputs.npz',**raw);np.savez_compressed(ROOT/'data/full_raw_outcomes.npz',pairs=pairs,outcomes=y)
    # Prior common panel replay is verified before any new predictions are scored.
    old=json.loads((h/'repeat_inputs/target_inputs.json').read_text());oi={tuple(p):i for i,p in enumerate(pairs)};ids=np.array([oi[tuple(p)] for p in old['pairs']])
    assert list(genes)==old['genes'];assert np.max(abs(single-np.array(old['singles'])))<1e-12
    old_y=load(h/'repeat_inputs/sealed_outcomes.npz')['outcomes'];assert np.max(abs(y[ids]-old_y))<1e-12
    np.savez_compressed(ROOT/'data/common_raw_inputs.npz',genes=genes,singles=single,pairs=pairs[ids]);np.savez_compressed(ROOT/'data/common_raw_outcomes.npz',pairs=pairs[ids],outcomes=y[ids])
    table=pd.read_excel(h/'source/41587_2017_BFnbt3834_MOESM267_ESM.xlsx',usecols=[0,1,2,3]);singles={};truth={}
    for label,a,b,v in table.itertuples(index=False,name=None):
        left,right=label.split('__')
        for gene,value in [(left,float(a)),(right,float(b))]:
            if gene in singles:assert singles[gene]==value
            singles[gene]=value
        if left!=right and np.isfinite(v):truth[tuple(sorted((left,right)))]=float(v)
    pg=sorted(g for g in singles if g not in ['0Safe','Safe','safe','NTC'] and np.isfinite(singles[g]));ps=np.array([[singles[g]] for g in pg]);pp=np.array(sorted(p for p in truth if set(p)<=set(pg)));py=np.array([[truth[tuple(p)]] for p in pp])
    np.savez_compressed(ROOT/'data/full_published_inputs.npz',genes=np.array(pg),singles=ps,pairs=pp);np.savez_compressed(ROOT/'data/full_published_outcomes.npz',pairs=pp,outcomes=py)
    common=np.array(old['pairs']);cy=np.array([[truth[tuple(p)]] for p in common]);np.savez_compressed(ROOT/'data/common_published_inputs.npz',genes=np.array(pg),singles=ps,pairs=common);np.savez_compressed(ROOT/'data/common_published_outcomes.npz',pairs=common,outcomes=cy)
    knowledge=load(args.knowledge);ki={g:i for i,g in enumerate(knowledge['genes'].astype(str))};mean=knowledge['embedding'].mean(0)
    kg=sorted(set(genes)|set(pg));missing=[g for g in kg if g not in ki]
    ke=np.array([knowledge['embedding'][ki[g]] if g in ki else mean for g in kg]);np.savez_compressed(ROOT/'data/target_knowledge.npz',genes=np.array(kg),embedding=ke,missing=np.array(missing))
    for name in ['common_raw','common_published','full_raw','full_published']:
        inp=load(ROOT/f'data/{name}_inputs.npz');out=load(ROOT/f'data/{name}_outcomes.npz');allids=np.arange(len(inp['pairs']))
        for seed in range(20261108,20261113):
            for fold,q in enumerate(split(inp['pairs'],allids,3,seed)):
                dest=ROOT/'data/splits'/name/str(seed)/str(fold);dest.mkdir(parents=True,exist_ok=True);tr=np.setdiff1d(allids,q)
                np.savez_compressed(dest/'calibration.npz',ids=tr,y=out['outcomes'][tr]);np.savez_compressed(dest/'query_outcomes.npz',ids=q,y=out['outcomes'][q]);dump(dest/'identities.json',dict(train=tr.tolist(),query=q.tolist()))
    receipt=dict(created_utc=now(),raw_pairs=len(pairs),raw_genes=len(set(pairs.ravel())),published_pairs=len(pp),published_genes=len(set(pp.ravel())),common_pairs=len(common),finite_raw_singles=len(genes),missing_knowledge=missing,normalization=norm,original_common_reconstruction_max_abs=0.,model_scores_computed=False,source_knowledge_sha256=sha(args.knowledge),source_workbook_sha256=sha(h/'source/41587_2017_BFnbt3834_MOESM267_ESM.xlsx'))
    dump(ROOT/'data/PREPARED.json',receipt);print(json.dumps(receipt,ensure_ascii=False),flush=True)
if __name__=='__main__':main()
