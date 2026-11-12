from common import *

def main():
    tasks=[];stats=[]
    for panel in ('full_raw','full_published'):
        inp=loadnpz(ROOT/'inputs'/f'{panel}_inputs.npz');truth=loadnpz(ROOT/'inputs'/f'{panel}_outcomes.npz')
        assert np.array_equal(inp['pairs'],truth['pairs'])
        genes,pairs,y=inp['genes'],inp['pairs'],truth['outcomes'];lut={g:i for i,g in enumerate(genes)};ed=np.array([[lut[a],lut[b]] for a,b in pairs])
        assert len(set(tuple(sorted(p)) for p in pairs))==len(pairs)
        for path in sorted((ROOT/'inputs/splits'/panel).glob('*/*/identities.json')):
            seed,fold=int(path.parent.parent.name),int(path.parent.name);ids=json.loads(path.read_text());tr,qu=np.array(ids['train']),np.array(ids['query'])
            assert len(set(tr)&set(qu))==0 and sorted(np.r_[tr,qu])==list(range(len(pairs)))
            dest=ROOT/'prepared'/panel/str(seed)/str(fold)
            savez(dest/'training.npz',genes=genes,singles=inp['singles'],edges=ed[tr],pairs=pairs[tr],y=y[tr],global_ids=tr,inner_fold=innerfold(pairs[tr]),query_edges=ed[qu],query_pairs=pairs[qu],query_ids=qu)
            savez(dest/'scoring_only.npz',query_ids=qu,y=y[qu])
            deg=np.bincount(ed[tr].ravel(),minlength=len(genes));st=dict(panel=panel,seed=seed,fold=fold,train=len(tr),query=len(qu),active_train=int(np.count_nonzero(deg)),degree_min_active=int(deg[deg>0].min()),degree_mean_active=float(deg[deg>0].mean()))
            if y.shape[1]==2:
                py=Project(ed[tr],len(genes)).p(y[tr]);st['training_partner_repeat_r']=float(np.corrcoef(py.T)[0,1])
            stats.append(st)
            for baseline in ('native','saturation'):tasks.append(dict(panel=panel,seed=seed,fold=fold,baseline=baseline))
    dump(ROOT/'TASKS.json',tasks);dump(ROOT/'TRAINING_PROPERTIES.json',stats)
    files=[ROOT/'PROTOCOL.md',ROOT/'TASKS.json',*sorted((ROOT/'inputs').rglob('*')),*sorted((ROOT/'code').rglob('*.py'))]
    dump(ROOT/'FREEZE.json',dict(utc=now(),status='fixed before new native Han fitting; previous Han outcomes known',files={str(p.relative_to(ROOT)):digest(p) for p in files if p.is_file()}))
    print('Prepared',len(tasks),'tasks; outer labels isolated; protocol and source hashes frozen.')
if __name__=='__main__':main()
