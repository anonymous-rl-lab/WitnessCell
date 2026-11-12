from common import *
import sys,tempfile
sys.path.insert(0,str(Path(__file__).parent/'vendor'))
import run_smoke as native
import baselines

def main():
    import fit_task
    stems=[f'{fit_task.key(f,t,n)}_f{k}' for k in range(3) for f in ('rank1','spectral') for t in fit_task.FRACTIONS for n in fit_task.NUS]
    assert len(set(stems))==48 and len(set(str(Path(s+'.npz')) for s in stems))==48
    rng=np.random.default_rng(3801);n=12;single=rng.normal(size=(n,2));ed=np.array([(i,j) for i in range(n) for j in range(i+1,n)]);pairs=np.array([[f'g{i}',f'g{j}'] for i,j in ed]);y=single[ed[:,0]]+single[ed[:,1]]+.1*rng.normal(size=(len(ed),2))
    with tempfile.TemporaryDirectory() as td:
        state=baselines.fit('native',single,ed[:50],y[:50],pairs[:50],20261108,Path(td)/'b.json');_,pred=native.baseline(single,ed[:50],y[:50],None);gap=float(np.max(abs(pred(ed[50:])-baselines.predict(state,single,ed[50:]))));assert gap==0
        state=baselines.fit('saturation',single,ed[:50],y[:50],pairs[:50],20261108,Path(td)/'s.json');z=baselines.predict(state,single,ed[50:]);assert z.shape==(16,2) and np.isfinite(z).all()
    p=Project(ed,n);a=rng.normal(size=(len(ed),2));ref=a-p.X@np.linalg.lstsq(p.X,a,rcond=None)[0];projection_gap=float(np.max(abs(ref-p.p(a))));assert projection_gap<1e-12
    dump(ROOT/'ENGINEERING_CHECKS.json',dict(native_baseline_adapter_max_abs=gap,projection_vs_lstsq_max_abs=projection_gap,saturation_smoke_finite=True,uses_Han_query_outcomes=False))
    print('Engineering checks passed: native baseline exact; projection agrees with independent least squares.')
if __name__=='__main__':main()
