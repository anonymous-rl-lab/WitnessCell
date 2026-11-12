"""Finite numerical checks accompany proofs; they do not establish cell mechanisms."""
from refit_global import *
def run():
    rng=np.random.default_rng(39024);N=100000;lo=rng.uniform(.01,3,N);hi=rng.uniform(.01,3,N);d=rng.uniform(-8,8,N);c=rng.uniform(-8,8,N)
    f=lambda x:np.clip(x,-lo,hi);eps=f(d+c)-f(d)-f(c);excess=np.max(abs(eps)-np.minimum(abs(d),abs(c)));assert excess<1e-12
    dp=abs(d);cp=abs(c);ep=f(dp+cp)-f(dp)-f(cp);exact=-np.maximum(0,np.minimum.reduce([dp,cp,hi,dp+cp-hi]));err=float(np.max(abs(ep-exact)));assert err<1e-12
    ncases=1000;worst=0.
    for _ in range(ncases):
        n=8;W=rng.uniform(0,1,(n,n));W*=.6/W.sum(axis=1).max();h=rng.uniform(.1,3,n);p=rng.uniform(0,1,n);q=rng.uniform(0,1,n);w=rng.uniform(.1,2,n)
        def solve(v):
            a=np.zeros(n)
            for i in range(200):
                b=np.minimum(h,v+W@a)
                if np.max(abs(a-b))<1e-13:return b
                a=b
            raise AssertionError('fixed point')
        a,b,ab=[float(w@solve(v)) for v in [p,q,p+q]];worst=max(worst,max(a,b)-ab,ab-a-b)
    assert worst<1e-10
    a=np.array([.2,.5,.8]);b=np.array([.3,.6,.9]);G=np.minimum(a[:,None]+b[None,:],1)-a[:,None]-b[None,:];H=np.eye(3)-np.ones((3,3))/3;assert np.linalg.matrix_rank(G)==3 and np.linalg.matrix_rank(H@G@H)==2
    # Endpoint-equivalent block descriptions give the same projected contrast.
    n=10;edges=np.array([(i,j) for i in range(n) for j in range(i+1,n)]);m=np.arange(n)<4;aa,bb,cc=.14,.08,.09;inside=m[edges].all(axis=1);cross=m[edges].sum(axis=1)==1;outside=~m[edges].any(axis=1);E=aa*inside+bb*cross+cc*outside;P=Project(edges,n);gap=float(np.max(abs(P.p(E[:,None])-(aa-2*bb+cc)*P.p(inside[:,None].astype(float)))));assert gap<1e-12
    dump(R/'results/THEORY_CHECKS.json',dict(status='PASS',clipping_samples=N,latent_bound_max_excess=float(excess),positive_clipping_formula_max_error=err,cooperative_networks=ncases,network_envelope_max_excess=worst,one_intermediate_raw_rank=3,one_intermediate_endpoint_removed_rank=2,endpoint_block_identity_max_error=gap,scope='Checks of specified mathematical constructions, not biological identification.'))
if __name__=='__main__':run()
