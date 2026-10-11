import numpy as np, math
from scipy import stats, optimize
from viralscan import sensitivity as s
rng=np.random.default_rng(0)
for d in (0.05,0.10,0.15,0.2,0.3):
    sub=rng.random((200000,90))<d
    # longest clean run >= 31
    ok=0
    # vectorised: compute run lengths
    run=np.zeros(200000,int); best=np.zeros(200000,int)
    for j in range(90):
        run=np.where(sub[:,j],0,run+1); best=np.maximum(best,run)
    print(d, "MC",(best>=31).mean(), "DP", round(s.fragment_capture_exact(d),4))
for T in (1,2,3,10,50):
    lam=optimize.brentq(lambda l: stats.poisson.sf(T-1,l)-0.95,1e-6,1e4)
    print("T",T,"exact lam",round(lam,3),"linear",round(s.LOD95_MOLECULES*T,3), "pdet(lam)",s.poisson_detection_probability(lam,T))
print(s.poisson_detection_probability(1e5,50), s.poisson_detection_probability(800,1))
print(s.lod95(0), s.lod95(1e6,capture=2.0), s.lod95(1e6,capture=0.5), s.lod95(1e6,capture=-1))
print(s.fragment_capture_exact(0.0), s.fragment_capture_exact(0.5), s.fragment_capture_exact(0.1,read_length=30))
