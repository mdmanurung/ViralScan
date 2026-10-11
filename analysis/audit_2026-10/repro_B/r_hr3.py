import numpy as np, warnings
warnings.filterwarnings("ignore")
from viralscan.scripts import hostresponse as hr
rng=np.random.default_rng(3)
n=3000; G=2000
depth=np.exp(rng.normal(8,0.5,n)).astype(int)
# NB-ish counts, heterogeneous mean; independent of y except through depth
mu=np.exp(rng.normal(-5,1.5,G))   # per-UMI proportion
frac=mu/mu.sum()
lam=depth[:,None]*frac[None,:]
cnt=rng.poisson(lam).astype(float)
tot=cnt.sum(1,keepdims=True)
Xn=np.log1p(cnt/tot*1e4)
# y depends on depth only (confounded) -> logistic
p=1/(1+np.exp(-(np.log(depth)-8)*3-(-3)))
y=(rng.random(n)<p).astype(int); print("n_pos",y.sum())
d=hr._genome_wide_differential(Xn,[f"g{i}" for i in range(G)],y,tot[:,0])
print("depth-confounded null: frac p<0.05 =",(d.p_value<0.05).mean(),"n fdr<0.05 =",(d.fdr<0.05).sum())
# pure null y
y2=(rng.random(n)<0.05).astype(int)
d=hr._genome_wide_differential(Xn,[f"g{i}" for i in range(G)],y2,tot[:,0])
print("pure null: frac p<0.05 =",(d.p_value<0.05).mean(),"n fdr<0.05 =",(d.fdr<0.05).sum())
lowexp=mu<np.quantile(mu,0.3)
print("pure null low-expr genes frac p<.05:",(d.set_index('gene').loc[[f"g{i}" for i in np.where(lowexp)[0]],'p_value']<0.05).mean())
