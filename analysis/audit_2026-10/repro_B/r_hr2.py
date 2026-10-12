import numpy as np, warnings, pandas as pd
warnings.filterwarnings("ignore")
from viralscan.scripts import hostresponse as hr
def T(n,f):
    try: print(n,"->",f())
    except Exception as e: print(n,"-> EXC",type(e).__name__,str(e)[:120])
T("bh nan", lambda: hr._bh_fdr([0.01,0.02,np.nan]))
T("bh ties", lambda: hr._bh_fdr([0.01,0.01,0.5]))
T("bh empty", lambda: hr._bh_fdr([]))
rng=np.random.default_rng(1)
n=100; y=np.r_[np.ones(20),np.zeros(80)].astype(int); depth=rng.integers(1000,5000,n)
X=rng.normal(size=(n,5)); X[:,0]=3.0   # zero variance
X[:,1]=y*1.0+rng.normal(scale=.5,size=n) # real signal
d=hr._genome_wide_differential(X,list("abcde"),y,depth)
print(d)
# NaN in X
X2=X.copy(); X2[0,2]=np.nan
T("diff NaN", lambda: hr._genome_wide_differential(X2,list("abcde"),y,depth)[["gene","p_value","fdr"]].to_dict("records"))
# n small
T("diff n=3", lambda: hr._genome_wide_differential(X[:3],list("abcde"),y[:3],depth[:3]).shape)
# df: partial r with 2 covariates but df n-2: compare with statsmodels OLS
import statsmodels.api as sm
cov=np.column_stack([np.ones(n),np.log1p(depth)])
Xc=np.column_stack([y,cov])
res=sm.OLS(X[:,1],Xc).fit(); print("OLS p(y)",res.pvalues[0], "| diff p", d.set_index('gene').loc['b','p_value'])
# small n df effect
n=14; y2=np.r_[np.ones(7),np.zeros(7)].astype(int); dp=rng.integers(1000,5000,n); x=rng.normal(size=(n,1))+0.8*y2[:,None]
dd=hr._genome_wide_differential(x,["g"],y2,dp); 
Xc=np.column_stack([y2,np.ones(n),np.log1p(dp)]); r=sm.OLS(x[:,0],Xc).fit()
print("small n: diff p",dd.p_value[0],"OLS p",r.pvalues[0])
# label ties
c=np.array([3,3,1,0,0,0.]); dep=np.array([3000,3000,1000,5000,5000,5000.])
lab=hr._virus_presence_label(c,dep,3,"cpm"); print("cpm label n_pos expected 2 got",lab.sum(), lab)
# direction at r==0
d0=hr._genome_wide_differential(np.zeros((20,1)),["z"],np.r_[np.ones(10),np.zeros(10)].astype(int),np.arange(20)+100.)
print(d0)
# evalue
print(hr._e_value(2.0), hr._e_value(0.5), hr._e_value(1.0))
