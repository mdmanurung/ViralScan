import numpy as np, anndata as ad, scipy.sparse as sp, tempfile, os, warnings, hashlib
warnings.filterwarnings("ignore")
from viralscan.scripts import umap as U
from viralscan.runconfig import RunConfig
print("px", U.px is not None)
def mk(n=400,G=300,viral=True,seed=0, hostdepth=True):
    rng=np.random.default_rng(seed)
    X=rng.poisson(rng.gamma(0.5,2,size=(1,G))*rng.uniform(1,5,size=(n,1))).astype(float)
    X[:,:2]=0
    if viral: X[:25,0]=2; X[10:20,1]=1
    a=ad.AnnData(sp.csr_matrix(X)); a.var_names=["v1","v2"]+[f"h{i}" for i in range(G-2)]; a.obs_names=[f"c{i}" for i in range(n)]
    return a
def run(a, found, tag):
    td=tempfile.mkdtemp(); os.makedirs(td+"/plots")
    U.config=RunConfig(output=td, min_counts=100, min_genes=20)
    try:
        U.umap(a,found)
        print(tag,"ok",sorted(os.listdir(td+"/plots")))
        return td
    except Exception as e:
        import traceback; print(tag,"EXC",type(e).__name__,str(e)[:200])
td=run(mk(),{"v1":50.,"v2":10.},"viral")
td2=run(mk(),{"v1":50.,"v2":10.},"viral again")
run(mk(viral=False),{},"noviral")
# reproducibility check: UMAP coords
def coords(a):
    td=tempfile.mkdtemp(); os.makedirs(td+"/plots"); U.config=RunConfig(output=td,min_counts=100,min_genes=20)
    U.umap(a,{"v1":50.,"v2":10.}); return open(td+"/plots/umap_binary.html").read()
import re
h1=coords(mk()); h2=coords(mk())
print("html equal across runs:", h1==h2)
# all cells filtered
run(mk(n=50),{"v1":5.},"filter-all?")
U.config=RunConfig(output=tempfile.mkdtemp(), min_counts=10**9, min_genes=20)
os.makedirs(U.config.output+"/plots",exist_ok=True)
print("all-filtered return:", U.umap(mk(),{"v1":5.}))
# neighbor enrichment edge
print(U.viral_neighbor_enrichment(np.random.rand(30,2),np.ones(30,int),10,50))
print(U.viral_neighbor_enrichment(np.random.rand(30,2),np.zeros(30,int),10,50))
print(U.calculate_k_neighbors(5), U.calculate_k_neighbors(100000))
try: print(U.viral_neighbor_enrichment(np.random.rand(8,2),np.r_[1,1,0,0,0,0,0,0],10,50))
except Exception as e: print("k>n EXC",type(e).__name__,str(e)[:100])
print("---- coords repro")
import re
n1=re.sub(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}","U",h1); n2=re.sub(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}","U",h2)
print("html equal modulo uuid:", n1==n2)
