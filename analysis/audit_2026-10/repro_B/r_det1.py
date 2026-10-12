import numpy as np, scipy.sparse as sp, pandas as pd, anndata as ad, warnings, traceback
warnings.filterwarnings("ignore")
from viralscan.scripts import detection as d
from viralscan.runconfig import RunConfig
import matplotlib; matplotlib.use("Agg")

def mk(X, genes, cells=None, dtype=None, sparse=False):
    X = np.asarray(X, dtype=dtype or float)
    if sparse: X = sp.csr_matrix(X)
    a = ad.AnnData(X)
    a.var_names = genes
    a.obs_names = cells or [f"c{i}" for i in range(a.n_obs)]
    return a
def T(name, f):
    try: print(name, "->", f())
    except Exception as e: print(name, "-> EXC", type(e).__name__, str(e)[:150])

# 1 detect_genes edge cases
T("detect empty", lambda: d.detect_genes(["h1","v1"], np.zeros((0,2)), ["v1"]))
T("detect all-zero", lambda: d.detect_genes(["h1","v1"], np.zeros((3,2)), ["v1"]))
T("detect accession absent", lambda: d.detect_genes(["h1","v1"], np.ones((3,2)), ["vX"]))
T("detect one UMI total thr1", lambda: d.detect_genes(["h1","v1"], np.array([[0,1],[0,0]]), ["v1"]))
T("detect NaN", lambda: d.detect_genes(["h1","v1"], np.array([[0,np.nan],[0,3]]), ["v1"]))
T("detect dup names", lambda: d.detect_genes(["v1","v1"], np.array([[5,0],[0,0]]), ["v1"]))
T("detect frac 0.4 thr1 group", lambda: d.detect_genes(["v1","v2"], np.array([[0.4,0.4],[0.4,0]]), ["v1","v2"], 1, {"V":["v1","v2"]}))
# int8 overflow
X=np.full((300,2),100,dtype=np.int8); 
T("detect int8 sum", lambda: d.detect_genes(["h","v"], X, ["v"],1))
X32=np.full((300,2),1.0,dtype=np.float32)
# 2 compute_stats
a=mk([[0,5,0],[0,0,0],[3,0,200]],["v1","v2","h"],sparse=True)
vs,pc=d.compute_stats(a,{}, {"V":["v1","v2"]}, ["v1","v2"], called_mask=[True,False,True], viral_count_matrix=a.X)
print(vs); print(pc)
# fractional tiny multimap value counted infected
a=mk([[1e-9,0,500],[0,0,500]],["v1","v2","h"])
vs,pc=d.compute_stats(a,{}, {"V":["v1"]}, ["v1"])
print("tiny frac infected:", vs["V"]["infected_cells"], vs["V"]["viral_molecules_total_est"])
# empty
a=mk(np.zeros((0,3)),["v1","v2","h"])
T("compute_stats 0 cells", lambda: d.compute_stats(a,{}, {"V":["v1"]}, ["v1"])[0])
a=mk([[0,0,0]],["v1","v2","h"])
T("compute_stats all zero 1 cell", lambda: d.compute_stats(a,{}, {"V":["v1"]}, ["v1"])[0])
# host_cells barcode suffix mismatch
a=mk([[5,0,500],[0,0,500],[2,0,600]],["v1","v2","h"],cells=["AAA","CCC","GGG"])
hc={"called":{"AAA-1","CCC-1","GGG-1"},"comparable":{"AAA-1","CCC-1","GGG-1"}}
vs,_=d.compute_stats(a,{}, {"V":["v1"]}, ["v1"], called_mask=np.ones(3,bool), host_cells=hc)
print("suffix mismatch:", {k:vs["V"][k] for k in ("infected_comparable","n_comparable_cells","pct_infected_comparable","infected_called","n_called_cells","pct_infected_called")})
# int32 overflow of total umi
X=sp.csr_matrix(np.array([[2_000_000_000,1,1],[2_000_000_000,0,0]],dtype=np.int32))
a=ad.AnnData(X); a.var_names=["v1","v2","h"]; a.obs_names=["a","b"]
print("sum int32 sparse:", d._sum_axis(a.X,1), d._sum_axis(a.X,1).dtype)
Xd=np.array([[2_000_000_000,1,1],[2_000_000_000,0,0]],dtype=np.int32)
print("sum int32 dense:", d._sum_axis(Xd,1), d._sum_axis(Xd,1).dtype)
Xf=np.full((1,20_000_000),1.0,dtype=np.float32)
print("sum f32 dense 2e7 ones:", d._sum_axis(Xf,1))
# 3 super_expressor empty / all-zero
import tempfile, os
cfgobj=RunConfig(); d.config=cfgobj
for nm,Xs in [("allzero",np.zeros((50,3))),("empty",np.zeros((0,3))),("one",np.array([[3,0,5.]]))]:
    a=mk(Xs,["v1","v2","h"])
    td=tempfile.mkdtemp(); os.makedirs(td+"/plots")
    T("super_expressor "+nm, lambda: d.super_expressor(a,"V",["v1"],td))
# 4 measure_positive_control expected=0
a=mk([[5,0,500]],["v1","v2","h"])
c=RunConfig(); 
class C: positive_control_gene="v1"; positive_control_expected_molecules=0; positive_control_scope=None
T("pos control expected 0", lambda: d.measure_positive_control(a,C()))
C.positive_control_expected_molecules=-5
T("pos control expected -5", lambda: d.measure_positive_control(a,C()))
# 5 build_sensitivity_outputs default index_viruses
T("below_gate_counts tuple", lambda: d.below_gate_counts(a,None,(),{}))
