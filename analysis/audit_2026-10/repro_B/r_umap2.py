import numpy as np, anndata as ad, scipy.sparse as sp, tempfile, warnings
warnings.filterwarnings("ignore")
from viralscan.scripts import umap as U
from viralscan.runconfig import RunConfig
rng=np.random.default_rng(0)
X=rng.poisson(2,size=(200,100)).astype(float); a=ad.AnnData(sp.csr_matrix(X)); a.var_names=[f"g{i}" for i in range(100)]
td=tempfile.mkdtemp()  # no plots/ dir (visual off)
U.config=RunConfig(output=td, min_counts=10, min_genes=5)
try: U.umap(a,{}); print("ok")
except Exception as e: print("EXC",type(e).__name__,str(e)[:150])
