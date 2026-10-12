import numpy as np, anndata as ad, scipy.sparse as sp, pandas as pd, os, sys, warnings
warnings.filterwarnings("ignore")
from viralscan.scripts import hostresponse as hr
S=sys.argv[1]
rng=np.random.default_rng(0)
n=200
bc=[f"AAAC{i:04d}" for i in range(n)]
depth=rng.integers(2000,6000,n)
host=rng.poisson(depth[:,None]/60.0, size=(n,60)).astype(float)
hn=[f"G{i}" for i in range(60)]
H=ad.AnnData(sp.csr_matrix(host)); H.var_names=hn; H.obs_names=bc
H.write_h5ad(f"{S}/host.h5ad")
# virus: ONE virus, two genes; each gene hits 12 distinct cells (>=1 UMI), union 24
V=np.zeros((n,2))
V[:12,0]=3; V[12:24,1]=3
Vd=ad.AnnData(sp.csr_matrix(V)); Vd.var_names=["vA_g1","vA_g2"]; Vd.obs_names=bc
Vd.write_h5ad(f"{S}/virus.h5ad")
from viralscan.virus_identity import GeneIdentity, VirusIdentityTable
# ORIGINAL (2 rows, by design: a legacy list of uncatalogued IDs says nothing about grouping):
#   open(f"{S}/acc.txt","w").write("vA_g1\nvA_g2\n")   and pass f"{S}/acc.txt" below
# The Run's identity table says both genes belong to ONE virus (post AUD-04 hostresponse pools them).
VirusIdentityTable(tuple(GeneIdentity(g,"X1","catalogued",True,"accession:X1","Virus A") for g in ("vA_g1","vA_g2"))).write_tsv(f"{S}/virus_identity.tsv")
prov=hr.run_hostresponse(f"{S}/virus.h5ad",f"{S}/host.h5ad",f"{S}/virus_identity.tsv",f"{S}/out",use_hvg=False,seeds=[0,1],n_stab_iter=5,detection_threshold=1)
print(pd.read_csv(f"{S}/out/hostresponse_status.tsv",sep="\t").to_string())
