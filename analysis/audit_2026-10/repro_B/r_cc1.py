import numpy as np, anndata as ad, tempfile, os, warnings
warnings.filterwarnings("ignore")
from pathlib import Path
from viralscan.scripts import cellcalling as cc, detection as d
from viralscan.runconfig import RunConfig
td=tempfile.mkdtemp()
open(f"{td}/cells.txt","w").write("AAA-1\nCCC-1\n")
a=ad.AnnData(np.array([[5,0,500],[0,0,500],[2,0,600.]])); a.var_names=["v1","v2","h"]; a.obs_names=["AAA-1","CCC-1","GGG-1"]
class C: called_cells_file=f"{td}/cells.txt"; cell_calling="external"; output=td+"/"; host_index="x"
m=cc.call_cells(a,C()); print("mask",m)
hc=cc.external_host_cells(C(),200.0); print("host_cells",hc)
vs,_=d.compute_stats(a,{}, {"V":["v1"]}, ["v1"], called_mask=m, host_cells=hc)
print({k:vs["V"][k] for k in ("n_called_cells","infected_called","pct_infected_called","n_comparable_cells","infected_comparable","pct_infected_comparable")})
# knee edge
print(cc.knee_cells(np.array([0,0,0.]),10), cc.knee_cells(np.array([]),10))
# external empty / zero matches
for bad in ["", "ZZZ\n"]:
    open(f"{td}/c2.txt","w").write(bad)
    try: cc.external_cells(a.obs_names,f"{td}/c2.txt")
    except Exception as e: print(type(e).__name__, str(e)[:80])
# read_emptydrops_mask with NA
open(f"{td}/e.tsv","w").write("barcode\tis_cell\ttotal\nAAA-1\tNA\t3\nCCC-1\tTRUE\t900\n")
print(cc.read_emptydrops_mask(a.obs_names,f"{td}/e.tsv"))
print(cc.host_cells_from_tsv(f"{td}/e.tsv",200))
