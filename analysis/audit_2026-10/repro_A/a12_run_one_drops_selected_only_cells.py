"""run_one drops cells whose UNIQUE-layer breadth is 0 even when the selected (allocated) layer carries >=2 productive groups."""
import numpy as np, anndata as ad, scipy.sparse as sp
from viralscan.scripts import gene_programs as gp
from viralscan.gene_programs import load_catalogue, summarise_programs
cat=[r for r in load_catalogue() if r["virus"]=="Epstein-Barr virus"]
prod=[r for r in cat if r["programme"]=="productive" and r["non_overlapping"]=="true"][:3]
names=[r["gene_id_bundled"] for r in prod]
print("markers:",names)
# cell0: unique=0 everywhere, allocated (X) positive on all three; cell1: unique positive on all three
unique=sp.csr_matrix(np.array([[0,0,0],[1,1,1]],dtype=np.float32))
alloc =sp.csr_matrix(np.array([[2,2,2],[1,1,1]],dtype=np.float32))
a=ad.AnnData(X=alloc,obs={"x":["c0","c1"]},var={"v":names}); a.obs_names=["c0","c1"]; a.var_names=names
a.layers["counts_unique_viral"]=unique
df=gp.run_one(a,["Epstein-Barr virus"],cat,2,"bundled")
print(df[["barcode","state","selected_state","productive_breadth","selected_productive_breadth"]])
s=summarise_programs(df,cat,viruses=["Epstein-Barr virus"],marker_resolution=df.attrs["marker_resolution"])
print("summary n_cells_productive:",int(s.n_cells_productive[0]),"n_cells_productive_selected_layer:",int(s.n_cells_productive_selected_layer[0]),"(true selected-layer productive cells = 2)")
