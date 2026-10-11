"""Which viruses can reach 'reactivating' via the real packaged catalogue and run_one()?

One cell per virus carrying ONLY its immediate-early productive marker(s).
Then: the same VZV cell when the later-productive markers are not in the index
(PROG-14 reduced marker set) -- reactivating is an absence claim.
"""
import anndata as ad
import numpy as np
import scipy.sparse as sp

from viralscan.gene_programs import load_catalogue
from viralscan.scripts.gene_programs import run_one, summarise_programs

cat = load_catalogue()


def adata_for(virus, only_genes=None, hot=()):
    ids = [r["gene_id_bundled"] for r in cat if r["virus"] == virus]
    if only_genes is not None:
        ids = [i for i in ids if i in only_genes]
    X = np.zeros((1, len(ids)), dtype=np.float32)
    for h in hot:
        X[0, ids.index(h)] = 3
    a = ad.AnnData(sp.csr_matrix(X))
    a.var_names = ids
    a.obs_names = ["cell1"]
    a.layers["counts_unique_viral"] = sp.csr_matrix(X)
    return a


print("== one IE marker only, full index ==")
for v in sorted({r["virus"] for r in cat}):
    ie = [r["gene_id_bundled"] for r in cat
          if r["virus"] == v and r["programme"] == "productive" and r["kinetic_class"] == "immediate_early"]
    if not ie:
        print(f"{v:28s} no productive IE marker")
        continue
    a = adata_for(v, hot=[ie[0]])
    df = run_one(a, [v], cat, 2, "bundled")
    print(f"{v:28s} IE={ie[0]:26s} state={df.state.tolist()} latency_not_observable={df.latency_not_observable.tolist()}")

print("== VZV: ORF62 only; later productive markers absent from the index ==")
v = "Varicella-zoster virus"
keep = {"VARICELLA_HHV3_gp63"}
a = adata_for(v, only_genes=keep, hot=["VARICELLA_HHV3_gp63"])
df = run_one(a, [v], cat, 2, "bundled")
print(df[["virus_name", "state", "productive_breadth"]].to_string())
s = summarise_programs(df, cat, viruses=[v], marker_resolution=df.attrs["marker_resolution"])
print(s[["n_cells_reactivating", "n_markers_resolved", "n_markers_unresolved", "caveat"]].to_string())
