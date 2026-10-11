"""run_one takes latency_observable/completeness from the FIRST catalogue row of a virus,
summarise_programs from the LAST; validate_catalogue does not require them to agree per virus."""
import anndata as ad
import numpy as np
import scipy.sparse as sp

from viralscan.gene_programs import load_catalogue, summarise_programs, validate_catalogue
from viralscan.scripts.gene_programs import run_one

cat = [dict(r) for r in load_catalogue() if r["virus"] == "Epstein-Barr virus"]
cat[-1]["latency_observable_in_rna"] = "false"  # one row disagrees
cat[-1]["panel_completeness"] = "partial"
print("validate_catalogue problems:", validate_catalogue(cat))

latent = next(r for r in cat if r["programme"] == "latent")
x = sp.csr_matrix(np.array([[5.0]], dtype=np.float32))
a = ad.AnnData(x)
a.var_names = [latent["gene_id_bundled"]]
a.obs_names = ["c1"]
a.layers["counts_unique_viral"] = x.copy()
df = run_one(a, ["Epstein-Barr virus"], cat, 2, "bundled")
s = summarise_programs(df, cat, viruses=["Epstein-Barr virus"], marker_resolution=df.attrs["marker_resolution"])
print("per-cell state:", df.state.tolist(), "latency_not_observable:", df.latency_not_observable.tolist())
print("summary n_cells_latent:", s.n_cells_latent.tolist(), "latency_observable_in_rna:", s.latency_observable_in_rna.tolist())
