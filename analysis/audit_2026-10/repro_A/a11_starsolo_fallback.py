"""Latent defect: with a row flagged available_in_starsolo=false, panel_form='bundled' still falls back to the bare STARsolo ID
(shared by >=7 genomes) -> would bind another genome's gene. Shipped TSV has zero such rows, so only reachable via catalogue=."""
from viralscan.gene_programs import resolve_markers, load_catalogue
cat=[dict(r) for r in load_catalogue() if r["virus"]=="Epstein-Barr virus"][:1]
cat[0].update(gene_id_bundled="EPSTEIN_unassigned_gene_3",gene_id_merged="EPSTEIN_unassigned_gene_3",gene_id_starsolo="unassigned_gene_3",available_in_starsolo="false")
for form in ("bundled","starsolo"):
    res,unres=resolve_markers("Epstein-Barr virus",{"unassigned_gene_3"},panel_form=form,catalogue=cat)
    print(form,"resolved",[m.var_name for m in res],"unresolved",len(unres))
