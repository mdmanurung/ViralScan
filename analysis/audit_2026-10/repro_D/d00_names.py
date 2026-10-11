"""Display names the identity table gives the programme-catalogue viruses."""
from viralscan import virus_identity as vi, virus_catalog as vc
from viralscan.gene_programs import load_catalogue
rows = vc.load_catalogue()
cat = vi._Catalogue(rows, {})
names = set(cat.names.values())
for v in sorted({r["virus"] for r in load_catalogue()}):
    print(repr(v), "in identity names:", v in names)
