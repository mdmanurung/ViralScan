#!/usr/bin/env python3
"""Write modelled single-CDS GTFs for the 3 PANEL-01 records whose GenBank entry has no CDS.

Coordinates come from the longest ATG ORF (Puumala L, CVA24) or from tblastn of the KFDV
polyprotein onto the genome (Alkhumra, whose record has frameshifts vs KFDV, so its span is
the aligned region, not a clean ORF).  Each model was checked by blastp/tblastn against an
annotated sibling (NC_005222.1 Hantaan L, NC_002058.3 poliovirus, NC_039218.1 KFDV).
Usage: python scripts/model_nocds_gtfs.py [outdir]   (default src/viralscan/data)
"""
import sys
from pathlib import Path

# acc, file stem, start, end, gene, product, note
MODELS = [
    ("JN860200.1", "Alkhumra_hemorrhagic_fever_virus_JN860200", 18, 10274, "POLY", "polyprotein",
     "modelled CDS: GenBank has none; span from tblastn of KFDV NC_039218.1 (96.7 % identity), record has frameshifts vs KFDV"),
    ("NC_005225.1", "Orthohantavirus_puumalaense_NC_005225", 37, 6507, "L", "RNA-dependent RNA polymerase",
     "modelled CDS: GenBank has none; longest ORF 2156 aa, 69 % identity to Hantaan L NC_005222.1"),
    ("D90457.1", "Coxsackievirus_A24_D90457", 751, 7395, "POLY", "polyprotein",
     "modelled CDS: GenBank has none; longest ORF 2214 aa, 81 % identity to poliovirus NC_002058.3"),
]

out = Path(sys.argv[1] if len(sys.argv) > 1 else "src/viralscan/data")
for acc, stem, s, e, gene, product, note in MODELS:
    gid = f"{acc}_{gene}"
    attrs = (f'gene_id "{gid}"; transcript_id "{gid}_t1"; gene_name "{gene}"; gene_biotype "protein_coding"; '
             f'product "{product}"; gene "{gene}"; note "{note}"; n_exons "1";')
    (out / f"{stem}.gtf").write_text(f"{acc}\tViralScan\texon\t{s}\t{e}\t.\t+\t0\t{attrs}\n")
    print("wrote", out / f"{stem}.gtf")
