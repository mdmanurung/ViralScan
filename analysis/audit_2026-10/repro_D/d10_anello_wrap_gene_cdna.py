"""The anellovirus GTF keeps the origin-wrapping gene's exons in transcript order
(2677:2824 then 1:80), but kb ref's own splitter (ngs_tools, the code `kb ref` runs) concatenates in
genomic order, so the indexed cDNA of KU243129.1's CDS is [1:80][2677:2824]."""
import random
import tempfile
from pathlib import Path

import ngs_tools as ngs

from viralscan.anellovirus import gtf_text_for, load_gene_table

random.seed(1)
acc = "KU243129.1"
seq = "".join(random.choice("ACGT") for _ in range(2824))
gene = next(r for r in load_gene_table() if r["accession"] == acc and r["origin_spanning"])
tmp = Path(tempfile.mkdtemp())
(tmp / "g.fa").write_text(f">{acc}\n{seq}\n")
(tmp / "g.gtf").write_text(gtf_text_for([acc], fasta_texts={acc: f">{acc}\n{seq}\n"}))
genes, txs = ngs.gtf.genes_and_transcripts_from_gtf(str(tmp / "g.gtf"), use_version=True)
out = ngs.fasta.split_genomic_fasta_to_cdna(str(tmp / "g.fa"), str(tmp / "cdna.fa"), genes, txs)
want = seq[2676:2824] + seq[0:80]  # transcript order for join(2677..2824,1..80)
got = None
with ngs.fasta.Fasta(out) as fa:
    for e in fa:
        if gene["transcript_id"] in e.header:
            got = e.sequence
print("transcript", gene["transcript_id"], "exons in GTF:", gene["exons"])
print("indexed cDNA == transcript-order CDS:", got == want)
print("indexed cDNA == genomic-order concat [1:80]+[2677:2824]:", got == seq[0:80] + seq[2676:2824])
