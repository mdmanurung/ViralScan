"""`viralscan --reference -fasta host_genome.fa -gtf host.gtf`: every GTF gene_id is passed as 'viral'
(menu._build_kb_ref -> prepare_reference_inputs(viral_gene_ids=gtf_gene_ids(gtf))), so host sequence is dust-masked
and later labelled viral by write_build_manifest_from_t2g. Needs dustmasker on PATH. Run with PYTHONPATH=src."""
import random, tempfile
from pathlib import Path
from viralscan.scripts.build_reference import prepare_reference_inputs
from viralscan.virus_identity import gtf_gene_ids
random.seed(1)
seq = "".join(random.choice("ACGT") for _ in range(600)) + "AC" * 200 + "".join(random.choice("ACGT") for _ in range(600))
with tempfile.TemporaryDirectory() as d:
    d = Path(d)
    (d / "host.fa").write_text(f">chr1 host genome\n{seq}\n")
    (d / "host.gtf").write_text('chr1\tens\texon\t1\t300\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001";\n')
    prepared, _ = prepare_reference_inputs(d / "host.fa", d / "host.gtf", d / "idx",
                                           viral_gene_ids=gtf_gene_ids(d / "host.gtf"), mask=True)
    out = "".join(l for l in prepared.read_text().splitlines() if not l.startswith(">"))
    print("host chr1 'viral'-masked bases:", sum(c in "Nn" for c in out), "of", len(out))
