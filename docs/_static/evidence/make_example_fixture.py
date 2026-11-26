"""Synthetic multi-gene viral fixture behind docs/_static/evidence/coverage_map_synthetic.png.

Usage: python make_example_fixture.py OUTDIR
Seeded random sequence (no real virus): a 4,000 bp genome with three genes, a 1,500 bp host
transcript, 40 cells with real 10x v3 on-list barcodes, 114 read pairs.
"""
import random, sys
from pathlib import Path
out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
rng = random.Random(7)
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
genome = rand(4000)
host = rand(1500)
genes = [("VG_alpha", 150, 1250, "+"), ("VG_beta", 1500, 2400, "-"), ("VG_gamma", 2700, 3800, "+")]
(out/"reference.fasta").write_text(f">VIRAL_GENOME\n{genome}\n>host_tx\n{host}\n")
(out/"viral.fasta").write_text(f">VIRAL_GENOME\n{genome}\n")
(out/"host.fasta").write_text(f">host_tx\n{host}\n")
gtf = []
for g, s, e, st in genes:
    a = f'gene_id "{g}"; transcript_id "{g}_tx"; gene_name "{g}";'
    gtf.append(f"VIRAL_GENOME\tViralScan\ttranscript\t{s}\t{e}\t.\t{st}\t.\t{a}")
    gtf.append(f"VIRAL_GENOME\tViralScan\texon\t{s}\t{e}\t.\t{st}\t.\t{a}")
a = 'gene_id "HOST_GENE"; transcript_id "host_tx"; gene_name "HOST_GENE";'
gtf.append(f"host_tx\tViralScan\ttranscript\t1\t1500\t.\t+\t.\t{a}")
gtf.append(f"host_tx\tViralScan\texon\t1\t1500\t.\t+\t.\t{a}")
(out/"combined.gtf").write_text("\n".join(gtf)+"\n")
(out/"host.gtf").write_text("\n".join(gtf[-2:])+"\n")
# molecules: weights favour alpha and gamma; 3'-biased starts
cells = ['AAACCCAAGAAACACT', 'AAACCCAAGAAACCAT', 'AAACCCAAGAAACCCA', 'AAACCCAAGAAACCCG', 'AAACCCAAGAAACCTG', 'AAACCCAAGAAACGAA', 'AAACCCAAGAAACGTC', 'AAACCCAAGAAACTAC', 'AAACCCAAGAAACTCA', 'AAACCCAAGAAACTGC', 'AAACCCAAGAAACTGT', 'AAACCCAAGAAAGAAC', 'AAACCCAAGAAAGACA', 'AAACCCAAGAAAGCCT', 'AAACCCAAGAAAGCGA', 'AAACCCAAGAAAGGAT', 'AAACCCAAGAAAGGTA', 'AAACCCAAGAAAGTCT', 'AAACCCAAGAAATAGG', 'AAACCCAAGAAATCCA', 'AAACCCAAGAAATGAG', 'AAACCCAAGAAATGCC', 'AAACCCAAGAAATGGG', 'AAACCCAAGAAATTCG', 'AAACCCAAGAAATTGC', 'AAACCCAAGAACAAGG', 'AAACCCAAGAACACCA', 'AAACCCAAGAACACGT', 'AAACCCAAGAACAGGA', 'AAACCCAAGAACCCGA', 'AAACCCAAGAACCCGG', 'AAACCCAAGAACCGAT', 'AAACCCAAGAACCGCA', 'AAACCCAAGAACCTCC', 'AAACCCAAGAACCTTT', 'AAACCCAAGAACGAGG', 'AAACCCAAGAACGCGT', 'AAACCCAAGAACGCTA', 'AAACCCAAGAACGTGC', 'AAACCCAAGAACGTGT', 'AAACCCAAGAACTACG', 'AAACCCAAGAACTCCT', 'AAACCCAAGAACTCTT', 'AAACCCAAGAACTGAC', 'AAACCCAAGAACTGAT', 'AAACCCAAGAACTGGT', 'AAACCCAAGAACTTCC', 'AAACCCAAGAACTTTC', 'AAACCCAAGAAGAAAC', 'AAACCCAAGAAGAACG', 'AAACCCAAGAAGAATG', 'AAACCCAAGAAGAGCA', 'AAACCCAAGAAGAGGC', 'AAACCCAAGAAGAGTA', 'AAACCCAAGAAGATAG', 'AAACCCAAGAAGATCT', 'AAACCCAAGAAGCCAC', 'AAACCCAAGAAGCCAT', 'AAACCCAAGAAGCCCG', 'AAACCCAAGAAGCCGC'][:40]
r1, r2 = [], []
n = 0
for ci, cb in enumerate(cells):
    for _ in range(rng.choice([1, 1, 2, 3])):
        g, s, e, st = rng.choices(genes, weights=[5, 2, 4])[0]
        tx = genome[s-1:e]
        if st == "-": tx = rc(tx)
        # 3' bias: read starts within the last 600 nt of the transcript
        lo = max(0, len(tx) - 600 - 90); start = rng.randint(lo, len(tx) - 90)
        seq = tx[start:start+90]
        umi = rand(12)
        n += 1
        r1.append(f"@m{n}/1\n{cb}{umi}\n+\n{'I'*28}\n")
        r2.append(f"@m{n}/2\n{seq}\n+\n{'I'*90}\n")
    # a few host reads
    hs = rng.randint(0, 1400)
    n += 1
    r1.append(f"@m{n}/1\n{cb}{rand(12)}\n+\n{'I'*28}\n"); r2.append(f"@m{n}/2\n{host[hs:hs+90]}\n+\n{'I'*90}\n")
(out/"R1.fastq").write_text("".join(r1)); (out/"R2.fastq").write_text("".join(r2))
(out/"called_cells.txt").write_text("\n".join(cells)+"\n")
print(n, "pairs")
