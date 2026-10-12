"""EBV LMP-2A/2B: native exon rows with a wrapped (part 1/2) transcript pass through gtf_normalise untouched.
Show what ngs_tools (kb ref path) builds for them."""
import re, sys, tempfile
from pathlib import Path
import ngs_tools as ngs
from viralscan.gtf_normalise import normalise_gtf_file
cache = Path.home()/".cache/viralscan/ncbi"
gtf = "src/viralscan/data/Epstein_Barr_virus_NC_007605.gtf"
new = normalise_gtf_file(gtf, cache)
sel = [l for l in new if 'HHV4_LMP-2A"' in l]
print("rows for LMP-2A after normalise:", len(sel), "rule tags:", {m for l in sel for m in re.findall(r'viralscan_norm "([^"]+)"', l)})
w = Path(tempfile.mkdtemp(dir="/tmp/claude-149488/-exports-archive-hg-funcgenom-research-mdmanurung-ViralScan/a0c3f363-e67e-4594-9357-4e376c79d18a/scratchpad/auditA"))
(w/"x.gtf").write_text("\n".join(sel)+"\n")
fa = (cache/"NC_007605.1"/"NC_007605.1.fasta").read_text()
(w/"g.fa").write_text(fa)
genes, tx = ngs.gtf.genes_and_transcripts_from_gtf(str(w/"x.gtf"), use_version=True)
ngs.fasta.split_genomic_fasta_to_cdna(str(w/"g.fa"), str(w/"c.fa"), genes, tx)
seqs = {e.name:e.sequence for e in ngs.fasta.Fasta(str(w/"c.fa"),"r")}
g = "".join(fa.split("\n",1)[1].split())
for k,s in seqs.items():
    print(k, len(s))
    ex = [(166040,166458),(58,272),(360,458),(540,788),(871,951),(1026,1196),(1280,1495),(1574,1682),(5408,5856)]
    in_order = "".join(g[a-1:b] for a,b in ex)
    genomic = "".join(g[a-1:b] for a,b in sorted(ex))
    print(" matches biological (exon1 first):", s==in_order, " matches genomic order (exon1 last):", s==genomic, "expected len", len(in_order))
