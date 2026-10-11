"""Verify repeat_gene_merges.tsv k-mer claims on the NORMALISED transcripts (copy fully contained in survivor?)."""
import csv, re, glob
from pathlib import Path
from viralscan.gtf_normalise import normalise_gtf_file
ATTR=re.compile(r'(\w+)\s+"([^"]*)"'); cache=Path.home()/".cache/viralscan/ncbi"
comp=str.maketrans("ACGTacgt","TGCAtgca")
def rc(s): return s.translate(comp)[::-1]
def kmers(s,k=31):
    out=set()
    for i in range(len(s)-k+1):
        w=s[i:i+k]
        if set(w)<=set("ACGT"): out.add(min(w,rc(w)))
    return out
merges=list(csv.DictReader(open("src/viralscan/data/repeat_gene_merges.tsv"),delimiter="\t"))
files={}
for p in glob.glob("src/viralscan/data/*.gtf"):
    t=Path(p).read_text()
    for m in merges:
        if f'gene_id "{m["gene_id"]}"' in t: files[m["gene_id"]]=p
tx={}  # per file: gene -> list of (seq, strand, [exons])
def seqs(p):
    if p in tx: return tx[p]
    d={}; ex={}
    for l in normalise_gtf_file(p,cache):
        c=l.split("\t")
        if len(c)==9 and c[2]=="exon":
            a=dict(ATTR.findall(c[8])); ex.setdefault((a["gene_id"],a["transcript_id"]),[]).append((c[0],c[6],int(c[3]),int(c[4])))
    gen={}
    for (g,t),e in ex.items():
        sq,st=e[0][0],e[0][1]
        if sq not in gen: gen[sq]="".join(Path(cache/sq/f"{sq}.fasta").read_text().split("\n",1)[1].split()).upper()
        s="".join(gen[sq][a-1:b] for _,_,a,b in sorted(e,key=lambda x:x[2]))
        d.setdefault(g,set()).update(kmers(rc(s) if st=="-" else s))
    tx[p]=d; return d
for m in merges:
    p=files.get(m["gene_id"]); 
    if not p: print(m["gene_id"],"not in any bundled GTF"); continue
    d=seqs(p); a=d.get(m["gene_id"]); b=d.get(m["merged_into"])
    if a is None or b is None: print(m["gene_id"],"missing",a is None,b is None); continue
    print(f'{m["virus"][:22]:22} {m["gene_id"][-14:]:14} copy={len(a):5d} survivor={len(b):5d} copy_in_survivor={len(a&b):5d} copy_only={len(a-b)} surv_only={len(b-a)}')
