"""Gene pairs (same GTF) whose normalised-transcript 31-mer sets are >=90% contained in each other but are NOT in repeat_gene_merges.tsv.
Flags curated markers. Pure read-only."""
import sys, re, collections, itertools
from pathlib import Path
sys.path.insert(0,"extras")
import build_gene_programs as G
from viralscan.gtf_normalise import normalise_gtf_file
from viralscan.repeat_merges import load_repeat_merges
ATTR=re.compile(r'(\w+)\s+"([^"]*)"'); cache=Path.home()/".cache/viralscan/ncbi"
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
merges=load_repeat_merges(); merged=set(merges)|set(merges.values())
def kmers(s,k=31):
    o=set()
    for i in range(len(s)-k+1):
        w=s[i:i+k]
        if set(w)<=set("ACGT"): o.add(min(w,rc(w)))
    return o
marker_ids=collections.defaultdict(set)
for virus,gtf in sorted(G.VIRUS_GTF.items()):
    p=f"src/viralscan/data/{gtf}"; rec=G.parse_gtf(p)
    marks={}
    for e in G.CURATED_MARKERS:
        if e["virus"]==virus:
            for g in G._resolve(rec,e["refseq_gene"]): marks[g]=e["programme"]
    ex=collections.defaultdict(list)
    for l in normalise_gtf_file(p):
        c=l.split("\t")
        if len(c)==9 and c[2]=="exon":
            a=dict(ATTR.findall(c[8])); ex[(a["gene_id"],a["transcript_id"])].append((c[0],c[6],int(c[3]),int(c[4])))
    genome={}; K=collections.defaultdict(set)
    for (g,t),e in ex.items():
        sq,st=e[0][0],e[0][1]
        if sq not in genome: genome[sq]="".join((cache/sq/f"{sq}.fasta").read_text().split("\n",1)[1].split()).upper()
        s="".join(genome[sq][a-1:b] for _,_,a,b in sorted(e,key=lambda x:x[2]))
        K[g]|=kmers(s)
    idx=collections.defaultdict(list)
    for g,ks in K.items():
        for k in ks: idx[k].append(g)
    sh=collections.Counter()
    for k,gs in idx.items():
        if len(gs)>1:
            for a,b in itertools.combinations(sorted(gs),2): sh[(a,b)]+=1
    for (a,b),n in sh.items():
        ca=n/len(K[a]); cb=n/len(K[b])
        if max(ca,cb)>=0.9:
            tag="MERGED" if (a in merged and b in merged) else ""
            mk=",".join(f"{x.split('_',2)[-1]}:{marks[x]}" for x in (a,b) if x in marks)
            print(f"{virus[:18]:18} {a.split('_',2)[-1][:18]:18} {b.split('_',2)[-1][:18]:18} shared={n} frac_a={ca:.2f} frac_b={cb:.2f} {tag} {('MARKER '+mk) if mk else ''}")
