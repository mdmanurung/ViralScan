"""For each curated marker: other genes sharing >=T of the marker's 31-mers (normalised transcripts), and whether the
build_overlap_groups partition (after PROG-13 merges) puts them in the same group. Exposes sequence sharing invisible to exon-coordinate overlap."""
import sys, re, collections, itertools
from pathlib import Path
sys.path.insert(0,"extras")
import build_gene_programs as G
from viralscan.gtf_normalise import normalise_gtf_file
from viralscan.repeat_merges import load_repeat_merges
ATTR=re.compile(r'(\w+)\s+"([^"]*)"'); cache=Path.home()/".cache/viralscan/ncbi"
comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
merges=load_repeat_merges()
def kmers(s,k=31):
    o=set()
    for i in range(len(s)-k+1):
        w=s[i:i+k]
        if set(w)<=set("ACGT"): o.add(min(w,rc(w)))
    return o
THR=float(sys.argv[1]) if len(sys.argv)>1 else 0.5
print("virus\tmarker\tprogramme\tother_gene\tfrac_of_marker_kmers_shared\tother_is_marker\tsame_overlap_group\tmerged_by_PROG13")
for virus,gtf in sorted(G.VIRUS_GTF.items()):
    p=f"src/viralscan/data/{gtf}"; rec=G.parse_gtf(p)
    ex=collections.defaultdict(list)
    for l in normalise_gtf_file(p):
        c=l.split("\t")
        if len(c)==9 and c[2]=="exon":
            a=dict(ATTR.findall(c[8])); ex[(a["gene_id"],a["transcript_id"])].append((c[0],c[6],int(c[3]),int(c[4])))
    genome={}; K=collections.defaultdict(set)
    for (g,t),e in ex.items():
        sq,st=e[0][0],e[0][1]
        if sq not in genome: genome[sq]="".join((cache/sq/f"{sq}.fasta").read_text().split("\n",1)[1].split()).upper()
        K[g]|=kmers("".join(genome[sq][a-1:b] for _,_,a,b in sorted(e,key=lambda x:x[2])))
    # apply merges to K and rec like main()
    for c_,s_ in merges.items():
        if c_ in rec and s_ in rec:
            g=rec.pop(c_); rec[s_]["blocks"]|=g["blocks"]; K[s_]|=K.pop(c_,set())
    grp=G.build_overlap_groups(rec)
    marks={}
    for e in G.CURATED_MARKERS:
        if e["virus"]==virus:
            for g in G._resolve(rec,e["refseq_gene"]): marks[g]=e["programme"]
    for m,prog in sorted(marks.items()):
        for o,ko in K.items():
            if o==m or not K[m]: continue
            f=len(K[m]&ko)/len(K[m])
            if f>=THR:
                print(virus,m.split("_",2)[-1],prog,o.split("_",2)[-1],f"{f:.2f}",marks.get(o,"-"),grp[m]==grp[o],"",sep="\t")
