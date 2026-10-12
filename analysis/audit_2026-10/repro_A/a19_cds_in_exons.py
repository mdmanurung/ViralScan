"""After normalisation, is every CDS/stop_codon row covered by the exons of the transcript it was retagged to?
(Collapsed HHV-6B join rows legitimately span introns -> reported separately.)"""
import glob, re, collections
from pathlib import Path
from viralscan.gtf_normalise import normalise_gtf_file
ATTR=re.compile(r'(\w+)\s+"([^"]*)"'); cache=Path.home()/".cache/viralscan/ncbi"
bad=collections.Counter(); ex_bad=[]
for p in sorted(glob.glob("src/viralscan/data/*.gtf")):
    ex=collections.defaultdict(list); cds=[]
    for l in normalise_gtf_file(p,cache):
        c=l.split("\t")
        if len(c)!=9: continue
        a=dict(ATTR.findall(c[8]))
        if c[2]=="exon": ex[a["transcript_id"]].append((int(c[3]),int(c[4])))
        elif c[2] in("CDS","stop_codon","start_codon") and "viralscan" not in c[1]:
            cds.append((a.get("transcript_id"),c[2],int(c[3]),int(c[4]),a.get("gene_id"),a.get("part")))
    for t,f,s,e,g,part in cds:
        if t not in ex: bad["no_exon_for_tx"]+=1; continue
        cov=any(es<=s and e<=ee for es,ee in ex[t])
        if not cov:
            # covered by union?
            u=sum(max(0,min(e,ee)-max(s,es)+1) for es,ee in ex[t])
            kind="partial_cover" if u>0 else "no_cover"
            bad[(kind,"part" if part else "")]+=1
            ex_bad.append((Path(p).name[:30],t,f,s,e,u))
print(dict(bad)); 
for x in ex_bad: print(x)
