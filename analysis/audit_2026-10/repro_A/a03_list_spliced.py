import glob, re, collections
from pathlib import Path
from viralscan.gtf_normalise import normalise_gtf_file
ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
cache = Path.home()/".cache/viralscan/ncbi"
rows=[]
for p in sorted(glob.glob("src/viralscan/data/*.gtf")):
    new = normalise_gtf_file(p, cache)
    ex = collections.defaultdict(list); tag={}
    for l in new:
        c=l.split("\t")
        if len(c)==9 and c[2] in("exon","transcript"):
            a=dict(ATTR.findall(c[8]))
            if c[2]=="transcript": tag[a["transcript_id"]]=(a.get("viralscan_norm"),a["gene_id"],c[6])
            else: ex[a["transcript_id"]].append((int(c[3]),int(c[4])))
    for t,(r,g,s) in tag.items():
        if r and "spliced" in r:
            e=sorted(ex[t]); gaps=[e[i+1][0]-e[i][1]-1 for i in range(len(e)-1)]
            rows.append((Path(p).name[:40],t,r,s,len(e),gaps))
for r in rows: print(*r, sep="\t")
print(len(rows))
