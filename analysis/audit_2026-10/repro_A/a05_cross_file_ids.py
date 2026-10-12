import glob, re, collections
from pathlib import Path
from viralscan.gtf_normalise import normalise_gtf_file
ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
cache = Path.home()/".cache/viralscan/ncbi"
gid=collections.defaultdict(set); tid=collections.defaultdict(set); tg=collections.defaultdict(set)
txnoexon=[]; 
for p in sorted(glob.glob("src/viralscan/data/*.gtf")):
    new = normalise_gtf_file(p, cache); f=Path(p).name
    txs=set(); exs=set()
    for l in new:
        c=l.split("\t")
        if len(c)!=9: continue
        a=dict(ATTR.findall(c[8]))
        if c[2]=="gene": gid[a["gene_id"]].add((f,c[0]))
        if c[2]=="transcript": tid[a["transcript_id"]].add((f,c[0])); txs.add(a["transcript_id"]); tg[a["transcript_id"]].add(a["gene_id"])
        if c[2]=="exon": exs.add(a["transcript_id"])
    txnoexon += [(f,t) for t in txs-exs]
print("gene_id in >1 file/seq:", {k:v for k,v in gid.items() if len({x[0] for x in v})>1 or len({x[1] for x in v})>1} and len({k for k,v in gid.items() if len(v)>1}))
print("example", [(k,v) for k,v in gid.items() if len(v)>1][:3])
print("transcript_id in >1 file:", len([k for k,v in tid.items() if len({x[0] for x in v})>1]), [(k,v) for k,v in tid.items() if len({x[0] for x in v})>1][:3])
print("transcript_id -> >1 gene:", [(k,v) for k,v in tg.items() if len(v)>1][:5])
print("transcript rows without exon:", len(txnoexon), txnoexon[:5])
