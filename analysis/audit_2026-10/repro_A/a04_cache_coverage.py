import glob
from pathlib import Path
from viralscan.scripts.ncbi_fetch import genbank_cache_path
cache = Path.home()/".cache/viralscan/ncbi"
miss=[]; seqs=set()
for p in glob.glob("src/viralscan/data/*.gtf"):
    for l in Path(p).read_text().splitlines():
        if l and l[0]!="#": seqs.add((l.split("\t",1)[0], Path(p).name))
byseq={}
for s,f in seqs: byseq.setdefault(s,f)
for s,f in sorted(byseq.items()):
    try: fp=genbank_cache_path(s,cache)
    except Exception as e: miss.append((s,f,"ERR "+type(e).__name__)); continue
    if not fp.is_file(): miss.append((s,f,"nofile"))
print(len(byseq),"seqnames;",len(miss),"without cached flatfile"); print(miss[:10])
