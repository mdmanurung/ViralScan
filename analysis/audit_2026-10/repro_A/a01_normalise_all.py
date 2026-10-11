"""Run gtf_normalise over every bundled GTF; sanity-check output. Read-only."""
import re, sys, collections, glob, time
from pathlib import Path
from viralscan.gtf_normalise import normalise_gtf_file, normalise_viral_gtf
ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
cache = Path.home()/".cache/viralscan/ncbi"
issues = collections.defaultdict(list); rules = collections.Counter(); nf=0
t=time.time()
for p in sorted(glob.glob("src/viralscan/data/*.gtf")):
    nf+=1
    raw = Path(p).read_text().splitlines()
    try:
        new = normalise_gtf_file(p, cache)
        nocache = normalise_gtf_file(p, "/nonexistent_cache")
        twice = normalise_viral_gtf(new)
    except Exception as e:
        issues["EXC"].append((p, repr(e))); continue
    if twice != new: issues["not_idempotent"].append(p)
    # raw lines all preserved (modulo retag)?
    rawset = [l for l in raw]
    newset = set(new)
    # compare cols 0-7 of every original non-comment row
    sig = lambda l: tuple(l.split("\t")[:8])
    ns = collections.Counter(sig(l) for l in new if l.count("\t")==8)
    for l in raw:
        if l.count("\t")==8 and sig(l) not in ns: issues["raw_row_lost"].append((p,l[:80]))
    # structural
    tx = {}; exons = collections.defaultdict(list); genes = {}; txgene = {}
    for l in new:
        if l.startswith("#") or not l: continue
        c = l.split("\t")
        if len(c)!=9: issues["bad_cols"].append((p,l[:80])); continue
        s,e = int(c[3]),int(c[4])
        if s<1 or s>e: issues["coord"].append((p,l[:100]))
        a = dict(ATTR.findall(c[8]))
        if c[2]=="gene": genes.setdefault(a.get("gene_id"),[]).append((c[0],c[6],s,e))
        if c[2]=="transcript":
            if a["transcript_id"] in tx: issues["dup_tx_id"].append((p,a["transcript_id"]))
            tx[a["transcript_id"]] = (c[0],c[6],s,e,a.get("gene_id")); rules[a.get("viralscan_norm","passthrough_tx")]+=1
        if c[2]=="exon":
            exons[a.get("transcript_id")].append((s,e,a.get("exon_number"),c[6],a.get("gene_id")))
        if c[2] in("transcript","exon") and not a.get("gene_id"): issues["no_gene_id"].append((p,l[:80]))
    for t_,ex in exons.items():
        if t_ not in tx and not t_ is None: issues["exon_without_tx_row"].append((p,t_))
        srt = sorted(ex)
        for (s1,e1,*_),(s2,e2,*_) in zip(srt,srt[1:]):
            if s2<=e1: issues["exon_overlap"].append((p,t_,(s1,e1),(s2,e2)))
        strand = ex[0][3]
        order = [x[2] for x in ex]  # as written
        if None in order: continue
        nums = list(map(int,order))
        if nums!=list(range(1,len(nums)+1)): issues["exon_num_seq"].append((p,t_))
        starts=[x[0] for x in ex]
        want = sorted(starts, reverse=(strand=="-"))
        if starts!=want: issues["exon_order_vs_strand"].append((p,t_))
        if t_ in tx and tx[t_][:2]!=("x","x"):
            _,_,ts,te,_=tx[t_]
            if min(starts)!=ts or max(x[1] for x in ex)!=te: issues["tx_span_ne_exons"].append((p,t_))
    # every gene row yields >=1 transcript with exons
    hastx = collections.Counter(v[4] for v in tx.values())
    for g in genes:
        if g not in hastx: issues["gene_no_tx"].append((p,g))
    # exons within gene span
    for t_,(sq,st,ts,te,g) in tx.items():
        if g in genes and not any(ts>=gs and te<=ge for (_,_,gs,ge) in genes[g]): issues["tx_outside_gene"].append((p,t_))
    # no-cache vs cache difference
    if nocache!=new: issues["cache_changes"].append(p)
print("files",nf,"secs",round(time.time()-t))
print("rules",dict(rules))
for k,v in issues.items(): print(k,len(v),v[:4])
