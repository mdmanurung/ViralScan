"""Transitive closure vs direct overlap, per virus, for curated markers. Applies repeat merges but NOT CO_TRANSCRIBED."""
import sys, collections, itertools
sys.path.insert(0,"extras")
import build_gene_programs as G
rep=__import__("viralscan.repeat_merges",fromlist=["x"]).load_repeat_merges()
def ov(a,b): return max((min(e,e2)-max(s,s2)+1 for s,e in a for s2,e2 in b),default=0)
print("virus\tn_marker_pairs_same_group(chained)\tdirect>=1bp\tdirect>=31\tdirect>=91\tlargest_component_genes\tpairs_in_largest\tdirect_pairs_in_largest(>=1)\t(>=91)")
for virus,gtf in sorted(G.VIRUS_GTF.items()):
    rec=G.parse_gtf(f"src/viralscan/data/{gtf}")
    for c,s in rep.items():
        if c in rec and s in rec:
            g=rec.pop(c); rec[s]["blocks"]|=g["blocks"]
    grp=G.build_overlap_groups(rec)
    mk={}
    for e in G.CURATED_MARKERS:
        if e["virus"]==virus:
            for gid in G._resolve(rec,e["refseq_gene"]): mk[gid]=e["programme"]
    pairs=[(a,b) for a,b in itertools.combinations(sorted(mk),2) if grp[a]==grp[b]]
    d=[ov(rec[a]["blocks"],rec[b]["blocks"]) for a,b in pairs]
    sizes=collections.Counter(grp.values()); big=sizes.most_common(1)[0][0]
    mem=[g for g in rec if grp[g]==big]
    bp=[ov(rec[a]["blocks"],rec[b]["blocks"]) for a,b in itertools.combinations(mem,2)]
    print(virus,len(pairs),sum(x>=1 for x in d),sum(x>=31 for x in d),sum(x>=91 for x in d),len(mem),len(bp),sum(x>=1 for x in bp),sum(x>=91 for x in bp),sep="\t")
