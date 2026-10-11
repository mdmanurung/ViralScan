"""Quantify overlap-group chaining vs a minimum shared-bp threshold (exon-overlap on the NORMALISED index).
Pairwise: genes i,j are linked if any block pair intersects >= T bp; groups = connected components.
T=1 reproduces build_overlap_groups (checked below)."""
import sys, collections, itertools
sys.path.insert(0,"extras")
import build_gene_programs as G
rep = __import__("viralscan.repeat_merges", fromlist=["x"]).load_repeat_merges()
def comps(records, T):
    ids=list(records); parent={g:g for g in ids}
    def f(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    # sweep: all blocks sorted; compare overlapping pairs
    bl=sorted((s,e,g) for g,r in records.items() for s,e in r["blocks"])
    active=[]
    for s,e,g in bl:
        active=[a for a in active if a[1]>=s]
        for s2,e2,g2 in active:
            if g2!=g and min(e,e2)-max(s,s2)+1>=T: parent[f(g)]=f(g2)
        active.append((s,e,g))
    return {g:f(g) for g in ids}
out=[]
Ts=(1,15,31,60,91)
cur_programme={}
for virus,gtf in sorted(G.VIRUS_GTF.items()):
    rec=G.parse_gtf(f"src/viralscan/data/{gtf}")
    for c,s in rep.items():
        if c in rec and s in rec:
            g=rec.pop(c); rec[s]["blocks"]|=g["blocks"]
    cur=G.build_overlap_groups(rec)
    c1=comps(rec,1)
    # sanity: T=1 partition == build_overlap_groups partition
    inv=lambda m: sorted(sorted(v) for v in __import__("itertools").groupby  and collections.defaultdict(list,{}).values()) 
    part=lambda m:{frozenset(k for k in m if m[k]==v) for v in set(m.values())}
    assert part(cur)==part(c1),(virus,"T=1 differs from build_overlap_groups")
    marks=collections.defaultdict(set)
    for e in G.CURATED_MARKERS:
        if e["virus"]==virus:
            for gid in G._resolve(rec,e["refseq_gene"]): marks[e["programme"]].add(gid)
    for T in Ts:
        m=comps(rec,T); sizes=collections.Counter(m.values())
        lat={m[g] for g in marks["latent"]}; prod={m[g] for g in marks["productive"]}
        # markers sharing a group with a marker of the other programme
        cross=len(lat&prod)
        allm=[g for s in marks.values() for g in s]
        grp=collections.Counter(m[g] for g in allm)
        nonov=sum(1 for g in allm if grp[m[g]]==1)
        out.append((virus,T,len(rec),len(sizes),max(sizes.values()),len(allm),nonov,cross))
print("virus\tT_bp\tgenes\tgroups\tlargest\tmarkers\tmarkers_alone_in_group\tgroups_with_both_programmes")
for r in out: print(*r,sep="\t")
