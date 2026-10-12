"""EBV: which pairwise overlaps (bp) chain a productive marker to a latent marker, at T=1 and T=31/60."""
import sys, collections
sys.path.insert(0,"extras")
import build_gene_programs as G
rec=G.parse_gtf("src/viralscan/data/Epstein_Barr_virus_NC_007605.gtf")
rep=__import__("viralscan.repeat_merges",fromlist=["x"]).load_repeat_merges()
for c,s in rep.items():
    if c in rec and s in rec: rec.pop(c)  # simplified
prog={}
for e in G.CURATED_MARKERS:
    if e["virus"]=="Epstein-Barr virus":
        for g in G._resolve(rec,e["refseq_gene"]): prog[g]=e["programme"]
def edges(T):
    bl=sorted((s,e,g) for g,r in rec.items() for s,e in r["blocks"]); act=[]; E={}
    for s,e,g in bl:
        act=[a for a in act if a[1]>=s]
        for s2,e2,g2 in act:
            if g2!=g:
                o=min(e,e2)-max(s,s2)+1
                if o>=T: k=tuple(sorted((g,g2))); E[k]=max(E.get(k,0),o)
        act.append((s,e,g))
    return E
import itertools
for T in (1,31,61,91):
    E=edges(T); adj=collections.defaultdict(list)
    for (a,b),o in E.items(): adj[a].append((b,o)); adj[b].append((a,o))
    lat=[g for g,p in prog.items() if p=="latent"]; prd=[g for g,p in prog.items() if p=="productive"]
    found=False
    for s in prd:
        # BFS to any latent
        prev={s:None}; q=collections.deque([s])
        while q:
            u=q.popleft()
            if u in lat and u!=s:
                path=[]; x=u
                while x: path.append(x); x=prev[x]
                path=path[::-1]
                print(f"T={T}: {s.split('_',2)[-1]}(prod) -> {u.split('_',2)[-1]}(latent):", " -> ".join(f"{p.split('_',2)[-1]}" + (f"[{dict(adj[p])[path[i+1]]}bp]" if i+1<len(path) else "") for i,p in enumerate(path)))
                found=True; break
            for v,o in adj[u]:
                if v not in prev: prev[v]=u; q.append(v)
        if found: break
    if not found: print(f"T={T}: no productive-latent chain")
