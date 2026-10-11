import sys, itertools
sys.path.insert(0,"extras")
import build_gene_programs as G
rec=G.parse_gtf("src/viralscan/data/Epstein_Barr_virus_NC_007605.gtf")
rep=__import__("viralscan.repeat_merges",fromlist=["x"]).load_repeat_merges()
for c,s in rep.items():
    if c in rec and s in rec: g=rec.pop(c); rec[s]["blocks"]|=g["blocks"]
grp=G.build_overlap_groups(rec)
mk={}
for e in G.CURATED_MARKERS:
    if e["virus"]=="Epstein-Barr virus":
        for gid in G._resolve(rec,e["refseq_gene"]): mk[gid]=e["programme"]
ov=lambda a,b:max((min(e,e2)-max(s,s2)+1 for s,e in a for s2,e2 in b),default=0)
n=lambda g:g.split("_",2)[-1]
for a,b in itertools.combinations(sorted(mk),2):
    if grp[a]==grp[b]: print(grp[a],n(a),mk[a],"|",n(b),mk[b],"| direct overlap bp:",ov(rec[a]["blocks"],rec[b]["blocks"]))
