"""Of the 63 stale overlap_group labels, how many change the partition of markers (vs pure relabelling)? Writes nothing to repo."""
import csv, itertools, sys
shipped=list(csv.DictReader(open("src/viralscan/data/gene_programs.tsv"),delimiter="\t"))
regen=list(csv.DictReader(open(sys.argv[1]),delimiter="\t"))
key=lambda r:(r["virus"],r["gene_id_bundled"],r["programme"])
A={key(r):r for r in shipped}; B={key(r):r for r in regen}
n=0; changed=[]
for v in sorted({k[0] for k in A}):
    ks=[k for k in A if k[0]==v]
    for x,y in itertools.combinations(ks,2):
        sa=A[x]["overlap_group"]==A[y]["overlap_group"]; sb=B[x]["overlap_group"]==B[y]["overlap_group"]
        if sa!=sb: changed.append((v,x[1],y[1],sa,sb))
print("marker pairs whose same-group status differs:",len(changed))
for c in changed: print(*c,sep="\t")
