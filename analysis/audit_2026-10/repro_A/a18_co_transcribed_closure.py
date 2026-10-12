import sys; sys.path.insert(0,"extras")
import build_gene_programs as G
rec=G.parse_gtf("src/viralscan/data/Human_herpesvirus_8_NC_009333.gtf")
grp=G.build_overlap_groups(rec)
names=G.CO_TRANSCRIBED["Human herpesvirus 8"][0]
ids=sorted({g for n in names for g in G._resolve(rec,n)},key=lambda g:grp[g])
print("co-transcribed ids:",ids,"groups before:",[grp[i] for i in ids])
for i in ids:
    mates=[g for g in rec if grp[g]==grp[i] and g not in ids]
    print(i,"original group-mates NOT in the co-transcribed set (left behind in old group):",mates)
