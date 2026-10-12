import collections
from viralscan import virus_catalog
from viralscan.virus_identity import assign_virus_keys, _Catalogue, _default_anello_genus
rows=virus_catalog.load_catalogue()
print("catalogue rows",len(rows), "cols", list(rows[0]))
ag=_default_anello_genus()
keys=assign_virus_keys(rows,ag)
c=collections.Counter(keys)
split=[k for k in c if "|" in k]
print("distinct keys",len(c),"split-by-strain keys",len(split))
bytax=collections.defaultdict(list)
for k in split: bytax[k.split("|")[0]].append((k.split("|",1)[1],c[k]))
for t,v in bytax.items(): print(t,v)
# duplicates of accession
acc=collections.Counter(r["accession"] for r in rows); print("dup accession rows",[a for a,n in acc.items() if n>1][:5])
av=collections.Counter(r["accession_version"] for r in rows); print("dup accession_version",[a for a,n in av.items() if n>1][:5])
# display names unique?
cat=_Catalogue(rows,ag); names=list(cat.names.values()); print("names",len(names),"unique",len(set(names)),"with [key] suffix",sum('[' in n for n in names))
# rows with empty taxid
print("rows no taxid",sum(1 for r in rows if not (r.get('taxid') or '').strip()))
