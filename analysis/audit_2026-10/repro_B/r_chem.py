import random, tempfile
from viralscan import chemistry as c
random.seed(1)
bcs=["".join(random.choice("ACGT") for _ in range(16)) for _ in range(200)]
td=tempfile.mkdtemp()
open(td+"/plain.txt","w").write("\n".join(bcs)+"\n")
open(td+"/cr.txt","w").write("\n".join(b+"-1" for b in bcs)+"\n")
seqs=[random.choice(bcs)+"".join(random.choice("ACGT") for _ in range(10)) for _ in range(500)]
for n in ("plain","cr"):
    r=c.match_rates([seqs],{"user":f"{td}/{n}.txt"})[0]
    d=c.classify(seqs,r,f"{td}/{n}.txt")
    print(n,r,d.chemistry,d.reason[:70])
# empty sample
print(c.classify([],{ "10xv2":0.0},None).reason[:80])
print(c.read_structure([]))
