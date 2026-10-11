import sys
from viralscan.runconfig import RunConfig
out=sys.argv[1]; extra=dict(a.split("=",1) for a in sys.argv[2:])
base=dict(output=out,index="i.idx",transcripts="t.txt",sample1="s1.fq",sample2="s2.fq",cores=2,gtf=None,fasta=None,visual=True,f1=None,reference=False,umap=False,technology="10xv3",whitelist=None,multimapping=True)
for k,v in extra.items(): base[k]={"true":True}.get(v,v)
RunConfig.from_snakemake_config(base).to_yaml(out+"config.yaml") if __import__("os").makedirs(out,exist_ok=True) is None else 0
