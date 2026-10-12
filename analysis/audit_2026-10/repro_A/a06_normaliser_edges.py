from viralscan.gtf_normalise import normalise_viral_gtf as N
def r(f,s,e,g,st="+",**a):
    x="".join(f' {k} "{v}";' for k,v in a.items())
    return f'S\tx\t{f}\t{s}\t{e}\t.\t{st}\t.\tgene_id "{g}";{x}'
def show(title, lines):
    print("##",title)
    try:
        for l in N(lines): print("  ",l.replace("\t"," | ")[:170])
    except Exception as ex: print("  EXC",type(ex).__name__,ex)
show("empty", [])
show("only comments", ["#x",""])
show("inverted gene 300..100", [r("gene",300,100,"G"), r("CDS",120,200,"G",protein_id="P")])
show("zero-start", [r("gene",0,100,"G")])
show("gene row w/o gene_id on CDS", [r("gene",1,100,"G"), 'S\tx\tCDS\t5\t50\t.\t+\t0\tprotein_id "P";'])
show("gap boundary 10 vs 11 (single protein)", [r("gene",1,300,"A"),r("CDS",11,100,"A",protein_id="P"),r("CDS",111,200,"A",protein_id="P"),
   r("gene",401,700,"B"),r("CDS",411,500,"B",protein_id="P"),r("CDS",512,600,"B",protein_id="P")])
show("two protein_ids at same CDS-rows, one missing protein_id", [r("gene",1,300,"A"),r("CDS",11,100,"A",protein_id="P"),r("CDS",201,260,"A")])
show("gene+transcript row but no exon", [r("gene",1,100,"G"), r("transcript",1,100,"G",transcript_id="T")])
show("CDS with no transcript_id, attr w/o trailing ;", ['S\tx\tgene\t1\t100\t.\t+\t.\tgene_id "G"', 'S\tx\tCDS\t5\t50\t.\t+\t0\tgene_id "G"'])
show("dup identical gene rows", [r("gene",1,100,"G"),r("gene",1,100,"G")])
show("unsorted input: CDS before gene", [r("CDS",11,100,"G",protein_id="P"),r("CDS",201,260,"G",protein_id="P"),r("gene",1,300,"G")])
show("gene_id with quoted semicolon/space", [r("gene",1,100,"G x;y")])
show("two genes same id on two strands", [r("gene",1,100,"G","+"),r("gene",1,100,"G","-"),r("CDS",5,50,"G","+")])
show("CDS outside gene row", [r("gene",1,100,"G"),r("CDS",500,600,"G",protein_id="P")])
