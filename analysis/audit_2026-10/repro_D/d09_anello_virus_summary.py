"""anello_align.virus_summary: docstring says fractions are read-weighted so a 3-read accession
cannot outvote a 3000-read one; identity / query coverage are unweighted medians of per-accession
medians, so exactly that happens."""
from viralscan.anello_align import virus_summary

base = dict(unique_reads=0, start_sites=1, splice_reads=0, homopolymer_fraction=0, complex_body_fraction=1,
            reagent_fraction=0, r1_tso_fraction="", breadth=0.1, breadth_unique=0.1)
rows = [
    dict(accession="A", reads=3000, median_identity=0.99, median_query_coverage=0.98, **base),
    dict(accession="B", reads=3, median_identity=0.80, median_query_coverage=0.40, **base),
    dict(accession="C", reads=3, median_identity=0.80, median_query_coverage=0.40, **base),
]
out = virus_summary(rows, {}, {"A": "V", "B": "V", "C": "V"})["V"]
print("alignment_median_identity =", out["alignment_median_identity"], "(read-weighted mean would be ~0.9896)")
print("alignment_median_query_coverage =", out["alignment_median_query_coverage"], "(3000 of 3006 reads have 0.98)")
