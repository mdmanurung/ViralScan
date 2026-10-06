# F-025: The HHV-6B call in the EBV sample is human telomere-repeat reads, not HHV-6B

**Status:** confirmed (read-level). 2026-10-05.
**Sample:** SRR12682296 (EBV, 10x v2), combined run `ebv_latest_ref_2026-09-27/runs/combined/SRR12682296`.
**Work dir:** `/exports/para-lipg-hpc/mdmanurung/ViralScan/hhv6b_in_ebv_2026-10-05/` (`run.sh`, Slurm job 25701629).

## Claim
ViralScan reports `NC_000898.1_gene1` (HHV-6B, whole-genome placeholder) at 8 molecules in 6 cells. The
bus file holds 9 barcode-UMI records (14 reads) in one EC, 109549, which is HHV-6B only. Every read that
can carry that EC is a (TTAGGG)n / (CCCTAA)n telomere repeat. Competitive minimap2 against GRCh38 + the
2,215-genome panel: no read aligns better to HHV-6B than to human. The best hits are subtelomeres
(chr1:10 k, chr5:10-11 k, chr18:80.26 M, chr22:50.8 M) or interstitial telomere-like loci (chr8:205 k,
chr9:2.8 M, MAPQ 60). Where HHV-6B ties, it is only on its DR telomere arrays (61-161, 8.2-8.6 k,
153.4 k, 161.5-162 k). No read lands on HHV-6B outside those arrays.

## Why it happens
HHV-6A/6B carry perfect and imperfect TTAGGG arrays at the ends of their direct repeats (the sequence
they use to integrate into human telomeres). GRCh38 subtelomeres are absent from the cDNA host
reference, so kallisto has no host k-mer to compete and the reads land in the HHV-6B-only EC.
The DEF-01 `artefact` read filter keeps all 23 extracted pairs: TTAGGG repeats pass its complex-body check.

## Implications
- Do not report HHV-6B in SRR12682296. No virus beyond EBV has a read-validated call in that sample.
- Same mechanism can inflate HHV-6A/6B (and any panel genome with telomere-like arrays, e.g. HHV-7)
  in every sample; the small HHV-6A "bleed" in SRR20710641 is untested for it.
- Candidate fixes: mask TTAGGG arrays in the panel, or add subtelomeres to the host D-list.
