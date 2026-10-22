# F-027: HHV-6B reads in the positional 3' flank class are not indexed by a CDS-only kallisto reference

**Status:** confirmed for one library (read-level, STAR alignment); preliminary as a general statement. 2026-10-22.
**Library:** GSM6416022 (HHV-6B, 10x 5' v2, reverse-stranded: host R2 reads are antisense, so `--rf-stranded` is the correct kallisto setting). Note: the "3' v2" comment inside the audit script is stale; the library is 5'.
**Reference:** panel viral contig `AF157706.1` (HHV-6B, 103 annotated `gene` rows). STAR alignments to GRCh38 + viral panel, primary reads only.

## Claim
The viral GTF supplies HHV-6B as gene-only rows with no exon rows and no 3' extension, so a kallisto index
built from CDS/exon sequence contains no sequence downstream of an annotated gene end. In GSM6416022,
8.9% of the viral reads that STAR places on HHV-6B fall in the positional 3' flank class (within 1 kb past an
annotated gene end, on the RNA strand of that gene), and the prominent coverage peak at about 125 kb (the
region downstream of U84/U85) lies in that class. STAR sees about 100x more reads there than minimap2 did.
Reads from this region are outside the indexed targets, so a CDS-only index cannot count them.

## Evidence
Read classes on `AF157706.1` (1,000 nt flank; RNA strand taken as opposite the read strand):

| Class | Reads | Share |
|---|---:|---:|
| Viral primary reads on AF157706.1 | 28,246,464 | 100% |
| Sense, overlapping an annotated gene | 19,320,601 | 68.4% |
| Antisense to an annotated gene | 5,440,651 | 19.3% |
| Positional 3' flank class | 2,504,792 | 8.9% |
| Intergenic (other) | 980,420 | 3.5% |

Other viral contigs in the same library are negligible (HHV-6A `NC_001664.4`: 21,172 reads; HHV-7: 149).

Index validation on the 125 kb reads (Track B pilot, flank index v2):

| Index | 125 kb reads pseudoaligned (of 2,849) |
|---|---:|
| Original CDS/gene-only index | 19 (0.7%) |
| First flank index (genes without exon rows not extended) | 19 |
| Flank index v2 (extension fixed) | 2,809 (98.6%) |

Provenance (not committed; local work area): `star_hostvirus/flank_audit/GSM6416022_classes.tsv`
(produced by `flank_audit.py`), `trackB_explore_2026-10/REPORT_trackB_explore_pilot.md` section 5,
`trackB_explore_2026-10/ref_flank_v2/` and its `validate_v2/` output.

## Scope
One library, one chemistry (10x 5' v2), one virus (HHV-6B, a library with a high viral read count).
The class labels depend on the annotation's gene coordinates and the 1 kb window chosen by the audit script.

## What this does NOT show
- It is not a measured loss fraction. 8.9% of viral reads fall in the flank class; that is a read
  classification, not the fraction of molecules or cells a CDS-only index would lose (flank reads might share a UMI with indexed reads -- a hypothesis, not measured -- and k-mer
  pseudoalignment of flank reads was not scored per cell). Class precedence is sense > flank > antisense:
  a read overlapping an annotated gene on its strand counts as sense, so flank reads lie wholly past the gene end.
- The downstream interval is not a validated UTR. It is positional (read lies past an annotated end on the gene's
  strand); no transcript-end, polyadenylation or 3' UTR annotation was used. It may also carry unannotated
  transcripts or readthrough.
- Nothing is inferred for other libraries, other chemistries or other viruses; the EBV library (GSM4796271) audit
  is a separate run.
- The antisense class (19.3%) and the intergenic class (3.5%) are not explained here.

## Implications
- A kallisto-only negative cannot exclude transcripts that lie outside the indexed targets. Report kallisto
  negatives as scoped to the indexed targets.
- Anellovirus scoped verdicts are STAR-based and are unaffected by this index gap.
- The Track A HHV-6B "scoped NEGATIVE" verdict is unchanged by this finding.
- Gene-only GTF rows should be extended (exon row plus flank) before index building; the fixed
  `extend_gtf` behaviour is what raised recovery from 19 to 2,809 reads.

## Comparison with Lareau et al. 2023 (Nature 623:608-615)
Source: public repository `caleblareau/hhv6-reactivation`, folder `hhv6-reference/`, pinned at commit
`dd2eb77d38527117bac7b1aab17f96b71234c886` (MIT licence; paper DOI 10.1038/s41586-023-06704-2). Files compared
(fetched 2026-10-22, sha256 recorded locally): `HHV6b_only.index.txt`, `AF157706_sequence.txt`,
`HHV6b_expression_annotations.tsv`.
- `HHV6b_only.index.txt`: 103 entries, all `AF157706.1_cds_<protein_id>_<n>` mapped to a gene label (DR genes
  appear twice, once per direct repeat).
- `AF157706_sequence.txt`: 103 FASTA records, all `gbkey=CDS` (NCBI CDS records of AF157706.1; 13 are spliced
  `join(...)` locations). No UTR, flank or genome sequence; no sign of masking (DR genes are ordinary entries).
- Coordinates: all 97 distinct CDS spans in ViralScan's HHV-6B annotation occur among the Lareau CDS records
  (same NCBI coordinates). ViralScan's gene rows end at or after the CDS end, so every read in the positional
  3' flank class here (past the ViralScan gene end) also lies outside every Lareau target.
Conclusion: the Lareau reference has the same blind spot as a CDS-only kallisto index: reads past annotated
coding ends, including the 125 kb peak region, cannot be represented. Not checked: the exact index build
command, the strand setting used in that study, and whether AF157706.1 is described there as strain Z29.
