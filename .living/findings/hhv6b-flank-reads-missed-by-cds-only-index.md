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
One library, one chemistry (10x 5' v2), one virus (HHV-6B, a high-abundance lytic/reactivated sample).
The class labels depend on the annotation's gene coordinates and the 1 kb window chosen by the audit script.

## What this does NOT show
- It is not a measured loss fraction. 8.9% of viral reads fall in the flank class; that is a read
  classification, not the fraction of molecules or cells a CDS-only index would lose (reads in the flank can
  share a UMI with indexed reads, and k-mer pseudoalignment of flank reads was not scored per cell).
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
**Not completed from the primary reference file.** The file `HHV6b_only.index.txt` could not be located under
the user archive (bounded `find`, depth 6-7, several roots, no match), so it was not characterised directly.
The only record is a handoff note stating that the HHV-6B reference used by Lareau et al. is CDS-only with 103
`_cds_` entries; this is second-hand and unverified.
If that holds, the two references are alike in kind: both index coding sequence only, and ViralScan's own
HHV-6B annotation also has 103 gene rows, so an identical gene count is expected and says nothing about flank
content. Both would therefore share the same blind spot for reads past annotated gene ends. Differences
that remain to be checked against the file: whether any entry carries UTR/flank sequence, whether
`KDM2A`/DR cross-homology masking is applied (the paper's workflow filters these), and the strand setting used.
Action: locate the file (or the paper's supplementary reference) and replace this section with a measured comparison.
