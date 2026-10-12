# HSV-1 LAT is already exons-only in the index; the intron record is dropped and panel_completeness stays partial (CAT-29)

Status: confirmed · Date: 2026-10-12 · PLAN: CAT-29, CAT-10

## Claim
User decision 2026-10-12: re-cut `HSV1-LATonly.fasta` to its two exons and drop the 1,956 nt intron record. The re-cut changes nothing in the shipped index, because the index already carries LAT as the spliced two-exon transcript. HSV-1 `panel_completeness` stays `partial` and `latency_observable_in_rna` stays `false`.

## Evidence
- The three vendored records (`extras/vendor_sources/hsv/reference-genomes/HSV1-LATonly.fasta`, sha256 `f47f3cfb...fd89`, gitignored) are byte-identical to the NC_001806.2 slices 118805-119464 (660), 119465-121420 (intron, 1,956) and 121421-127151 (5,731).
- The bundled GTF models `HUM_HERP1_HHV1gp00s02` as two exons, 118805-119464 and 121421-127151 (no intron row). The built panel (`viral_ref_final/build/cdna.fa`, `HUM_HERP1_unassigned_transcript_62`) is **6,391 nt = 660 + 5,731**, so the indexed record contains no intron. `panel.t2g` lists the genomic span 118805-127151 only as coordinates.
- Nothing builds from the vendored FASTA (`extras/vendor_sources/VENDOR_MANIFEST.tsv`, VENDOR_REPORT section 6). It stays byte-identical to upstream; no derived FASTA was written.
- Canonical 31-mers (the `panel_kmer_sharing.canonical_kmers` encoding), 79 HSV-1 genes of the bundled GTF (CDS rows used where a gene has no exon row), TRL copy `s01` excluded:

| k-mer set | distinct 31-mers | in UL9/UL41/UL44/UL54 (the 4 productive markers) | in ICP0 (`p17`, `p76`) | in ICP34.5 (`p16`, `p77`) | in 465,769 human cDNA transcripts |
|---|---:|---:|---:|---:|---:|
| exons (660 + 5,731) | 5,765 | 0 | 1,879 (32.6 %), 1,852 (32.1 %) | 700 (12.1 %) | 0 |
| intron (1,956) | 1,825 | 0 | 716 (39.2 %) | 0 | 0 |

  - Exon 1 holds 630 and exon 2 holds 5,135 of the exon 31-mers; 30 more span the splice junction.
  - The intron's 1,825 31-mers occur in no exon, so dropping the record removes 1,825 independent LAT-intron k-mers and nothing else. They would count as HSV-1 latency signal while being unspliced, non-polyadenylated sequence that ICP0 antisense reads share by canonical k-mer.
  - Outside the terminal and internal repeats (positions 1-9,300 and 117,000-128,000 masked), 595 of the exon 31-mers occur elsewhere in the genome; none of the intron's do.
- The ICP0 and ICP34.5 sharing is the known antisense overlap of LAT with those genes, not a new overlap: ICP0 and ICP34.5 are not catalogue markers, and the repeat-merge rule (`repeat_gene_merges.tsv`) already folds TRL `s01` into `s02`.

## Re-evaluation of `panel_completeness`
The re-cut adds no second latency anchor and, if anything, removes the stable LAT introns that dominate 3' scRNA-seq LAT signal (Ouwendijk et al. 2024, PMID 38445887). So the shipped catalogue wording is unchanged: `gene_programs.tsv` HSV-1 `LAT` row says the indexed record is the spliced transcript, which excludes the intron; `partial`, `latency_observable_in_rna=false`, `latent` unreachable. Pinned by `tests/test_gene_programs_capturability.py::test_hsv1_lat_is_exons_only_and_the_sole_latency_anchor`.

## Caveats
- The host check is against Ensembl GRCh38 cDNA only. A genome-wide check (the F-005 lesson: non-coding and intergenic host sequence is absent from a cDNA index) was not run; this is a necessary, not a sufficient, host screen.
- "Tier 2" has no catalogue field yet: the `tier` and `persistence_class` columns of `virus_catalog.tsv` are empty on all 4,437 rows (CAT-10 is open). HSV-1 is in CAT-10's named set of nine human herpesviruses; set the column there, not for this row alone.
- k = 31 exact matches only, which is what kallisto uses to assign reads.

Tags: hsv-1, LAT, CAT-29, kmer, panel_completeness
