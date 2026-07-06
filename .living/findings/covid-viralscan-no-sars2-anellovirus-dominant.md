# covid_viralscan survey — no SARS-CoV-2, anellovirus dominant

**ID**: F-005
**Status**: ✅ SARS-CoV-2=0 CONFIRMED | ❌ **ANELLOVIRUS MAGNITUDE CLOSED: host-homology artifact** (2026-07-06)
**Date**: 2026-07-02 (invalidated→re-run→validated); 2026-07-03 anello magnitude under review; 2026-07-06 DECISIVE VERDICT

## ❌ VERDICT: ANELLOVIRUS SIGNAL IS HOST-HOMOLOGY ARTIFACT (2026-07-06)

Read-origin decisive test (job 25151971, commit `fee3397`): aligned 5M x213-g R2 reads to the
combined GRCh38+anellovirus STAR genome (`combined_GRCh38_2024A_serratus_plus_anellovirus`,
2313 viral contigs, `--outFilterMultimapNmax 50`).

**Result: viral-primary reads = 0 / 4,500,299 total primary-aligned (0.0%)**

Every read that ViralScan assigns to anellovirus has its STAR primary alignment on GRCh38, not
on any viral contig. This means the reads originate from **GRCh38 non-coding/intronic/intergenic
sequence that shares homology with the anellovirus panel** — regions absent from the cDNA-only
host reference, so kallisto/kb cannot suppress them as host reads.

**Consequence:**
- The ~90% Alphatorquevirus prevalence ("~90% of real cells at ≥5 UMI") is **NOT genuine
  anellovirus infection** — it is a cDNA-reference-homology artifact.
- **Do NOT cite the 90% figure in the manuscript.** Remove any pending TTV paragraph.
- The bulk GSE128078 pilot (job 25151978) confirms the same pattern: total_viral_rpm ~900,000
  per sample (90% of reads "viral") with Alphatorquevirus dominant — same artifact in bulk RNA-seq.
- SARS-CoV-2 = 0 stands and is unaffected. The negative control (SARS-CoV-1 = 0) also stands.

**Root cause**: ViralScan uses a cDNA-only host reference. Reads from GRCh38 non-coding regions
(introns, intergenic) that share sequence similarity with viral references are not counted as
host-mapping and appear as viral signal. `--multimap-method host-conservative` cannot correct
this because the host cDNA does not span those non-coding regions.

**Note on job exit code**: job 25151971 failed (exit 1:0) because `samtools view` exits with
code 1 when it cannot add a PG line due to duplicate `NC_002076.2` in the BAM header (known
dedup issue from the reference build). The analysis output is valid; the STAR run and awk
analysis completed and printed results before the bash `set -e` triggered on the samtools exit.

## ⚠️ ANELLOVIRUS MAGNITUDE (PRIOR REVIEW NOTE, 2026-07-03 — now resolved above)

**SARS-CoV-2 = 0 stands** (single genome, confirmed in both kb and STARsolo). But the
**anellovirus magnitude** ("~90% of real cells", 1.16M UMI) is UNDER REVIEW pending a
read-level check. A STARsolo combined-ref cross-check (job 25140486) gave **0 viral UMI**
for the whole panel — but that is a **GTF artifact, NOT evidence against ViralScan**: the
anellovirus gene records in `viral_genome.gtf` have gene+CDS but **no `exon`** records
(only 2,292/4,650 gene_ids have exons), so STARsolo GeneFull builds no countable interval
for them (EBV counted fine in P22.6 because its GTF has exons). What STARsolo DID validate:
cell-calling agrees (STARsolo 19,920 ⊂ CellRanger 28,922 ≈ emptyDrops 30,849), 90.2% valid
barcodes + 79.7% genome mapping (confirms the whitelist fix; explains kb's 6.4% as cDNA-only).
Open question (decisive test = align a read subsample to a viral-only STAR index vs GRCh38):
do the anello reads map uniquely to viral contigs (→ real, ViralScan right) or also to GRCh38
(→ host homology, kb cDNA-only over-call)? Do not cite the ~90% number until resolved.

## ✅ RE-RUN VALIDATION (2026-07-02): corrected whitelist works end-to-end

Re-ran covid quant (job 25140008) with the correct barcode whitelist — CellRanger's raw
barcode universe (2,974,869 barcodes; raw R1 match 68.1% vs 0.4% for v3). The matrix is now
real, validated against CellRanger's called cells for x213 (same sample as the CellRanger run):

- **CellRanger's 28,922 real cells now overlap the ViralScan matrix at 100%** (was 0.5% broken).
- Per-barcode UMI recovered: x216 max 30,125 / 4,985 barcodes ≥1000 UMI (was max 5,311 / 195).
- **Viral signal concentrates in real cells**: among the 28,922 CellRanger cells, 99.7% are
  viral+ (≥1 UMI); among empty droplets only 13.3%. Real cells median 398 total UMI vs empty 1.
  → **The viruses ARE in non-empty droplets**, ~7.5× enriched vs empty. (≥1 UMI is permissive —
  much of the 99.7% is ambient anellovirus; meaningful per-virus rates need a ≥2–5 UMI threshold
  + proper cell-calling, the motivation for the report-both-denominators feature.)

**Net**: the F-005 SARS-CoV-2=0 result holds; the per-cell anellovirus story is now on a valid
matrix and answerable. Fix committed in `slurm_viralscan_quant.sh` (WHITELIST → CellRanger-derived).

### Denominator demonstration + emptyDrops validation (2026-07-03)

Corrected re-run (25140008) COMPLETE. `emptyDrops` (DropletUtils, isolated conda env
`viralscan_celltools`) on the x213 matrix: **30,849 cells** vs CellRanger's 28,922 — **81.4%
of CellRanger cells recovered, Jaccard 0.65** (divergence expected: 6.4% pseudoalignment → a
sparse partial matrix, so emptyDrops' knee shifts). Validates emptyDrops AND that CellRanger
cells are the better ground truth when available.

**The denominator swing (Alphatorquevirus / Torque teno, x213):**
| Denominator | ≥1 UMI | ≥2 UMI | ≥5 UMI |
|---|---|---|---|
| all 578,938 barcodes | 99.8% | 56.1% | **37.6%** |
| CellRanger cells (28,922) | 99.7% | 98.6% | **90.0%** |
| emptyDrops cells (30,849) | 99.8% | 99.2% | 93.6% |

→ Over real cells, **~90% carry a genuine TTV load (≥5 UMI)** — anellovirus is near-ubiquitous
in this patient's cells; the all-barcode denominator reads 37.6% (diluted by empties). Same
mechanism as the HSV-1 P22.5 "25× discrepancy." Corrected summary also recovered far more total
signal (Alphatorquevirus 1,167,103 UMI vs 57,138 pre-fix) and more viruses (HHV-6, CeHV, EBV
EBNA-2, molluscum). SARS-CoV-2 / SARS-CoV-1 still 0 on the corrected matrix.

## ⛔ INVALIDATION (2026-07-02): wrong 10x barcode whitelist  _(kept for the record — resolved above)_

The per-cell results below are **not trustworthy**. The covid quant used the 10x **v3**
whitelist (`10x_version3_whitelist.txt.gz`), but the data's barcodes do not match it:

- `bustools correct`: **96.5% of BUS records were "uncorrected"** (off-whitelist).
- ViralScan matrix: 163,203 barcodes but **median 1 UMI/barcode, max 5,311, only 195 > 1000 UMI**
  — essentially all empty droplets. CellRanger called **28,922 real cells** from the same FASTQs.
- Only **156/28,922 (0.5%)** CellRanger cells are in the v3 whitelist. Raw R1 barcodes match
  the v3 whitelist **0.4%** (RC 0%), and match **no** bundled ngs_tools whitelist > ~5%
  (best: `10x_version4`/GEM-X 4.7%), yet match CellRanger's called cells **56.5%** directly.

Conclusion: this library's chemistry (likely **GEM-X 5′** or another CellRanger-auto-detected
set) is not covered by the v3 whitelist, so ~96% of reads were dropped and the resulting
per-cell matrix is noise. **Re-run required with the correct whitelist** (source the GEM-X/5′
barcode list CellRanger used, or pass CellRanger's whitelist to `viralscan -w`) before any
per-cell or "% infected" claim. Bulk-level "SARS-CoV-2 = 0 UMI" is the only relatively robust
takeaway, and even that should be re-confirmed post-fix.

---
_Original (now-invalidated) preliminary write-up follows:_

## Claim

The two 10x 5′ v3 GEX libraries `LUM-SJ-x213-g` and `LUM-SJ-x216-g`, scanned against a
combined human + SARS-CoV-2 + Serratus/anellovirus reference, show **no SARS-CoV-2
(NC_045512.2) signal** in either sample. The dominant viral signal is
**Alphatorquevirus (Torque teno virus / anellovirus)**.

## Evidence

`multimap_evidence.tsv` (host-conservative multimapping), both samples:

| Target | x213-g | x216-g |
|---|---|---|
| SARS-CoV-2 (NC_045512.2_gene1) | 0.0 UMI, 0 cells, `not_detected` | 0.0 UMI, 0 cells, `not_detected` |
| SARS-CoV-1 (sarsp1, neg. control) | 0.0 UMI, 0 cells, `not_detected` | 0.0 UMI, 0 cells, `not_detected` |

`viral_summary.tsv` — top detections:
- x213: Alphatorquevirus 57,138 UMI / 9,611 cells (5.9%); Beta/Gamma/Samektorque + Anelloviridae trace; MPXV_gp132 1 UMI.
- x216: Alphatorquevirus 85,471 UMI / 13,414 cells (6.6%); trace HHV-1/2/6 (≤3 UMI); MPXV_gp132 1 UMI.

The negative control being clean (SARS-CoV-1 = 0) alongside a strong, consistent
anellovirus signal indicates the specificity machinery is working, not that viral detection
is globally failing.

## Caveats (must resolve before firm conclusions)

- **Pseudoalignment is low: 6.4 % (x213) / 8.9 % (x216)** — consistent across samples;
  77–83 M reads aligned (usable) but the low rate is unexplained. Candidates: deep-library
  intronic/ambient fraction, chemistry/whitelist edge, D-list masking. No known-good baseline
  from the same index was locatable this session.
- **total_cells is high: 163,203 / 203,892 barcodes** at `min_counts=1000` on 1.2 B / 0.9 B
  read libraries — pulls in ambient barcodes, deflates `pct_infected`. Stage 4
  (`summarize_survey.py` CellRanger called-cell overlap) is needed to separate real cells from
  ambient before trusting per-cell fractions.
- Clinical expectation for these samples (were they expected COVID-positive?) is unconfirmed.

## Implications

If the libraries were expected COVID-positive, the low pseudoalignment + zero SARS-CoV-2
warrants a data/chemistry check before concluding true-negative. If they are general PBMC/
blood samples, anellovirus dominance with no SARS-CoV-2 is biologically unremarkable.

Related: [[raw-count-thresholds-confound-with-sequencing-depth]] (depth/threshold caveats),
pipeline bugs fixed en route logged in [[learnings.md]].
