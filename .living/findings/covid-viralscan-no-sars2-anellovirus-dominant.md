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

## ❗ UPDATE 2026-07-07: the `--genome-dlist` fix is LARGELY INEFFECTIVE (~15% removal)

Built a covid-matched genome-D-list index (`kallisto index --d-list genome.fa` on the existing
covid `cdna.fa` — identical 470,468 targets, D-list k-mers 621,183 → **2,809,623**) and
re-quantified both covid samples against it over the CellRanger cells (jobs 25167175/25167177).

**Result — anellovirus barely dropped:**
- x213-g: 975,301 → 835,971 UMI (**85.7% retained**); prevalence 99.68% → 99.48% of 28,922 CR cells.
- x216-g: 1,389,299 → 1,153,155 UMI (**83.0% retained**); prevalence 89.64% → 89.54% of 19,183 CR cells.
- Verified: right index used (config), SARS-CoV-2 stayed 0.

**Reconciliation with the read-origin test (NOT a contradiction):** STAR (mismatch-tolerant
alignment) puts 0/4.5M anello reads on viral contigs → reads are host. kallisto `--d-list` masks
only **exact** genome k-mers; imperfect host↔anellovirus homology means the anello-matching k-mers
a host read uses often are not *exactly* in GRCh38, so the D-list can't catch them. Exact-k-mer
masking removes ~15%; mismatch alignment removes ~100%.

**Consequences (revises the remediation plan):**
- The anellovirus signal remains a host-homology artifact (read-origin test is decisive and
  unchanged). Do NOT cite the ~90% figure.
- **`--genome-dlist` (the B5/bulk fix built + validated on disk) does NOT clean it** — only ~15%.
  The bulk GSE128078 (B5) analysis cannot rely on it. A real fix needs mismatch-tolerant host
  removal (align to full GRCh38 via STAR and drop host-origin reads) or excluding anellovirus
  from the panel. See [[learnings]] 2026-07-07.
- The genome-D-list index still helps other viruses whose host-homology is exact, and does no
  harm (SARS-CoV-2 still 0), but it is not the anellovirus remedy it was assumed to be.

## ✅ UPDATE 2026-07-07: STAR host-filter is the RELIABLE fix — removes ~95% of the artifact

Re-ran both covid samples with **`viralscan --host-filter starsolo --host-index
references/starsolo/human_GRCh38_2024A`** (STAR aligns to the full GRCh38 genome, mismatch-tolerant;
only unmapped reads reach the viral kallisto quant). Job 25175116 → `results_hostfilter/`.

**Result — anellovirus collapsed:**
- x213-g: 975,301 → **50,205 UMI (5.2% retained)**; prevalence 99.68% → 63.77% of CR cells. STAR
  flushed ~90% of reads as host (1.20 B input, only ~9.8% unmapped survived to viral quant).
- x216-g: 1,389,299 → **69,291 UMI (5.0% retained)**; prevalence 89.64% → 62.60%.
- SARS-CoV-2 = 0 throughout.

**This is the reliable-detection method** (mismatch-tolerant host removal, not exact-k-mer d-list):
STAR host-filter removes **95%** vs the d-list's 15%, at whole-dataset scale — confirming the
read-origin verdict AND the fix. Both the code path (`host_filter.py::_starsolo_filter`) and the
GRCh38 STAR index already exist; no code change was needed, the covid/bulk runs simply weren't using it.

**Residual (~5%, under characterization)**: 50k/69k UMI spread thinly over ~63% of cells (~2.7 UMI/cell).
Phase 2 = `viralscan evidence` coverage-breadth on the survivors (job 25181135) to decide genuine
low-level anellovirus (anelloviruses are ubiquitous in blood) vs residual artifact. Verdict pending.

**Recommended workflow / B5 correction**: use `--host-filter starsolo` (GRCh38 STAR genome) for
anellovirus and any homology-prone virus; drop the `--genome-dlist` plan for B5. See [[decisions]]
2026-07-07 (reliable anellovirus detection) and [[learnings]] 2026-07-07.

## ✅ UPDATE 2026-07-07 (Phase 2): coverage breadth shows the 5% residual is ALSO artifact — no real virus

Ran `viralscan evidence` on the host-filtered survivors (extract by CB/UMI → minimap2 → `samtools
coverage`; job 25181135 — NOTE: this was a manual remediation job after `viralscan evidence` job
25180994 crashed at `samtools sort` due to duplicate NC_002076.2 in the BAM header; the hf_align
script is not yet committed to `covid_viralscan/scripts/`). **Max coverage breadth: 3.41%
(KP343825.1, Gammatorquevirus, x216); most anellovirus contigs ≤2.5%; most other viral contigs
≤0.15%.** The residual reads pile deep+narrow on single loci — e.g. NC_001479.1 (EMCV, a
picornavirus — not anellovirus; likely GRCh38/IRES-homology locus) 114k/204k reads at
**1.99% breadth, 841×/1621× depth**. Tells:
- NC_001479.1 depth doubles (841× → 1621×) from x213 to x216 with **no increase in covered bases
  (156/7835 = 1.99107% in both)** → saturation at a fixed host-homologous locus, not infection.
- Anellovirus contigs (KP343825.1, MW455365.1, etc.) show 2.15–3.41% breadth — all far below
  the ≥10–20% expected for real infection of a 2.8 kb genome.
- The anellovirus-assigned reads largely **do not align to anellovirus genomes** under mismatch-tolerant
  minimap2 → at read level they are not anellovirus.

**Verdict — no genuine viral infection is supported in these samples.** Three converging methods
(two STAR-based + one independent minimap2 breadth analysis) agree: read-origin (0/4.5M on viral
contigs), STAR host-filter (95% removed), coverage-breadth (≤3.4% on any contig — all below the
≥5–10% threshold for real detection). SARS-CoV-2 = 0 remains the only trustworthy result.

✅ **Reproducibility note (2026-07-15 — RESOLVED)**: the decisive coverage.tsv was produced by
manual hf_align job 25181135 (not the documented `viralscan evidence` pipeline, which crashed at
`samtools sort` due to duplicate NC_002076.2). Fix committed 2026-07-15 as
`covid_viralscan/scripts/slurm_evidence_rerun.sh` (re-aligns existing `viral_reads.fasta` to
`viral_genome.dedup.fa`). The breadth figures above are now reproducible from the committed script.

## ✅ UPDATE 2026-07-15: CellTypist enrichment (T6) — anelloviruses in epithelial/plasma cells; EBV absent from B cells

CellTypist cell-type labels (PBMC immune atlas) + per-cell viral UMI from host-filtered results
(file: `results_hostfilter/celltypist_enrichment.tsv`).

**EBV (T6 tripwire)**: HHV4_EBNA-2 has **5 positive cells total** — all 5 in Epithelial cells
(OR = inf, p = 0.117, FDR = 1.0). Zero B cells (Memory B cells 0/1083; no other B cell
categories positive). **Conclusion: EBV is absent at biologically meaningful levels in this
COVID PBMC cohort.** The 5 putative EBV cells are noise/ambient, not B cell infection.

**Anellovirus cell-type bias (consistent with EVE-artifact hypothesis)**:

| Virus | OR in Epithelial cells | FDR |
|---|---|---|
| Alphatorquevirus | 3.92 | 0.0 |
| Betatorquevirus | 4.19 | 4e-102 |
| Samektorquevirus | 3.76 | 2e-05 |
| Gammatorquevirus | 3.25 | 4e-40 |

Additional: Alphatorquevirus enriched in **Plasma cells** (OR 3.08, FDR 5e-13). This is
consistent with the EVE-artifact mechanism: the known anellovirus EVE integrations are in
intronic regions of genes broadly expressed in epithelial and B-lineage tissues (e.g. NALCN/chr13,
LINC02742/chr11). Plasma cells are activated B cells with high transcriptional output — more
intronic pre-mRNA → more EVE reads passing the cDNA-only host filter.

HHV-1 (gp00p39): 4 positive cells, no significant enrichment. HHV-6b: 2 cells (summary). All
below the threshold for biological interpretation.

**Reliable-detection recipe (validated end-to-end)**: (1) `--host-filter starsolo` + GRCh38 STAR genome
(mismatch-tolerant host removal, ~95%); (2) require **breadth of coverage** via `viralscan evidence`
(a real call spreads across the genome, not a single homologous locus — this alone catches the artifact);
(3) anchor on SARS-CoV-2 = 0. Count-based UMI detection alone fabricated ~100%-of-cells prevalence.

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

→ ⛔ **SUPERSEDED** — pre-artifact-verdict claim: **~90% carry a genuine TTV load (≥5 UMI)**.
The 2026-07-07 host-filter + 2026-07-07 coverage-breadth analysis established these UMI are
artifact (host-homology reads, not anellovirus); do not cite this figure. See the ✅ UPDATE
2026-07-07 (Phase 2) section above and review 2026-07-15 findings F1–F4.
Anellovirus is near-ubiquitous in this patient's cells; the all-barcode denominator reads 37.6% (diluted by empties). Same
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
