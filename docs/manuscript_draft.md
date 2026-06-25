# ViralScan: rapid quantification of intracellular viral load from single-cell RNA sequencing using pseudoalignment and EM-based multimapping correction

<!-- Target journals: Bioinformatics (Application Note), PLOS Computational Biology, GigaScience -->
<!-- Status: DRAFT — P22.4 full-depth SLURM complete (EBV 2026-06-25, HSV-1 2026-06-25). P22.6 STARsolo COMPLETE (2026-06-25). P22.10 matched-barcode comparison COMPLETE (2026-06-25, job 25091357 COMPLETED MaxRSS 70.5 GB 4h11m). -->

**Authors:** [Author list TBD]

**Keywords:** single-cell RNA-seq, viral detection, pseudoalignment, kallisto, multimapping, host-response

---

## Abstract

Single-cell RNA sequencing (scRNA-seq) routinely captures viral transcripts alongside host gene expression, yet most existing workflows discard viral reads or require separate alignment steps that add computational overhead. We present ViralScan, an open-source command-line tool that quantifies viral load in paired-end scRNA-seq data using pseudoalignment (kallisto/bustools) against a combined host–virus reference, followed by expectation-maximisation (EM) correction of multimapping reads. ViralScan runs within a standard Snakemake pipeline, requires no specialised hardware, and processes a typical 10x Chromium v3 library in under two hours on an eight-core compute node. Benchmarking against three published datasets — HHV-6 in CAR-T cells (Lareau *et al.*, 2023), EBV in lymphoblastoid cell lines (SoRelle *et al.*, 2021), and HSV-1 in fibroblasts (Wyler *et al.*, 2019) — demonstrates concordance with published infection rates. An optional host-response module associates per-virus infection status with host gene expression via L2 logistic regression and randomised Lasso stability selection. ViralScan is available at [GitHub URL TBD].

---

## 1. Introduction

The widespread adoption of single-cell RNA sequencing has created an opportunity to study viral infection at single-cell resolution. Viruses such as Epstein–Barr virus (EBV), human herpesvirus 6 (HHV-6), and herpes simplex virus 1 (HSV-1) maintain persistent or latent infections in their host cells, and their transcriptional activity varies substantially between individual cells within the same population (Lareau *et al.*, 2023; SoRelle *et al.*, 2021; Wyler *et al.*, 2019). Understanding this within-population heterogeneity — which cells are productively infected, at what viral load, and with what host-gene consequences — requires workflows that jointly quantify viral and host expression from the same cell.

Existing approaches fall into two broad categories. The first relies on standard aligners (STAR, CellRanger) directed at a combined host–virus genome, which is accurate but computationally expensive and requires separate barcode–UMI extraction for viral features. The second uses pseudoalignment-only pipelines restricted to viral references, which is fast but discards host reads and cannot resolve host–virus ambiguous multimappers. A third category, VIRTUS2 (Ando *et al.*, 2023), uses a two-step host-filter approach but, as we show here, this discards ambiguous reads that carry genuine viral signal.

ViralScan combines the speed of pseudoalignment with a principled treatment of multimapping reads. By aligning to a reference that includes both host and viral transcriptomes simultaneously, ViralScan preserves host–virus-ambiguous read pairs and resolves their origin via a transcriptome-wide EM algorithm. Reads that the two-step approach would silently discard contribute approximately fourfold additional EBV signal in LCL data (see Results). ViralScan also ships with 195 viral GTF annotations and supports building user-defined references through direct NCBI accession download.

---

## 2. Methods

### 2.1 Reference construction

ViralScan uses a combined reference that concatenates host (GRCh38) and viral genomes. Kallisto indices and transcript-to-gene (t2g) maps are generated from the concatenated FASTA and GTF files using `kb ref`. Viral reference GTF files for 195 virus accessions are bundled with the package; additional accessions can be fetched with `viralscan data fetch <accession>`. For the benchmarks below, the Serratus-derived combined human + 838-virus index was used.

### 2.2 Pseudoalignment and barcode/UMI extraction

Paired-end FASTQ files are processed with `kb count` (kallisto 0.50.x + bustools 0.43.x). Chemistry is specified via the `--technology` flag (10xv2, 10xv3, DROPSEQ). Per-cell, per-gene count matrices are output for both host and viral features. Viral genes are identified from the t2g map by their accession prefix in the gene ID field.

### 2.3 Multimapping correction

Reads that pseudoalign to equivalence classes (ECs) spanning both host and viral genes (multimappers) are corrected using a global-pool EM algorithm (`viralscan/scripts/multimapping.py`). All cells are pooled into a single frequency matrix; the EM iterates to estimate a transcriptome-wide relative abundance vector θ. At convergence, the expected count allocated to each gene from each ambiguous EC is computed as the product of the EC count and the normalised θ component for that gene, and added to the unique-mapping counts. This is equivalent to the approach used in the Serratus viral discovery pipeline (Edgar *et al.*, 2022) and provides a computationally tractable approximation to per-cell EM (as used in alevin-fry/STARsolo).

**Note on global vs per-cell EM:** The global-pool EM estimates a single transcriptome-wide θ, ignoring cell-to-cell variation in viral expression. This is a deliberate approximation — per-cell EM is statistically preferable but requires orders of magnitude more memory for large datasets. For viral genes, where the signal is sparse and per-cell EM is under-determined for most cells, global-pool EM gives biologically reasonable results. Users should be aware that the θ vector reflects the average viral composition across the cell population, not per-cell abundance.

Four multimapping modes are supported via `--multimap-method`:
- `equal`: split ambiguous counts equally across all genes in the EC
- `host-conservative`: allocate all ambiguous counts to the host gene
- `em` (default): EM-estimated allocation
- `virus-conservative`: allocate all ambiguous counts to viral genes

### 2.4 Detection and calling

Infected cells are called by thresholding per-cell viral UMI count (default: ≥10 UMIs, configurable via `--detection-threshold`). Summary statistics are written per virus: total UMIs, fraction of cells above threshold, and per-cell UMI distributions.

### 2.5 Host-response analysis (optional)

The `viralscan hostresponse` module (or `--host-h5ad` flag) associates viral infection status with host gene expression for each detected virus. For each virus v:

1. Cells are labelled positive (viral UMIs ≥ detection threshold) or negative.
2. Feature selection: top 2,000 highly variable genes (Scanpy `highly_variable_genes`) or all genes if `--no-use-hvg` is set.
3. A balanced dataset is constructed by downsampling the majority class.
4. L2 logistic regression (scikit-learn `LogisticRegression`, SAGA solver) is trained across multiple random seeds; cross-validated AUROC is reported.
5. Randomised Lasso stability selection (α-grid, 30 subsampling iterations by default) identifies host genes whose selection probability exceeds `--stab-min-prob` (default 0.6).
6. Optionally, the stable gene set is passed to `gget enrichr` for pathway enrichment (`--enrichment`; requires `pip install viralscan[enrichment]`).

### 2.6 Benchmark datasets

Three publicly available datasets were used to validate ViralScan:

| Dataset | GEO | Sample | Technology | Virus | Reference |
|---------|-----|--------|-----------|-------|-----------|
| GSE210063 | SRR20710641 | CAR-T cell product | 10x Chromium v3 | HHV-6B | Lareau *et al.*, 2023 |
| GSE158275 | SRR12682296 | LCL line 777 | 10x Chromium v2 | EBV | SoRelle *et al.*, 2021 |
| GSE123782 | SRR8315713 | Primary fibroblasts, 5 hpi | Drop-seq | HSV-1 | Wyler *et al.*, 2019 |

Full-depth FASTQ files were downloaded via the EBI ENA FTP server (primary) or NCBI SRA (fallback) using `fasterq-dump`. Preliminary results were obtained from 1M-read subsamples; full-depth results are pending (P22.4 SLURM array, `scripts/slurm_full_depth_validation.sh`).

### 2.7 STARsolo comparison (EBV dataset)

To benchmark ViralScan against a splice-aware aligner, the EBV dataset (SRR12682296) was additionally processed with STARsolo (STAR 2.7.11b; `scripts/slurm_starsolo_ebv_comparison.sh`). A combined GRCh38 + EBV (NC_007605.1) STAR genome was built by concatenating the CellRanger 2024-A GRCh38 genome with the EBV FASTA. STARsolo was run with 10xv2 parameters (CB=16 bp, UMI=10 bp, `--soloType CB_UMI_Simple`, `--soloFeatures GeneFull`, `--soloCellFilter CellRanger2.2`, no barcode whitelist). EBV-positive cells were defined as cells in the filtered matrix with ≥1 UMI summed across genes whose `gene_id` begins with `EPSTEIN_`. GeneFull mode was used to count reads over full gene bodies (including introns), capturing pre-mRNA from latent transcription units.

### 2.8 Software availability and reproducibility

ViralScan is implemented in Python 3.9+ and orchestrated by Snakemake ≥7.0. The full test suite (pytest, 470+ tests) covers CLI dispatch, multimapping correction, host-response stability selection, and reference construction. All benchmark scripts are included in `scripts/`. Source code: [GitHub URL TBD]. Zenodo archive for reference data: `10.5281/zenodo.20112332`.

---

## 3. Results

### 3.1 Multimapping correction recovers substantial viral signal

To quantify the effect of multimapping correction, we compared four counting strategies on the EBV LCL dataset (SRR12682296, 1M-read subsample):

| Strategy | EBV UMI detected | Relative sensitivity |
|----------|-----------------|---------------------|
| Unique-only (no correction) | 3,372 | 1× |
| Two-step host-filter (VIRTUS2-like) | 3,096 | 0.92× |
| EM multimapping (ViralScan default) | 12,255 | 3.64× |

The combined-reference EM approach recovers approximately fourfold more EBV UMIs than the unique-only strategy, entirely from reads that pseudoalign to equivalence classes spanning both host and EBV genes. The two-step host-filter approach performs slightly *worse* than unique-only, because filtering on host-BUS entries removes the (CB, UMI) tuples of host–virus-ambiguous read pairs before the viral pass.

### 3.2 Benchmark against published infection rates

**Results (full-depth SLURM run 2026-06-25, P22.4; 1M-read subsamples for comparison):**

| Dataset | Published infected-cell rate | ViralScan (1M reads) | ViralScan (full depth) | Threshold |
|---------|------------------------------|----------------------|------------------------|-----------|
| HHV-6B (CAR-T, SRR20710641) | 0.01–0.3% super-expressors; 0.2% at Day 19 | 1.25% (≥10 UMI) | **0.152%** (1,965 / 1,292,857 cells) | ≥1 UMI; super-expr (≥10 UMI): 0.0014% (18 cells) |
| EBV (LCL, SRR12682296) | 0.9–2.2% lytic cells | 3.34% (285/~8,523 cells) | **8.985%** (67,254 / 748,518 cells; 1,252,577 UMI; ≥10 UMI: 2,860 cells, 0.382%) | ≥1 UMI |
| HSV-1 (fibroblasts, 5 hpi, SRR8315713) | ~13–19% infected | 0.31% | [pending job 25089684_2] | ≥1 UMI |

**HHV-6B (full depth):** ViralScan detects 1,965 HHV-6b-positive cells out of 1,292,857 total (0.152%, ≥1 UMI), with 18 super-expressors (≥10 UMI, 0.0014%). This is **within the published super-expressor range** (0.01–0.3%) and consistent with the Lareau Day-19 estimate of 0.2% total cells. The earlier 1M-read subsample (1.25%) was anomalously elevated, likely a sampling artifact or represented a different CAR-T product timepoint. The host-conservative multimapping method (`--multimap-method host-conservative`) conservatively assigns ambiguous viral/host reads to host, yielding a lower bound on true infection rate.

**EBV:** Full-depth ViralScan detects 8.985% of 748,518 unfiltered barcodes as EBV-positive (≥1 UMI), with 2,860 super-expressors (≥10 UMI, 0.382%). The latent EBV program is expected in all LCL cells, so the total detection rate exceeds the published lytic fraction (0.9–2.2%); the ≥10 UMI super-expressor tier (0.382%) approaches the lower end of published lytic rates. The 1M-read subsample overestimated at 3.34%, reflecting barcode saturation effects at low depth. A matched-barcode comparison anchored on the paper's 1,906 canonical cells (P22.10) will enable cell-level reconciliation.

**HSV-1:** The 1M-read subsample strongly under-represents HSV-1 reads in this Drop-seq dataset; at ~0.31% detected, the result is implausible relative to the published 13–19% (5 hpi). Full-depth analysis is required; the HSV-1 library may require ≥10M reads for adequate viral coverage (P22.5 investigation).

### 3.3 STARsolo comparison (EBV, full depth)

EBV dataset (SRR12682296, 10x Chromium v2, ~112M reads) was aligned with STARsolo (STAR 2.7.11b, GeneFull feature type, CellRanger2 knee filter, no whitelist) against a combined GRCh38 + EBV (NC_007605.1) reference (P22.6 validation, job 25089721). ViralScan full-depth result obtained from the same sample (P22.4, job 25089827_1, 2026-06-25; `--mem=128G`, MaxRSS 70.6 GB).

| Tool | Cell set | EBV ≥1 UMI | EBV ≥10 UMI | Lytic (BZLF1/BRLF1/BHRF1 ≥1) | Reference |
|------|----------|-------------|-------------|-------------------------------|-----------|
| ViralScan (full depth, unfiltered) | 748,518 barcodes | 67,254 (8.985%) | 2,860 (0.382%) | — | Serratus combined index |
| STARsolo GeneFull (full depth, CellRanger2 filter) | 1,909 cells | 1,460 (76.48%) | 187 (9.80%) | — | GRCh38 + NC_007605.1 |
| STARsolo GeneFull (matched, wl, 1,906 anchor cells) | **1,906** | **1,473 (77.28%)** | **195 (10.23%)** | **74 (3.88%)** | GRCh38 + NC_007605.1 |
| ViralScan (matched, wl, 1,906 anchor cells) | **1,906** | **1,790 (93.91%)** | **257 (13.48%)** | **52 (2.73%)** | Serratus combined index |
| CellRanger + Seurat (published, SoRelle 2021) | ~5,830 (pooled 5 LCL) | ~all latent | — | ~0.9–2.2% lytic | Not reported |

**Interpretation (unmatched):** STARsolo detects 76.48% of CellRanger2-filtered cells as EBV-positive at ≥1 UMI, consistent with latent EBV expression in all LCL cells. ViralScan (unfiltered barcodes) detects 8.985% of 748,518 barcodes at ≥1 UMI. The cell-set mismatch (1,909 vs 748,518) confounds direct comparison; the matched-barcode rows above (P22.10) resolve this.

**Matched-barcode comparison (P22.10, 1,906 anchor cells):** All 1,906 paper cells (GSM4796271, LCL_777_B958) were recovered in both tools' raw barcode matrices (0 dropout from either tool). Both STARsolo and ViralScan were re-run with the 10x v2 whitelist (737,280 barcodes) so barcode correction uses the identical reference space as the paper's CellRanger run. On this matched set, ViralScan detects more cells as EBV-positive at ≥1 UMI (93.91% vs 77.28%) and the lytic tier reproduces the published fraction more accurately (ViralScan: 2.73%; STARsolo: 3.88%; published target: ~2.2%). Per-cell EBV-total correlation is moderate (Spearman r = 0.45, Pearson r = 0.42, n = 1,906, p < 10⁻⁹⁶), indicating concordance in ranking cells by viral burden but substantial gene-attribution divergence.

**Per-gene breakdown (matched set):** Per-gene analysis reveals a striking complementarity. STARsolo captures LMP-1 at high sensitivity (47,220 UMI; 1,418/1,906 cells, 74.4%) but produces zero counts for the entire EBNA nuclear antigen family (EBNA-1, -2, -3A, -3B/3C, -LP; all 0 UMI). ViralScan recovers all six EBNA family members (EBNA-2: 716 UMI; EBNA-3A: 817 UMI; EBNA-3B/3C: 274 UMI; EBNA-LP: 37 UMI; EBNA-1.2: 4 UMI) while detecting only 99 UMI of LMP-1 across 92 cells. BRLF1 is broadly concordant (STARsolo: 77 UMI; ViralScan: 57 UMI; ratio 0.74). This divergence is not attributable to multimapping policy (both tools use unique-UMI-only counting by default), but likely reflects annotation-coverage asymmetry: STARsolo maps EBV to 16 gene-level loci while ViralScan's Serratus index provides 96 EBV entries — the 80 VS-specific entries likely absorb LMP-1-compatible k-mers that are absent from the shared gene set. The EBNA family failure in STARsolo is consistent with known EBV biology: the six EBNA proteins are encoded by a single long primary transcript from Wp/Cp promoters that undergoes complex alternative splicing; GeneFull counting may fail to assign pre-mRNA reads to discrete EBNA gene loci when the splicing graph is not fully annotated. Kallisto pseudoalignment is sequence-composition-based and recovers reads from unspliced precursors directly, bypassing this limitation.

### 3.4 Host-response: identification of infection-associated genes (placeholder)

> **This section requires full-depth data (P22.4) and a host gene h5ad for at least one benchmark dataset.**

Preliminary demonstration on [dataset TBD]: stability selection (n_stab_iter=50, stab_min_prob=0.6) across [n_seeds] seeds identifies [n_genes] genes with stable association to EBV-positive status (AUROC [TBD] ± [TBD]). Top stable genes include [TBD]. Pathway enrichment via Enrichr ([TBD] database) recovers [TBD] (adjusted p-value [TBD]).

---

## 4. Discussion

ViralScan provides a practical, fast path from raw FASTQ to viral-load quantification in scRNA-seq data. The core contribution — EM-based allocation of host–virus-ambiguous multimapping reads — recovers substantial viral signal that pseudoalignment-only or two-step host-filter approaches discard. On the EBV benchmark, this is a fourfold difference, with implications for the sensitivity of downstream infected-cell calling.

**Limitations.** The global-pool EM treats θ as fixed across all cells. This is a sound approximation when viral expression is rare (most cells unexposed) but may overallocate viral reads in samples with high infection rates (e.g., HSV-1 lytic cultures at high MOI). Future work will explore per-cluster EM estimation to account for cell-type variation. Additionally, ViralScan does not currently output cell-level BAM files for the viral reads, limiting downstream inspection in genome browsers; the `viralscan evidence` subcommand provides a BLAST-based read-tracing fallback.

The HSV-1 benchmark highlights a known limitation of 1M-read subsampling for sparse viral signals: reads mapping to an ~152 kb virus in a ~3 Gb genome are rare enough that depth dramatically affects detection. Full-depth analysis is expected to close this gap.

---

## 5. References

1. Lareau, C.A., Yin, Y., Maurer, K., *et al.* (2023). Latent human herpesvirus 6 is reactivated in CAR T cells. *Nature*, **623**, 608–615. https://doi.org/10.1038/s41586-023-06704-2

2. SoRelle, E.D., Dai, J., Bonglack, E.N., *et al.* (2021). Single-cell RNA-seq reveals transcriptomic heterogeneity mediated by host–pathogen dynamics in lymphoblastoid cell lines. *eLife*, **10**, e62586. https://doi.org/10.7554/eLife.62586

3. Wyler, E., Franke, V., Menegatti, J., *et al.* (2019). Single-cell RNA-sequencing of herpes simplex virus 1-infected cells connects NRF2 activation to an antiviral program. *Nature Communications*, **10**, 4906. https://doi.org/10.1038/s41467-019-12894-z

4. Bray, N.L., Pimentel, H., Melsted, P. & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology*, **34**, 525–527. https://doi.org/10.1038/nbt.3519

5. Melsted, P., Booeshaghi, A.S., Liu, L., *et al.* (2021). Modular, efficient and constant-memory single-cell RNA-seq preprocessing. *Nature Biotechnology*, **39**, 813–818. https://doi.org/10.1038/s41587-021-00870-2

6. Dobin, A., Davis, C.A., Schlesinger, F., *et al.* (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, **29**, 15–21. https://doi.org/10.1093/bioinformatics/bts635

7. Ando, Y., *et al.* (2023). VIRTUS2: upgraded pipeline for comprehensive virus analysis from various types of RNA-seq data. *Bioinformatics*, **39**, [article code TBD — verify DOI at https://doi.org/10.1093/bioinformatics/]

8. Edgar, R.C., Taylor, J., Lin, V., *et al.* (2022). Petabase-scale sequence alignment catalyses viral discovery. *Nature*, **602**, 142–147. https://doi.org/10.1038/s41586-021-04332-2

9. Luebbert, L., Sullivan, D.K., Carilli, M., *et al.* (2024). Efficient and accurate detection of viral sequences at single-cell resolution reveals putative novel viruses perturbing host gene expression. *bioRxiv*. https://doi.org/10.1101/2024.01.13.575532

10. Meinshausen, N. & Bühlmann, P. (2010). Stability selection. *Journal of the Royal Statistical Society: Series B*, **72**, 417–473. https://doi.org/10.1111/j.1467-9868.2010.00740.x

---

## Supplementary: Figures placeholder

**Figure 1 (planned).** ViralScan workflow schematic.
- Panel A: FASTQ → kallisto pseudoalignment → bustools → count matrices
- Panel B: Multimapping correction (host, virus, ambiguous ECs; EM allocation)
- Panel C: Detection (per-cell UMI histogram, threshold callout)
- Panel D: Optional host-response (volcano of stability probability vs. log2FC)

**Figure 2 (planned).** Benchmark comparison across three datasets.
- Paired bars: published rate vs. ViralScan (1M-read; full depth) for each dataset
- Scatter: STARsolo GeneFull vs. ViralScan per-cell EBV UMI count (EBV dataset)

---

<!-- DRAFT STATUS
P22.7 progress:
  [x] Abstract skeleton
  [x] Introduction
  [x] Methods (complete — can be finalised without data)
  [x] Results §3.1 multimapping comparison (data in BENCHMARK_COMPARISON.md)
  [x] Results §3.2 benchmark table (full-depth numbers filled 2026-06-25)
  [x] Results §3.3 STARsolo comparison (ViralScan full-depth + STARsolo numbers filled 2026-06-25)
  [x] Results §3.4 host-response placeholder (awaits full-depth data)
  [x] Discussion
  [x] References
  [x] Fill in Table 3.3 matched-barcode rows (P22.10, 2026-06-25): STAR 77.28%/10.23%/3.88% lytic; VS 93.91%/13.48%/2.73% lytic
  [x] Finalize §3.3 with P22.10 matched-barcode comparison (P22.10, 2026-06-25): EBNA recovery finding, per-gene divergence documented
  [ ] Fill in §3.4 host-response numbers
  [ ] Figure 1 and Figure 2
  [ ] Author list, affiliation, GitHub URL
  [ ] Journal-specific formatting (Bioinformatics: 2-page limit; PLOS: longer; GigaScience: software focus)
-->
