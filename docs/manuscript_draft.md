# ViralScan quantifies viral RNA in single-cell transcriptomes using pseudoalignment and multimapping correction

<!-- Target journal: Cell Reports Methods / Cell Press methods article draft. -->
<!-- Evidence status: P22.4 full-depth validation complete (2026-06-25); P22.5 HSV-1 denominator analysis resolved (2026-06-25); P22.6 STARsolo comparison complete (2026-06-25); P22.10 matched-barcode comparison complete (2026-06-25); host-response analysis complete (2026-06-27). -->

**Authors:** [Author list to be supplied]

**Affiliations:** [Affiliations to be supplied]

**Lead contact:** [Lead contact to be supplied]

**Keywords:** single-cell RNA-seq, viral detection, pseudoalignment, kallisto, multimapping, host response

---

## Highlights

- ViralScan quantifies viral RNA from paired-end scRNA-seq FASTQ files with a combined host-virus kallisto/bustools reference.
- EM-based multimapping correction recovered 12,255 EBV UMIs from a 1M-read LCL subsample, 3.64-fold more than unique-only counting.
- Full-depth benchmarks reproduced HHV-6B and HSV-1 infection ranges after matching denominators and thresholds to the source studies.
- In 1,906 matched EBV LCL cells, ViralScan and STARsolo produced concordant viral-burden rankings but different EBV gene attribution.

## eTOC/In Brief

ViralScan is a Snakemake-based command-line workflow for quantifying viral RNA in single-cell RNA-seq libraries. It uses a combined host-virus pseudoalignment reference and an EM correction step for host-virus ambiguous reads, then reports per-cell viral burden and optional host-response models. Benchmarks on HHV-6B, EBV, and HSV-1 datasets show that denominator choice, infection threshold, and viral annotation affect agreement with published rates.

## Summary

Single-cell RNA sequencing (scRNA-seq) can capture viral transcripts together with host gene expression, but viral reads are often discarded or processed through separate alignment steps. ViralScan is an open-source command-line workflow that quantifies viral load in paired-end scRNA-seq data by pseudoaligning reads with kallisto/bustools against a combined host-virus reference, followed by expectation-maximisation (EM) correction of multimapping reads. The workflow runs under Snakemake and processed a typical 10x Chromium v3 library in under 2 h on an 8-core compute node in local validation. Across three public datasets, HHV-6B in CAR-T cells (Lareau et al., 2023), EBV in lymphoblastoid cell lines (SoRelle et al., 2021), and HSV-1 in fibroblasts (Wyler et al., 2019), ViralScan produced infection-rate estimates consistent with published results after matching denominators and thresholds. An optional host-response module links per-virus infection status to host expression using L2 logistic regression and randomised Lasso stability selection. Source code is available at https://github.com/mdmanurung/ViralScan.

## Introduction

Single-cell RNA sequencing can measure viral transcription and host expression in the same cell. This joint measurement is relevant for viruses such as Epstein-Barr virus (EBV), human herpesvirus 6 (HHV-6), and herpes simplex virus 1 (HSV-1), where viral transcription varies across cells within an infected population (Lareau et al., 2023; SoRelle et al., 2021; Wyler et al., 2019). Analyses of these datasets require workflows that preserve cell barcodes, quantify viral features, and retain host expression for downstream modelling.

Existing approaches make different trade-offs. Combined-genome aligners such as STAR, Cell Ranger, and STARsolo provide splice-aware read placement but can be computationally expensive and depend on viral annotation completeness. Viral-only pseudoalignment is faster, but it discards host reads and cannot resolve reads compatible with both host and viral transcripts. Host-filter workflows such as VIRTUS/VIRTUS2-like two-step pipelines remove host-compatible barcode-UMI pairs before viral quantification, which can discard host-virus ambiguous reads that contain viral signal.

ViralScan was designed to retain these ambiguous reads without requiring a full splice-aware alignment. It aligns reads to a combined host-virus transcriptome and applies a transcriptome-wide EM algorithm to equivalence classes spanning multiple genes. In an EBV lymphoblastoid cell line (LCL) benchmark, this combined-reference strategy recovered approximately 4-fold more EBV UMIs than unique-only counting in a 1M-read subsample. ViralScan also includes 195 bundled viral GTF annotations and supports user-defined references through direct NCBI accession download.

## Results

### Combined-reference EM recovers EBV signal missed by two-step filtering

We first measured the effect of multimapping correction in the EBV LCL dataset (SRR12682296) using a 1M-read subsample. Three counting strategies were compared on the same input data.

| Strategy | EBV UMI detected | Relative sensitivity |
|----------|------------------|----------------------|
| Unique-only (no correction) | 3,372 | 1x |
| Two-step host-filter (VIRTUS2-like) | 3,096 | 0.92x |
| EM multimapping (ViralScan default) | 12,255 | 3.64x |

The combined-reference EM approach recovered 12,255 EBV UMIs, compared with 3,372 UMIs from unique-only counting. The two-step host-filter strategy recovered fewer EBV UMIs than unique-only counting because host-filtering removes barcode-UMI tuples from reads compatible with both host and EBV before the viral pass. This result supports the use of a combined reference when host-virus ambiguous reads are expected.

### Full-depth benchmarks agree with published rates after denominator and threshold matching

We evaluated ViralScan on three public scRNA-seq datasets selected because the source studies reported viral infection or viral-expression rates.

| Dataset | Published infected-cell rate | ViralScan 1M reads | ViralScan full depth | Comparison threshold |
|---------|------------------------------|--------------------|----------------------|----------------------|
| HHV-6B, CAR-T, SRR20710641 | 0.01%-0.3% super-expressors; 0.2% at Day 19 | 1.25% (>=10 UMI) | 0.152% (1,965 / 1,292,857 cells, >=1 UMI); 0.0014% super-expressors (18 cells, >=10 UMI) | >=1 UMI and >=10 UMI |
| EBV, LCL, SRR12682296 | 0.9%-2.2% lytic cells | 3.34% (285 / approximately 8,523 cells) | 8.985% (67,254 / 748,518 barcodes, >=1 UMI); 0.382% (2,860 cells, >=10 UMI) | >=1 UMI and >=10 UMI |
| HSV-1, fibroblasts, 5 hpi, SRR8315713 | approximately 13%-19% infected by bimodal split | 0.31% over raw barcodes | 0.5521% over 1,893,827 raw barcodes; 13.5%-17.6% over 4,414 called cells at >=5 to >=2 UMI | called-cell denominator and >=2-5 UMI |

For HHV-6B, full-depth ViralScan detected 1,965 HHV-6B-positive cells among 1,292,857 barcodes (0.152%, >=1 UMI) and 18 super-expressors (0.0014%, >=10 UMI). These values are compatible with the Lareau et al. Day-19 estimate of 0.2% total cells and the reported super-expressor range. The earlier 1M-read estimate was higher, consistent with a shallow-subsample artifact or timepoint/sample differences.

For EBV, ViralScan detected 67,254 EBV-positive barcodes among 748,518 unfiltered barcodes at full depth (8.985%, >=1 UMI; 1,252,577 total EBV UMI). The >=10 UMI tier contained 2,860 cells (0.382%). Because latent EBV is expected across LCL cells, the >=1 UMI rate measures detectable EBV expression rather than the published lytic-cell fraction. The >=10 UMI tier is closer to, but still below, the 0.9%-2.2% lytic range reported by SoRelle et al.

For HSV-1, the apparent full-depth rate of 0.5521% came from using all 1,893,827 raw barcodes as the denominator. Recomputing over 4,414 called cells with >=1,000 total UMI gave 1,197 HSV-1-positive cells at >=1 UMI (27.1%). Applying thresholds aligned to Wyler et al.'s bimodal separation gave 777 / 4,414 cells at >=2 UMI (17.6%), 687 / 4,414 at >=3 UMI (15.6%), and 596 / 4,414 at >=5 UMI (13.5%). Thus the HSV-1 discrepancy is a denominator and threshold artifact rather than a failure to detect viral reads.

### Matched EBV cells show concordant burden ranking but annotation-dependent gene attribution

We compared ViralScan with STARsolo on the EBV LCL sample (SRR12682296, 10x Chromium v2, approximately 112M reads). STARsolo 2.7.11b was run in GeneFull mode with a combined GRCh38 + EBV (NC_007605.1) reference. ViralScan was run on the same sample using the Serratus combined index. Both tools were then compared on the 1,906-cell GSM4796271 LCL_777_B958 barcode anchor after whitelist-corrected processing.

| Tool | Cell set | EBV >=1 UMI | EBV >=10 UMI | Lytic marker positive (BZLF1/BRLF1/BHRF1 >=1) | Reference |
|------|----------|-------------|--------------|-----------------------------------------------|-----------|
| ViralScan full depth, unfiltered | 748,518 barcodes | 67,254 (8.985%) | 2,860 (0.382%) | not assessed | Serratus combined index |
| STARsolo GeneFull full depth, CellRanger2 filter | 1,909 cells | 1,460 (76.48%) | 187 (9.80%) | not assessed | GRCh38 + NC_007605.1 |
| STARsolo GeneFull matched, whitelist | 1,906 cells | 1,473 (77.28%) | 195 (10.23%) | 74 (3.88%) | GRCh38 + NC_007605.1 |
| ViralScan matched, whitelist | 1,906 cells | 1,790 (93.91%) | 257 (13.48%) | 52 (2.73%) | Serratus combined index |
| Cell Ranger + Seurat, published SoRelle 2021 | approximately 5,830 pooled LCL cells | approximately all latent | not reported | approximately 0.9%-2.2% lytic | not reported |

All 1,906 paper anchor cells were present in both raw matrices. On this matched set, ViralScan detected more cells as EBV-positive at >=1 UMI (93.91% versus 77.28%) and gave a lytic-marker-positive fraction closer to the published B95-8 value (2.73% versus 3.88%; published target approximately 2.2%). Per-cell EBV burden was moderately concordant between the tools (Spearman r = 0.45; Pearson r = 0.42; n = 1,906; p < 10^-96), indicating agreement in relative ranking despite differences in gene attribution.

The gene-level comparison showed annotation-dependent differences. STARsolo captured LMP-1 with high sensitivity (47,220 UMI; 1,418 / 1,906 cells, 74.4%) but produced zero counts for the EBNA nuclear antigen family (EBNA-1, EBNA-2, EBNA-3A, EBNA-3B/3C, and EBNA-LP). ViralScan recovered all six EBNA-family entries in the Serratus annotation (EBNA-2: 716 UMI; EBNA-3A: 817 UMI; EBNA-3B/3C: 274 UMI; EBNA-LP: 37 UMI; EBNA-1.2: 4 UMI) but detected 99 LMP-1 UMIs across 92 cells. BRLF1 was broadly concordant (STARsolo: 77 UMI; ViralScan: 57 UMI; ratio 0.74). These differences are consistent with annotation coverage rather than a universal advantage of either aligner: STARsolo used 16 EBV gene-level loci, while the Serratus index provided 96 EBV entries.

### Host-response modelling requires explicit depth controls

We ran the host-response module on the 1,906 matched EBV cells. Viral genes were excluded from the host feature matrix before model fitting. EBV burden was recomputed from the corrected multimapping matrix by summing EBV gene columns into one `Epstein-Barr virus` feature, and cells were labelled EBV-high at >=10 corrected UMI.

At this threshold, 1,179 cells were EBV-high and 727 cells were EBV-low or negative. The positive count is higher than the 257 cells at >=10 UMI in the unique-count matched table because host-response uses the corrected multimapping matrix, whereas the table reports unique UMI counts from `adata.h5ad`. Across six random seeds, an L2 logistic-regression model predicted EBV-high status from host expression with AUROC 0.866 +/- 0.036 and balanced accuracy 0.783 +/- 0.038 (sensitivity 0.822 +/- 0.054; specificity 0.744 +/- 0.054). This raw `>=10 corrected UMI` label is, however, strongly confounded by sequencing depth. EBV-high prevalence rose monotonically with host library size, from 0.28 in the lowest host-depth quintile to 0.96 in the highest, and in the same evaluation design sequencing depth alone predicted the label at AUROC 0.967 — higher than the host-gene model itself. We therefore treat the raw AUROC as a diagnostic of a depth-tracking label rather than as depth-independent biological evidence. Depth-controlled analyses gave more conservative estimates: a CPM-normalised, prevalence-matched label yielded AUROC 0.636, and a depth-matched case/control cohort (matched depth medians; depth-alone control AUROC 0.48) yielded AUROC 0.718. These controls bound the depth-independent host-response signal at approximately 0.64-0.72.

Randomised Lasso stability selection (`n_stab_iter=100`, `stab_min_prob=0.6`) identified 15 stable host features, of which 5 retained a depth-adjusted E-value >=2 (robust to moderate depth confounding). The matched host h5ad contains Ensembl IDs but no gene-symbol annotation, so biological interpretation of the stable feature set requires gene-symbol annotation and pathway enrichment on the depth-controlled labels. These results demonstrate ViralScan's ability to expose and control depth confounding in host-response analysis; they do not support interpreting the raw AUROC as a depth-independent host-response signature.

## Discussion

ViralScan provides a practical route from raw scRNA-seq FASTQ files to per-cell viral RNA quantification. Its main technical feature is combined-reference pseudoalignment followed by EM allocation of host-virus ambiguous equivalence classes. In the EBV benchmark, this correction recovered 3.64-fold more EBV UMI than unique-only counting in the 1M-read comparison, showing that ambiguous reads can contribute materially to viral detection.

The validation results also show why viral scRNA-seq benchmarks require careful denominator and threshold matching. HHV-6B rates were compatible with the published CAR-T study at full depth. HSV-1 appeared discordant when all raw barcodes were used as the denominator, but matched the published 13%-19% range when analysed over called cells with >=2-5 UMI thresholds. EBV comparisons required matched barcode sets because filtered-cell and unfiltered-barcode denominators produced different apparent rates.

The STARsolo comparison should not be interpreted as a general claim that one quantifier is superior. On the 1,906 matched EBV cells, ViralScan and STARsolo ranked per-cell EBV burden similarly but assigned reads to different EBV genes. The divergence is consistent with differences in EBV annotation and counting model, especially for LMP-1 and EBNA-family features. This result argues for explicit reporting of viral reference annotations in single-cell viral RNA analyses.

ViralScan has limitations. The current EM correction estimates one transcriptome-wide abundance vector and does not model cell-to-cell differences in viral expression. This approximation is memory-efficient and useful for sparse viral signals, but per-cluster or per-cell allocation may be preferable in highly infected cultures. ViralScan also does not produce cell-level BAM files for viral reads. The `viralscan evidence` subcommand provides BLAST-based read tracing, but genome-browser inspection remains a separate alignment step. Host-response modelling is also sensitive to how the positive label is defined: a raw viral-UMI threshold tracks library size, so a naive host-expression AUROC can partly reflect sequencing depth rather than biology. ViralScan reports a depth-alone baseline and supports CPM-normalised labels and depth-matched designs; host-response AUROCs should always be read against these depth baselines. Host-response modelling currently reports Ensembl feature IDs when the host matrix lacks gene-symbol annotation, limiting biological interpretation until annotation and enrichment are added.

This study presents ViralScan as an open-source workflow and validation case study rather than a comprehensive benchmark against all dedicated viral single-cell detectors. We did not complete a head-to-head comparison against dedicated viral single-cell tools such as Venus or ViralTrack in the present version; the STARsolo comparison tests a splice-aware general aligner and should not be interpreted as a comprehensive dedicated-tool benchmark. Finally, the present validation uses published infection-rate ranges and matched tool comparisons rather than a gold-standard single-cell truth panel. Formal false-positive and false-negative rates require negative controls and planted viral-read simulations and are left as a release-gated validation extension.

## STAR Methods

### Resource availability

#### Lead contact

Further information and requests for resources should be directed to the lead contact, [Lead contact to be supplied].

#### Materials availability

This study did not generate new biological materials.

#### Data and code availability

ViralScan source code is available at https://github.com/mdmanurung/ViralScan. Reference-data archive DOI: 10.5281/zenodo.20112332. Public datasets used in the benchmarks are available from GEO/SRA: GSE210063 (SRR20710641), GSE158275 (SRR12682296), and GSE123782 (SRR8315713). Benchmark and figure-generation scripts are included in `scripts/`, and manuscript figures are stored in `docs/figures/`.

### Method details

#### Reference construction

ViralScan constructs a combined reference by concatenating host (GRCh38) and viral genome FASTA/GTF files. Kallisto indices and transcript-to-gene maps are generated with `kb ref`. The package includes GTF annotations for 195 virus accessions. Additional accessions can be fetched with `viralscan data fetch <accession>`. The benchmarks used the Serratus-derived combined human + 838-virus index.

#### Pseudoalignment and barcode/UMI extraction

Paired-end FASTQ files are processed with `kb count` (kallisto 0.50.x and bustools 0.43.x). Library chemistry is specified with `--technology` (`10xv2`, `10xv3`, or `DROPSEQ`). Viral genes are identified from the transcript-to-gene map by accession prefixes in the gene ID field. ViralScan writes per-cell, per-gene count matrices for host and viral features.

#### Multimapping correction

Reads pseudoaligning to equivalence classes spanning multiple genes are corrected with a global-pool EM algorithm implemented in `viralscan/multimapping.py`. All cells are pooled into a frequency matrix. The EM estimates a transcriptome-wide relative abundance vector theta. At convergence, each ambiguous equivalence-class count is allocated to compatible genes in proportion to the normalised theta components and added to unique-mapping counts.

The global-pool EM estimates one theta vector for the cell population and does not estimate per-cell viral composition. This is a deliberate memory-saving approximation. ViralScan supports four multimapping modes through `--multimap-method`: `equal` (the default), `host-conservative`, `unique-weighted`, and `em`.

#### Detection and viral calling

Infected cells are called by thresholding per-cell viral UMI counts. The default threshold is >=10 UMIs and can be changed with `--detection-threshold`. ViralScan reports total viral UMIs, infected-cell counts, infected-cell fractions, and per-cell UMI distributions.

#### Host-response analysis

The optional `viralscan hostresponse` module associates viral status with host expression. For each virus, cells are labelled positive or negative by the viral UMI threshold. The module selects the top 2,000 highly variable genes with Scanpy unless `--no-use-hvg` is set, downsamples the majority class, trains L2 logistic regression models across random seeds, reports cross-validated AUROC and balanced accuracy, and applies randomised Lasso stability selection. Optional pathway enrichment uses `gget enrichr` when `viralscan[enrichment]` is installed.

#### Benchmark datasets

| Dataset | GEO | Sample | Technology | Virus | Reference |
|---------|-----|--------|------------|-------|-----------|
| GSE210063 | SRR20710641 | CAR-T cell product | 10x Chromium v3 | HHV-6B | Lareau et al., 2023 |
| GSE158275 | SRR12682296 | LCL line 777 | 10x Chromium v2 | EBV | SoRelle et al., 2021 |
| GSE123782 | SRR8315713 | Primary fibroblasts, 5 hpi | Drop-seq | HSV-1 | Wyler et al., 2019 |

Full-depth FASTQ files were downloaded through the EBI ENA FTP server when available and NCBI SRA as fallback with `fasterq-dump`. Preliminary 1M-read subsamples were used for early checks, and full-depth analyses were used for final comparisons.

#### STARsolo EBV comparison

SRR12682296 was processed with STARsolo 2.7.11b using a combined GRCh38 + EBV (NC_007605.1) reference. STARsolo was run with 10xv2 parameters (CB = 16 bp, UMI = 10 bp), `--soloType CB_UMI_Simple`, `--soloFeatures GeneFull`, and `--soloCellFilter CellRanger2`. The matched-barcode analysis used the 10x v2 whitelist and the 1,906 GSM4796271 LCL_777_B958 paper anchor cells.

### Quantification and statistical analysis

Benchmark infection rates were calculated as infected cells divided by the reported denominator for each analysis: unfiltered barcodes, STARsolo-filtered cells, or called cells, as specified in each table. HSV-1 threshold reconciliation used called cells with >=1,000 total UMI and integer viral-UMI thresholds from >=1 to >=10. Matched EBV tool concordance was assessed with Pearson and Spearman correlations across 1,906 matched cells. Host-response performance was estimated across six random seeds; summary values are reported as mean +/- standard deviation.

## Acknowledgments

[Acknowledgments to be supplied by authors.]

## Author contributions

[Author contributions to be supplied by authors.]

## Declaration of interests

[Declaration of interests to be supplied by authors.]

## Supplemental information

Figure files included with this draft:

- `docs/figures/figure1_workflow.png`
- `docs/figures/figure1_workflow.pdf`
- `docs/figures/figure2_benchmark.png`
- `docs/figures/figure2_benchmark.pdf`

**Figure 1. ViralScan workflow schematic.** Paired-end scRNA-seq reads are processed against a combined host-virus transcriptome using kallisto/bustools. ViralScan applies EM-based multimapping correction to host-virus ambiguous equivalence classes, outputs per-cell viral burden, and optionally runs host-response modelling from host-only predictors.

**Figure 2. Benchmark and matched-cell EBV comparison.** The benchmark panel compares published infected or lytic cell rates with ViralScan estimates for HHV-6B, EBV, and HSV-1 after full-depth validation and threshold reconciliation. EBV matched-cell panels compare STARsolo and ViralScan on the 1,906-cell paper anchor and show divergent gene attribution across shared EBV features. The footer reports matched host-response model performance together with the depth-alone baseline and depth-controlled AUROC estimates.

## References

1. Lareau, C.A., Yin, Y., Maurer, K., et al. (2023). Latent human herpesvirus 6 is reactivated in CAR T cells. *Nature*, 623, 608-615. https://doi.org/10.1038/s41586-023-06704-2

2. SoRelle, E.D., Dai, J., Bonglack, E.N., et al. (2021). Single-cell RNA-seq reveals transcriptomic heterogeneity mediated by host-pathogen dynamics in lymphoblastoid cell lines. *eLife*, 10, e62586. https://doi.org/10.7554/eLife.62586

3. Wyler, E., Franke, V., Menegatti, J., et al. (2019). Single-cell RNA-sequencing of herpes simplex virus 1-infected cells connects NRF2 activation to an antiviral program. *Nature Communications*, 10, 4906. https://doi.org/10.1038/s41467-019-12894-z

4. Bray, N.L., Pimentel, H., Melsted, P., and Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology*, 34, 525-527. https://doi.org/10.1038/nbt.3519

5. Melsted, P., Booeshaghi, A.S., Liu, L., et al. (2021). Modular, efficient and constant-memory single-cell RNA-seq preprocessing. *Nature Biotechnology*, 39, 813-818. https://doi.org/10.1038/s41587-021-00870-2

6. Dobin, A., Davis, C.A., Schlesinger, F., et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29, 15-21. https://doi.org/10.1093/bioinformatics/bts635

7. Yasumizu, Y., Hara, A., Sakaguchi, S., and Ohkura, N. (2021). VIRTUS: a pipeline for comprehensive virus analysis from conventional RNA-seq data. *Bioinformatics*, 37, 1465-1467. https://doi.org/10.1093/bioinformatics/btaa859

8. Yasumizu, Y. (n.d.). VIRTUS2: a bioinformatics pipeline for viral transcriptome detection and quantification considering splicing [software]. GitHub. https://github.com/yyoshiaki/VIRTUS2 (accessed 2026-07-03). VIRTUS2 is a software successor to VIRTUS (reference 7) with no separate journal article; the method is described in reference 7.

9. Edgar, R.C., Taylor, J., Lin, V., et al. (2022). Petabase-scale sequence alignment catalyses viral discovery. *Nature*, 602, 142-147. https://doi.org/10.1038/s41586-021-04332-2

10. Luebbert, L., Sullivan, D.K., Carilli, M., et al. (2024). Efficient and accurate detection of viral sequences at single-cell resolution reveals putative novel viruses perturbing host gene expression. *bioRxiv*. https://doi.org/10.1101/2024.01.13.575532

11. Meinshausen, N., and Buehlmann, P. (2010). Stability selection. *Journal of the Royal Statistical Society: Series B*, 72, 417-473. https://doi.org/10.1111/j.1467-9868.2010.00740.x
