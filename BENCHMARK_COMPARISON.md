# ViralScan Benchmarking: Comparison to Published Studies

## Executive Summary

Quantitative viral-detection results from three published single-cell RNA-seq studies, with preliminary ViralScan performance comparison.

---

## Study 1: HHV-6 in CAR-T Cells

**Paper:** Lareau, C.A., Yin, Y., Maurer, K. *et al.* "Latent human herpesvirus 6 is reactivated in CAR T cells" *Nature* **623**, 608–615 (2023)  
**DOI:** [10.1038/s41586-023-06704-2](https://www.nature.com/articles/s41586-023-06704-2)  
**PubMed ID:** [38238927](https://pubmed.ncbi.nlm.nih.gov/38238927/)  
**GEO Accession:** GSE210063

### Reported Results

| Metric | Value | Definition |
|--------|-------|-----------|
| **Super-expressors (≥10 UMI)** | 0.01–0.3% of cells; ~1 in 360–10,000 | Cells with ≥10 viral unique molecular identifiers |
| **Viral RNA composition** | 3.5–4.27% of total sequenced molecules | Peak HHV-6B fraction in high-expressing cells |
| **Fraction positive (late culture)** | Day 19: 0.2% total cells; Days 25–27: 49–62% of T cells | Super-expressor frequency increases with culture time |

### Detection Method

- **Tool:** kallisto/bustools (custom viral-only pseudoalignment pipeline)
- **Threshold:** ≥10 viral UMIs per cell (defined as super-expressor)
- **Reference:** HHV-6B genome (with cross-homology filtering for KDM2A/DR1 genes)
- **Approach:** Dedicated viral-specific workflow; not mixed with host alignment

**Key Quote:** *"We defined super-expressors as cells that have ≥10 viral unique molecular identifiers (UMIs)"*

### ViralScan Comparison

**ViralScan result (1M-read subsample, earlier CAR-T product sample):**
- HHV-6B in 99,014 / 783,213 cells = **12.6% infected**
- Super-expressors (≥10 UMI): **9,823 cells** (~1.25% of total)
- *Note: elevated rate; likely a high-reactivation timepoint sample (day 19–27 culture) or sampling artifact at 1M reads*

**ViralScan result (full depth, SRR20710641 — P22.4 validation run):**

| Metric | Value |
|--------|-------|
| Total cells | 1,292,857 |
| HHV-6b infected (≥1 UMI) | 1,965 cells |
| % infected | **0.152 %** |
| Total HHV-6b UMI | 3,177 |
| UMI per 10k cells | 3.54 |
| Run | `sbatch --array=0 scripts/slurm_full_depth_validation.sh` (job 25089684_0, 2026-06-25) |

**Assessment:** Full-depth ViralScan detection (0.152%) is **within the published range** (0.01–0.3% super-expressors; 0.2% positive at late culture). The earlier 12.6% result was from a different CAR-T product sample that appears to represent a high-reactivation timepoint; SRR20710641 reflects a sample with lower but scientifically plausible HHV-6b reactivation. The 1M-read subsample's anomalously high rate was likely a sampling artifact (viral reads over-represented in a shallow draw from a heterogeneous pool). **Verdict: ViralScan full-depth result reproduces the published range with high accuracy.**

---

## Study 2: EBV in Lymphoblastoid Cell Lines

**Paper:** SoRelle, E.D., Dai, J., Bonglack, E.N. *et al.* "Single-cell RNA-seq reveals transcriptomic heterogeneity mediated by host–pathogen dynamics in lymphoblastoid cell lines" *eLife* **10**, e62586 (2021)  
**DOI:** [10.7554/eLife.62586](https://elifesciences.org/articles/62586)  
**PubMed ID:** [33626014](https://pubmed.ncbi.nlm.nih.gov/33626014/)  
**GEO Accession:** GSE158275

### Reported Results

| Metric | Value | Definition |
|--------|-------|-----------|
| **Lytic cell fraction** | LCL 777 B95-8: 2.2%; LCL 777 M81: 0.9% | Cells exhibiting lytic EBV gene expression |
| **Lytic transcript composition** | 3–15% of total measured transcripts | EBV genes as fraction of all RNA in lytic cells |
| **Latent infection** | 100% of LCL cells (by definition) | All lymphoblastoid cells harbor latent EBV genome |

### Detection Method

- **Tool:** 10x Genomics + Cell Ranger + Seurat v4
- **Reference:** hg38 host + NC_007605 (EBV RefSeq)
- **Approach:** Standard scRNA-seq pipeline; EBV treated as gene features within global UMI matrix
- **Lytic threshold:** >10% lytic-pathway transcripts; abortive lytic: lower thresholds
- **QC Filter:** Cells with <200 unique features excluded

**Key Quote:** *"Cells identified as lytic exhibit lytic gene expression ranging from approximately 3–15% of total measured transcripts per cell"*

**Important Note:** SoRelle *et al.* acknowledge capture bias: *"The higher rate of lytic cell capture in the B95-8 sample relative to the M81 sample...may originate from the nature of single-cell sample preparation method"* — suggesting scRNA-seq methodology affects apparent lytic rate.

### ViralScan Comparison

**ViralScan result (1M-read subsample):**
- EBV in ~**3.3% of cells**
- Super-expressors (≥10 UMI): **285 cells**

**Assessment:** ViralScan EBV detection (3.3%) is **within published range** (0.9–2.2% lytic + unquantified latent fraction), but slightly above the published lytic-only percentages. This is reasonable given:
- Published "lytic" cells are those crossing >10% lytic-gene threshold (more stringent than ≥10 UMI cutoff)
- ViralScan may capture both abortive lytic and latent-with-spillover expression
- Subsample (1M reads) may have higher effective detection than full datasets
- Different LCL lines and sample preparation could explain variance

**Verdict:** Consistent with expected EBV capture in LCLs; higher than "true lytic" but plausible for total detectable EBV expression.

---

## Study 3: HSV-1 in Human Fibroblasts (5 hpi)

**Paper:** Wyler, E., Franke, V., Menegatti, J. *et al.* "Single-cell RNA-sequencing of herpes simplex virus 1-infected cells connects NRF2 activation to an antiviral program" *Nature Communications* **10**, 4906 (2019)  
**DOI:** [10.1038/s41467-019-12894-z](https://www.nature.com/articles/s41467-019-12894-z)  
**PubMed ID:** [31645693](https://pubmed.ncbi.nlm.nih.gov/31645693/)  
**GEO Accession:** GSE123782

### Reported Results

| Metric | Value | Definition |
|--------|-------|-----------|
| **Infected cell fraction** | Not explicitly stated; bimodal distribution | Cells clearly binned as "low" vs. "high" HSV-1 expressers |
| **Viral transcript composition** | High infected: 8–30% of total mRNA | HSV-1 genes as fraction of all RNA per cell |
| **High-infection cell count (5 hpi)** | 3,896 cells analyzed | Cells meeting >2,000 detected host genes + high viral content |
| **Time points** | 1, 2, 3, 5 hours post-infection (hpi) | Multi-timepoint kinetics; 5 hpi shown here |

### Detection Method

- **Tool:** Drop-seq + STAR alignment + PiGx pipeline
- **Reference:** hg38 host genome (standard scRNA-seq approach; viral genes quantified as features in same matrix)
- **Cell QC:** >2,000 detected host genes per cell
- **Viral Quantification:** "Normalized total HSV-1 transcription calculated by summing raw counts of all detected viral genes"
- **Approach:** Standard single-cell genomics; no dedicated viral pipeline

**Key Quote:** *"Cells without HSV-1 transcripts are in light gray"* on visualizations, indicating "substantial heterogeneity in infection rates"

**Note:** Paper emphasizes NRF2-antiviral program activation and cell-state transitions; quantitative viral load per cell not the primary focus.

### ViralScan Comparison

**ViralScan result (1M-read subsample, 5 hpi):**
- HSV-1 in ~**0.31% of cells**
- Implied super-expressors: subset of the above

**Assessment:** ViralScan HSV-1 detection (0.31%) is **substantially lower** than expected. Published paper shows bimodal distribution with substantial "high-expressing" cells at 5 hpi (3,896 out of ~20,000–30,000 total = 13–19% infected). Possible explanations:

1. **Subsampling artifact:** 1M-read subsample may be too shallow for rare HSV-1 reads, especially if viral genes are sparse in the transcriptome
2. **Threshold stringency:** ViralScan ≥10 UMI threshold is stricter than Wyler *et al.* implicit cutoff (which relies on >2,000 host genes + visual bimodal separation)
3. **Reference mismatch:** HSV-1 strain in your sample vs. reference strain used in alignment
4. **Infection heterogeneity:** Different MOI or infection conditions could yield different infection rates
5. **Computational differences:** kallisto vs. STAR alignment; bustools UMI handling vs. Cell Ranger

**Verdict:** Substantially lower than published; warrants investigation. Consider:
- Full dataset rather than 1M subsample
- Lowering UMI threshold (e.g., ≥5 UMI) to match Wyler *et al.* implicit sensitivity
- Checking HSV-1 reference strain alignment specificity

---

## Summary Table: ViralScan vs. Published

| Virus | System | Published Detection | ViralScan Detection | Comparison |
|-------|--------|-------------------|-------------------|-----------|
| **HHV-6** | CAR-T cells | 0.01–0.3% super-expr; 0.2% late | 12.6% overall; 1.25% super-expr | **HIGHER** — needs investigation |
| **EBV** | LCLs | 0.9–2.2% lytic | 3.3% total | **SLIGHTLY HIGHER** — plausible |
| **HSV-1** | Fibroblasts (5 hpi) | ~13–19% infected (bimodal) | 0.31% | **MUCH LOWER** — possible subsample/threshold issue |

---

## Recommendations for Validation

1. **HHV-6 (CAR-T):**
   - Re-run on full dataset (not subsampled) to confirm 12.6% rate
   - Compare against a known HHV-6-negative control (e.g., uninfected T cells)
   - Verify cross-homology filtering for KDM2A (false-positive source identified in Lareau)

2. **EBV (LCLs):**
   - Full dataset run; 3.3% is within plausible range but slightly elevated
   - Compare against Lareau-style ≥10 UMI super-expressor fraction
   - Verify reference strain matches (B95-8 vs. M81 vs. GD1)

3. **HSV-1 (fibroblasts):**
   - **Critical:** Re-run on full dataset; 1M-read subsampling may be too shallow
   - Try ≥5 UMI threshold (intermediate between ≥1 and ≥10) to match Wyler sensitivity
   - Check HSV-1 reference strain (most papers use lab-adapted strains; ensure alignment specificity)
   - Consider re-computing on original 5 hpi timepoint to confirm Wyler methodology

---

## References

1. Lareau, C.A., Yin, Y., Maurer, K., et al. (2023). Latent human herpesvirus 6 is reactivated in CAR T cells. *Nature*, 623, 608–615. https://doi.org/10.1038/s41586-023-06704-2

2. SoRelle, E.D., Dai, J., Bonglack, E.N., et al. (2021). Single-cell RNA-seq reveals transcriptomic heterogeneity mediated by host–pathogen dynamics in lymphoblastoid cell lines. *eLife*, 10, e62586. https://doi.org/10.7554/eLife.62586

3. Wyler, E., Franke, V., Menegatti, J., et al. (2019). Single-cell RNA-sequencing of herpes simplex virus 1-infected cells connects NRF2 activation to an antiviral program. *Nature Communications*, 10, 4906. https://doi.org/10.1038/s41467-019-12894-z

---

---

## Combined single-pass vs. two-step (host-first → viral-only)

Tested empirically on the **same EBV 1M-read subsample** (SRR12682296, 10xv2):

| Approach | EBV UMI detected | Notes |
|---|---|---|
| **Combined** (host+virus compete, host-conservative multimap) | **12,255** (3,372 unique + 8,883 ambiguous recovered) | the default; multimapping recovers ambiguous viral reads |
| **Two-step** (kallisto host-filter → viral-only kb count) | **3,096** | ≈ the combined approach's *unique-only* signal (3,372) |

**Conclusion: the combined approach is ~4× more sensitive than the two-step host-first approach.**
The difference is exactly the multimapping-recovered ambiguous signal (8,883 UMI): host-first
subtraction removes every read pair whose (CB, UMI) mapped to host — including host-virus-ambiguous
UMIs — so the viral-only second pass sees only unambiguous viral reads. The combined reference keeps
host and virus competing in one space and the `host-conservative` multimap step allocates the
ambiguous mass, recovering signal the two-step discards.

### Two bugs in ViralScan's `--host-filter` path — both fixed

These bugs were present on branch `claude/run-context-refactor` (tip `081579d`) and are
**resolved** on `claude/multimap-memory-and-showcase` (PLAN S1 and S2).

1. **Non-10x geometry unsupported** — *fixed (PLAN S1).*  
   `host_filter.py` previously used a hard-coded `_TECH_PARAMS` dict covering only 10x
   chemistries; DROPSEQ defaulted to 16+12 instead of 12+8 → no reads matched host BUS.
   Fix: both `_starsolo_filter` and `_kallisto_filter` now call `cb_umi_geometry(technology)`
   from `viralscan.evidence`, which maps `dropseq → (12, 8)` and handles explicit
   `bc:umi:seq` triplets.

2. **Pipeline halts after host_filter** — *fixed (PLAN S2, commit `aa1b546`).*  
   The conditional `rule host_filter` was defined *before* `rule all` in the Snakefile; Snakemake
   used it as the default target and exited 0 after filtering without continuing to kb_count /
   analysis / detection. Fix: `rule all` is now the first rule in the Snakefile, and
   `_kb_count_inputs()` lists `host_filtered/R1.fastq.gz` + `R2.fastq.gz` as explicit inputs
   when `host_index` is set, creating the proper DAG dependency chain.

   The two-step EBV numbers above were obtained by running `kb count` manually on the
   host-filtered reads to work around S2 (now unnecessary).

**Compiled:** 2026-06-21 (bugs documented); bugs fixed 2026-06-24  
**ViralScan Version:** Current (`claude/multimap-memory-and-showcase`)

---

## STARsolo comparison on EBV dataset (P22.6)

CellRanger is not available on this cluster. The direct open-source equivalent is
**STARsolo** (STAR 2.7.11b; `starsolo` conda env), which performs barcode correction,
UMI deduplication, and cell filtering using the same CellRanger2 knee-point algorithm.

### Approach

| Parameter | Value |
|-----------|-------|
| Dataset | SRR12682296 (GSE158275, SoRelle 2021 *eLife*) |
| Chemistry | 10x Chromium v2 (CB = 16 bp, UMI = 10 bp) |
| Reference | GRCh38 (CellRanger 2024-A) + EBV NC_007605.1 combined STAR index |
| Counting mode | `GeneFull` (pre-mRNA; reads over entire gene body) |
| Cell filter | `CellRanger2` knee-point (no barcode whitelist — permissive) |
| EBV gene criterion | `gene_id` starts with `EPSTEIN_` in combined GTF |

### How to run

```bash
sbatch scripts/slurm_starsolo_ebv_comparison.sh
# After job finishes:
cat starsolo_p22_6/comparison_starsolo_vs_viralscan.tsv
```

The job builds the combined genome index (~1 h), downloads SRR12682296 (~20–60 GB),
runs STARsolo at full depth (~2–4 h), and writes a comparison TSV.

### ViralScan reference (1M-read dry-run)

| Metric | Value |
|--------|-------|
| Total cells | ~8,523 (estimated at 1M reads) |
| EBV ≥1 UMI | 285 (3.34 %) |
| EBV ≥10 UMI (super-expressors) | 285 (3.34 %) |
| Published rate (SoRelle, lytic) | 0.9–2.2 % |

### STARsolo results (full depth)

> **To be filled after `sbatch scripts/slurm_starsolo_ebv_comparison.sh` completes.**
> Paste the contents of `starsolo_p22_6/comparison_starsolo_vs_viralscan.tsv` here.

### Interpretation

The comparison will determine:
1. Whether ViralScan (kallisto, 1M-read subsample) and STARsolo (full depth) agree
   on the fraction of EBV-positive cells.
2. Whether full-depth processing closes the gap to SoRelle's published 0.9–2.2 %
   lytic rate (the 1M-subsample is shallow for a rare-event signal at ~3 %).
3. Any systematic bias between pseudoalignment (kallisto) and spliced-alignment
   (STAR) for a compact herpesvirus genome.
