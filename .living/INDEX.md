<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-06

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 2 sections | 2026-07-06 | scRNA-seq associations: control depth AND %mito, and define labels depth-independently, Cross-validation: fit feature selection inside the split |
| decisions.md | 0 entries (large — read selectively) | 2026-07-06 | — |
| last-session.md | 15 entries | 2026-07-06 | 2026-07-06 — SH2.3: detection-level sibling cross-mapping warning (commit `f6786b2`), 2026-07-06 — Task 1: close-out depth-robust hostresponse module (commit `a1ad4a3`), 2026-07-06 (continued) — Tasks 2/3: cluster jobs submitted + manuscript covid section written, Task 2A/2B — STARsolo re-run + TTV read-origin (submitted), Task 2D — Manuscript covid specificity section (commit `76b769c`) |
| learnings.md | 15 entries (large — read selectively) | 2026-07-06 | cDNA-only host reference causes false-positive viral signal from GRCh38 non-coding reads, samtools view exits 1 on duplicate BAM header entry (NC_002076.2), covid_viralscan/results/ is gitignored — SURVEY_SUMMARY.md not tracked, summarize_survey.py --cellranger-outs skipped: script expects one barcode set for all samples, 64 GB references/ was untracked but NOT gitignored |
| log/ | 28 sessions | 2026-07-06 | viralscan (28) |
| findings/ | 4 findings across 5 topics | 2026-07-06 | unsupervised-feature-selection-leakage-is-often-negligible, reference_strategy_2x2, raw-count-thresholds-confound-with-sequencing-depth, host-transcriptome-encodes-viral-infection-state, covid-viralscan-no-sars2-anellovirus-dominant |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-06 (heuristic)

## Tag clusters

- **scrna-seq** (5 entries) — L-6, L-7, L-8, L-9, L-14
- **bioinformatics** (3 entries) — L-5, L-11, L-12
- **confounding** (3 entries) — L-6, L-7, L-8
- **gotcha** (3 entries) — L-11, L-13, L-14
- **verify-by-artifact** (3 entries) — L-11, L-12, L-14
- **barcodes** (2 entries) — L-4, L-14

## Most recent (10)

- [2026-07-06] L-1: cDNA-only host reference causes false-positive viral signal from GRCh38 non-coding reads
- [2026-07-06] L-2: samtools view exits 1 on duplicate BAM header entry (NC_002076.2)
- [2026-07-06] L-3: covid_viralscan/results/ is gitignored — SURVEY_SUMMARY.md not tracked
- [2026-07-06] L-4: summarize_survey.py --cellranger-outs skipped: script expects one barcode set for all samples
- [2026-07-06] L-15: Synthetic depth-proxy test: "fragile" is too strict; use "not robust"
- [2026-07-02] L-6: Controlling for a composite covariate that CONTAINS the tested feature is circular
- [2026-07-02] L-7: Mitochondrial genes need a %mito control before claiming them as biology
- [2026-07-02] L-11: kb ref / kallisto index have two silent FASTA-vs-GTF contracts
- [2026-07-02] L-12: bustools correct silently needs an UNCOMPRESSED whitelist
- [2026-07-02] L-13: Two panel dirs — bulk scan defaulted to the incomplete one

## By tag

- `scrna-seq`: L-6, L-7, L-8, L-9, L-14
- `bioinformatics`: L-5, L-11, L-12
- `confounding`: L-6, L-7, L-8
- `gotcha`: L-11, L-13, L-14
- `verify-by-artifact`: L-11, L-12, L-14
- `barcodes`: L-4, L-14
- `bustools`: L-12, L-14
- `causal-inference`: L-7, L-8
- `covid_viralscan`: L-3, L-4
- `git`: L-3, L-5
- `gitignore`: L-3, L-5
- `kallisto`: L-11, L-12
- `kb-python`: L-11, L-12
- `review`: L-6, L-9
- `summarize_survey`: L-3, L-4
- `whitelist`: L-12, L-14
- `10x`: L-14
- `5-prime`: L-14
- `NC_002076.2`: L-2
- `anellovirus`: L-1
- `bam-header`: L-2
- `bug`: L-10
- `bulk`: L-13
- `bulk-rnaseq`: L-1
- `cDNA-reference`: L-1
- `cellranger`: L-4
- `circularity`: L-6
- `classifier`: L-9
- `conventions`: L-10
- `count-data`: L-8
- `covariate`: L-6
- `cross-validation`: L-9
- `depth-confound`: L-15
- `duplicate-contig`: L-2
- `empty-droplets`: L-14
- `evalue`: L-15
- `false-positive`: L-1
- `feature-selection`: L-9
- `gem-x`: L-14
- `gene-selection`: L-7
- `gse128078`: L-13
- `gzip`: L-12
- `host-homology`: L-1
- `hostresponse`: L-15
- `hvg`: L-9
- `ingest`: L-5
- `kb-ref`: L-11
- `label-definition`: L-8
- `large-files`: L-5
- `leakage`: L-9
- `mito`: L-6
- `mitochondrial`: L-7
- `mycelium`: L-10
- `ngs_tools`: L-11
- `paths`: L-13
- `qc`: L-7
- `reference-build`: L-11
- `reference-panel`: L-13
- `references`: L-5
- `results`: L-3
- `samtools`: L-2
- `sequencing-depth`: L-8
- `set-e`: L-2
- `silent-correctness`: L-13
- `slurm-exit-code`: L-2
- `snakemake`: L-12
- `specificity`: L-1
- `star`: L-1
- `starsolo`: L-4
- `statistics`: L-6
- `synthetic-data`: L-15
- `testing`: L-15
- `thresholding`: L-8
- `tooling`: L-10
- `yaml`: L-10

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
