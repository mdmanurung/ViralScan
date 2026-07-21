<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-20

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 2 sections | 2026-07-06 | scRNA-seq associations: control depth AND %mito, and define labels depth-independently, Cross-validation: fit feature selection inside the split |
| decisions.md | 0 entries (large — read selectively) | 2026-07-17 | — |
| last-session.md | 24 entries | 2026-07-15 | 2026-07-15 (post-compaction) — T5, T6, SH2.6 closed; EVE job 25237061 running, Pending (as of session end), 2026-07-15 — EVE (Endogenous Viral Element) characterisation analysis, 2026-07-15 — aifi-scrna-pipeline skill pack installed, Pending (as of session end) |
| learnings.md | 22 entries (large — read selectively) | 2026-07-17 | cDNA-only host reference causes false-positive viral signal from GRCh38 non-coding reads, samtools view exits 1 on duplicate BAM header entry (NC_002076.2), covid_viralscan/results/ is gitignored — SURVEY_SUMMARY.md not tracked, summarize_survey.py --cellranger-outs skipped: script expects one barcode set for all samples, 64 GB references/ was untracked but NOT gitignored |
| log/ | 45 sessions | 2026-07-20 | viralscan (45) |
| findings/ | 4 findings across 5 topics | 2026-07-17 | covid-viralscan-no-sars2-anellovirus-dominant, unsupervised-feature-selection-leakage-is-often-negligible, reference_strategy_2x2, raw-count-thresholds-confound-with-sequencing-depth, host-transcriptome-encodes-viral-infection-state |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-20 (heuristic)

## Tag clusters

- **scrna-seq** (5 entries) — L-6, L-7, L-8, L-9, L-14
- **anellovirus** (4 entries) — L-1, L-18, L-20, L-21
- **covid** (4 entries) — L-16, L-18, L-20, L-21
- **verify-by-artifact** (4 entries) — L-11, L-12, L-14, L-16
- **barcodes** (3 entries) — L-4, L-14, L-16
- **bioinformatics** (3 entries) — L-5, L-11, L-12

## Most recent (10)

- [2026-07-21] L-25: Committed profiling artifacts go stale — re-profile current code before trusting them (ViralScan's EM hotspot was already fixed; the checked-in cProfile pointed at a dead path)
- [2026-07-20] L-24: Reorganize a sprawling repo non-destructively with a gitignored symlink view (relative links, tracked manifest+builder, disposable view)
- [2026-07-20] L-23: Tutorials must be verified against BOTH executed code AND the live CLI parser (RunConfig field names ≠ argparse flags; [skip-ci] docs rot silently)
- [2026-07-17] L-21: [2026-07-17]
- [2026-07-17] L-22: [2026-07-17]
- [2026-07-15] L-17: [2026-07-15]
- [2026-07-15] L-18: [2026-07-15]
- [2026-07-15] L-19: [2026-07-15]
- [2026-07-15] L-20: [2026-07-15]
- [2026-07-07] L-16: A ViralScan quant "redo" silently reuses cached outputs; and how to re-run only the detection tail
- [2026-07-06] L-1: cDNA-only host reference causes false-positive viral signal from GRCh38 non-coding reads
- [2026-07-06] L-2: samtools view exits 1 on duplicate BAM header entry (NC_002076.2)
- [2026-07-06] L-3: covid_viralscan/results/ is gitignored — SURVEY_SUMMARY.md not tracked

## By tag

- `scrna-seq`: L-6, L-7, L-8, L-9, L-14
- `anellovirus`: L-1, L-18, L-20, L-21
- `covid`: L-16, L-18, L-20, L-21
- `verify-by-artifact`: L-11, L-12, L-14, L-16
- `barcodes`: L-4, L-14, L-16
- `bioinformatics`: L-5, L-11, L-12
- `confounding`: L-6, L-7, L-8
- `eve`: L-18, L-20, L-21
- `gotcha`: L-11, L-13, L-14
- `host-homology`: L-1, L-18, L-20
- `bustools`: L-12, L-14
- `causal-inference`: L-7, L-8
- `cellranger`: L-4, L-16
- `covid_viralscan`: L-3, L-4
- `git`: L-3, L-5
- `gitignore`: L-3, L-5
- `kallisto`: L-11, L-12
- `kb-python`: L-11, L-12
- `mycelium`: L-10, L-19
- `review`: L-6, L-9
- `snakemake`: L-12, L-16
- `summarize_survey`: L-3, L-4
- `whitelist`: L-12, L-14
- `10x`: L-14
- `5-prime`: L-14
- `NC_002076.2`: L-2
- `aav2`: L-21
- `aifi`: L-19
- `anndata`: L-22
- `bam-header`: L-2
- `blast`: L-21
- `bug`: L-10
- `bulk`: L-13
- `bulk-rnaseq`: L-1
- `cDNA-reference`: L-1
- `cell-calling`: L-16
- `cell-type-enrichment`: L-20
- `celltypist`: L-20
- `ci`: L-22
- `circularity`: L-6
- `classifier`: L-9
- `cluster`: L-17
- `convention-pack`: L-19
- `conventions`: L-10
- `count-data`: L-8
- `covariate`: L-6
- `cpu-limit`: L-17
- `cross-validation`: L-9
- `dependencies`: L-22
- `depth-confound`: L-15
- `duplicate-contig`: L-2
- `emcv`: L-21
- `empty-droplets`: L-14
- `emptydrops`: L-16
- `endogenous-viral-elements`: L-18
- `environment-yml`: L-22
- `evalue`: L-15
- `false-positive`: L-1
- `feature-selection`: L-9
- `gem-x`: L-14
- `gene-selection`: L-7
- `grchr38-homology`: L-21
- `gse128078`: L-13
- `gzip`: L-12
- `hostresponse`: L-15
- `hvg`: L-9
- `ingest`: L-5
- `kb-ref`: L-11
- `label-definition`: L-8
- `large-files`: L-5
- `leakage`: L-9
- `medium-partition`: L-17
- `mito`: L-6
- `mitochondrial`: L-7
- `ngs_tools`: L-11
- `no-deps`: L-22
- `packaging`: L-22
- `panel-screen`: L-21
- `partition`: L-17
- `paths`: L-13
- `phase-a-b`: L-21
- `plasma-cells`: L-20
- `pyproject`: L-22
- `qc`: L-7
- `qos`: L-17
- `quant`: L-16
- `redo`: L-16
- `reference-build`: L-11
- `reference-panel`: L-13
- `references`: L-5
- `release-hygiene`: L-22
- `results`: L-3
- `samtools`: L-2
- `scrna`: L-19
- `sentinel`: L-16
- `sequencing-depth`: L-8
- `set-e`: L-2
- `silent-correctness`: L-13
- `skill-installation`: L-19
- `slurm`: L-17
- `slurm-exit-code`: L-2
- `specificity`: L-1
- `star`: L-1
- `starsolo`: L-4
- `statistics`: L-6
- `synthetic-data`: L-15
- `testing`: L-15
- `thresholding`: L-8
- `tooling`: L-10
- `transposable-elements`: L-18
- `viralscan`: L-16
- `yaml`: L-10

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
