<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-10-10

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| HANDOFF_2026-09-29.md | 12 entries | 2026-09-29 | Read these first, not this file, State at handoff, The thing to decide next: F-017, Other open rows, Validating the index |
| HANDOFF_2026-09-30.md | 14 entries | 2026-10-05 | 0. Orientation (read first), 1. State at a glance, 2. Environments and commands, Governance re-pin (needed whenever a pinned file changes), 3. What this session did (commits, oldest first) |
| HANDOFF_2026-10-01.md | 8 entries | 2026-10-05 | 1. Running right now, 2. Done this session, MECH-A closed (Virus Identity table), Datasets the user supplied: scoped (3 parallel agents), EXPL-HPV16, done: the HPV16 positive control passes |
| HANDOFF_2026-10-03.md | 5 entries | 2026-10-05 | Closed in this pass, Advanced, still `[~]`, Running, Waiting on user decisions, Gotchas found |
| conventions.md | 2 sections | 2026-07-06 | scRNA-seq associations: control depth AND %mito, and define labels depth-independently, Cross-validation: fit feature selection inside the split |
| decisions.md | 43 entries (large — read selectively) | 2026-10-09 | 3.0 default-selection design settled by user grill, CMP-06 compares implementations on one reference, DEF-03 contradiction = resolved viral set vs manifest viral set, `--strand` stays opt-in; new manifest options omitted when unset, Parallel implementers with fixed file ownership |
| last-session.md | 24 entries | 2026-07-15 | 2026-07-15 (post-compaction) — T5, T6, SH2.6 closed; EVE job 25237061 running, Pending (as of session end), 2026-07-15 — EVE (Endogenous Viral Element) characterisation analysis, 2026-07-15 — aifi-scrna-pipeline skill pack installed, Pending (as of session end) |
| learnings.md | 63 entries (large — read selectively) | 2026-10-09 | cDNA-only host reference causes false-positive viral signal from GRCh38 non-coding reads, samtools view exits 1 on duplicate BAM header entry (NC_002076.2), covid_viralscan/results/ is gitignored — SURVEY_SUMMARY.md not tracked, summarize_survey.py --cellranger-outs skipped: script expects one barcode set for all samples, 64 GB references/ was untracked but NOT gitignored |
| log/ | 80 sessions | 2026-10-10 | viralscan (80) |
| findings/ | 4 findings across 28 topics | 2026-10-07 | twostep-hpv77-call-is-host-reads-star-missed, evidence-replay-used-primary-ec-numbering, hhv6a-residual-in-hhv6b-sample-is-6b, hhv6b-call-in-ebv-sample-is-telomere-repeat, per-sample-run-cost, +23 more |
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-10-10 (heuristic)

## Tag clusters

- **scrna-seq** (5 entries) — L-6, L-7, L-8, L-9, L-14
- **anellovirus** (4 entries) — L-1, L-18, L-20, L-21
- **covid** (4 entries) — L-16, L-18, L-20, L-21
- **verify-by-artifact** (4 entries) — L-11, L-12, L-14, L-16
- **barcodes** (3 entries) — L-4, L-14, L-16
- **bioinformatics** (3 entries) — L-5, L-11, L-12

## Most recent (10)

- [2026-10-09] L-63: `build_bundled_panel_ref.py` reconciliation can be run offline and stopped before `kb ref`
- [2026-10-09] D-43: `AB303562` recorded as the first `index_exclusions.tsv` row; builder step 4b fetches uncovered anellovirus rows (PANEL-01 WP4)
- [2026-10-09] D-42: The 3 no-CDS PANEL-01 records are modelled as single-exon CDS GTFs, not excluded (WP4)
- [2026-10-08] L-62: `viral_ref_final/build/viral.fa` is not in this checkout; synthetic-test pitfall
- [2026-10-08] D-41: Candidates with under 5 % of k-mers outside the panel are excluded as `kmer_twin_of` (PANEL-01 WP1b)
- [2026-10-08] D-40: WP1b k-mer sharing: group = sibling_group, else genus, else species; no NCBI fetch without NCBI_EMAIL (PANEL-01)
- [2026-10-07] D-39: One frozen emptyDrops call per sample for every cell-level analysis (CELLS-01)
- [2026-10-06] L-61: The old overlapping-window capture formula overstated recall by tens of points
- [2026-10-06] L-60: Stop hook flagged `prompt.MD` on a question-only session
- [2026-10-06] L-59: em-cell stores a dense all-gene theta per barcode and OOMs at full depth

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

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 "$(cat .mycelium/plugin-root)/skills/core/scripts/recall_lessons.py" --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
