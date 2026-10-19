<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-10-19

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| HANDOFF_2026-09-29.md | 12 entries | 2026-10-10 | Read these first, not this file, State at handoff, The thing to decide next: F-017, Other open rows, Validating the index |
| HANDOFF_2026-09-30.md | 14 entries | 2026-10-05 | 0. Orientation (read first), 1. State at a glance, 2. Environments and commands, Governance re-pin (needed whenever a pinned file changes), 3. What this session did (commits, oldest first) |
| HANDOFF_2026-10-01.md | 8 entries | 2026-10-05 | 1. Running right now, 2. Done this session, MECH-A closed (Virus Identity table), Datasets the user supplied: scoped (3 parallel agents), EXPL-HPV16, done: the HPV16 positive control passes |
| HANDOFF_2026-10-03.md | 5 entries | 2026-10-05 | Closed in this pass, Advanced, still `[~]`, Running, Waiting on user decisions, Gotchas found |
| conventions.md | 2 sections | 2026-07-06 | scRNA-seq associations: control depth AND %mito, and define labels depth-independently, Cross-validation: fit feature selection inside the split |
| decisions.md | 72 entries (large — read selectively) | 2026-10-19 | 3.0 default-selection design settled by user grill, CMP-06 compares implementations on one reference, DEF-03 contradiction = resolved viral set vs manifest viral set, `--strand` stays opt-in; new manifest options omitted when unset, Parallel implementers with fixed file ownership |
| last-session.md | 24 entries | 2026-07-15 | 2026-07-15 (post-compaction) — T5, T6, SH2.6 closed; EVE job 25237061 running, Pending (as of session end), 2026-07-15 — EVE (Endogenous Viral Element) characterisation analysis, 2026-07-15 — aifi-scrna-pipeline skill pack installed, Pending (as of session end) |
| learnings.md | 94 entries (large — read selectively) | 2026-10-19 | cDNA-only host reference causes false-positive viral signal from GRCh38 non-coding reads, samtools view exits 1 on duplicate BAM header entry (NC_002076.2), covid_viralscan/results/ is gitignored — SURVEY_SUMMARY.md not tracked, summarize_survey.py --cellranger-outs skipped: script expects one barcode set for all samples, 64 GB references/ was untracked but NOT gitignored |
| log/ | 101 sessions | 2026-10-19 | viralscan (100), translated-proteins (1) |
| findings/ | 4 findings across 29 topics | 2026-10-12 | hsv1-lat-exon-recut, twostep-hpv77-call-is-host-reads-star-missed, evidence-replay-used-primary-ec-numbering, hhv6a-residual-in-hhv6b-sample-is-6b, hhv6b-call-in-ebv-sample-is-telomere-repeat, +24 more |
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-10-19 (heuristic)

## Tag clusters

- **trackA** (19 entries) — D-68, D-69, D-70, D-71, D-72
- **anello** (7 entries) — L-92, L-93, D-68, D-69, D-70
- **hostresponse** (5 entries) — L-15, L-86, L-90, L-94, D-67
- **scrna-seq** (5 entries) — L-6, L-7, L-8, L-9, L-14
- **anellovirus** (4 entries) — L-1, L-18, L-20, L-21
- **covid** (4 entries) — L-16, L-18, L-20, L-21

## Most recent (10)

- [2026-10-19] L-94: Cross-library hostresponse controls: hsv1 has matched 37C uninfected libs, hhv6b has one clean control and no batch null
- [2026-10-19] D-72: Grilling round 3 approved; Track A handoff refreshed
- [2026-10-19] L-93: Real viral 3' reads rarely carry poly-A tails; planter default (1/3) was far too high
- [2026-10-19] D-71: Grilling rounds 1-2 approved (user: "i approve" x2)
- [2026-10-19] L-92: Anello plant scorer keyed on raw CB: hnscc recovery under-counted; tail geometry caps 3' datasets
- [2026-10-19] D-70: hhv6b real-library anello array released; scorer fixed (no threshold change)
- [2026-10-19] L-91: Anello v2 plant scoring: only hhv6b passes the STAR exact@100 gate
- [2026-10-19] D-69: Anello real-library array held after plant gates
- [2026-10-19] L-90: hostresponse `--label cpm` picks the wrong cells when viral UMI are sparse; clean controls are scarce in infected libraries
- [2026-10-19] D-68: Track A anellovirus settle: plan approved, PREREG frozen, tools built

## By tag

- `trackA`: L-81, L-82, L-83, L-84, L-85, L-86, L-89, L-90, L-91, L-92, L-93, L-94, D-66, D-67, D-68, D-69, D-70, D-71, D-72
- `anello`: L-89, L-91, L-92, L-93, D-68, D-69, D-70
- `hostresponse`: L-15, L-86, L-90, L-94, D-67
- `scrna-seq`: L-6, L-7, L-8, L-9, L-14
- `anellovirus`: L-1, L-18, L-20, L-21
- `covid`: L-16, L-18, L-20, L-21
- `verify-by-artifact`: L-11, L-12, L-14, L-16
- `barcodes`: L-4, L-14, L-16
- `bioinformatics`: L-5, L-11, L-12
- `confounding`: L-6, L-7, L-8
- `deletion-gates`: L-83, L-84, L-85
- `eve`: L-18, L-20, L-21
- `gotcha`: L-11, L-13, L-14
- `host-homology`: L-1, L-18, L-20
- `anndata`: L-22, L-68
- `bustools`: L-12, L-14
- `causal-inference`: L-7, L-8
- `cellranger`: L-4, L-16
- `ci`: L-22, L-68
- `covid_viralscan`: L-3, L-4
- `deletion`: L-82, D-66
- `git`: L-3, L-5
- `gitignore`: L-3, L-5
- `kallisto`: L-11, L-12
- `kb-python`: L-11, L-12
- `mycelium`: L-10, L-19
- `planter`: L-92, L-93
- `provenance`: L-83, L-85
- `review`: L-6, L-9
- `slurm`: L-17, L-81
- `snakemake`: L-12, L-16
- `summarize_survey`: L-3, L-4
- `whitelist`: L-12, L-14
- `10x`: L-14
- `5-prime`: L-14
- `NC_002076.2`: L-2
- `aav2`: L-21
- `aifi`: L-19
- `argparse`: L-68
- `bam-header`: L-2
- `blast`: L-21
- `bug`: L-10
- `bulk`: L-13
- `bulk-rnaseq`: L-1
- `cDNA-reference`: L-1
- `cell-calling`: L-16
- `cell-type-enrichment`: L-20
- `celltypist`: L-20
- `circularity`: L-6
- `classifier`: L-9
- `cluster`: L-17
- `controls`: L-94
- `convention-pack`: L-19
- `conventions`: L-10
- `count-data`: L-8
- `covariate`: L-6
- `cpu-limit`: L-17
- `cross-validation`: L-9
- `deletion-window`: L-81
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
- `gates`: L-91
- `gem-x`: L-14
- `gene-selection`: L-7
- `governance`: L-80
- `grchr38-homology`: L-21
- `grilling`: D-71
- `gse128078`: L-13
- `gzip`: L-12
- `handoff`: D-72
- `hhv6b`: L-82
- `hvg`: L-9
- `ingest`: L-5
- `inventory`: L-80
- `kb-ref`: L-11
- `label-definition`: L-8
- `labels`: L-90
- `large-files`: L-5
- `leakage`: L-9
- `macos`: L-68
- `medium-partition`: L-17
- `mito`: L-6
- `mitochondrial`: L-7
- `ngs_tools`: L-11
- `no-deps`: L-22
- `normalization`: L-86
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
- `scorer`: L-92
- `scrna`: L-19
- `sentinel`: L-16
- `sequencing-depth`: L-8
- `set-e`: L-2
- `silent-correctness`: L-13
- `skill-installation`: L-19
- `slurm-exit-code`: L-2
- `specificity`: L-1
- `star`: L-1
- `starsolo`: L-4
- `statistics`: L-6
- `synthetic-data`: L-15
- `taxonomy`: L-80
- `testing`: L-15
- `thresholding`: L-8
- `tonsil`: D-66
- `tooling`: L-10
- `transposable-elements`: L-18
- `validation`: L-89
- `viralscan`: L-16
- `yaml`: L-10

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 "$(cat .mycelium/plugin-root)/skills/core/scripts/recall_lessons.py" --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
