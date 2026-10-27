<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-10-26
Last audit: 2026-10-27

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| HANDOFF_2026-09-29.md | 12 entries | 2026-10-10 | Read these first, not this file, State at handoff, The thing to decide next: F-017, Other open rows, Validating the index |
| HANDOFF_2026-09-30.md | 14 entries | 2026-10-05 | 0. Orientation (read first), 1. State at a glance, 2. Environments and commands, Governance re-pin (needed whenever a pinned file changes), 3. What this session did (commits, oldest first) |
| HANDOFF_2026-10-01.md | 8 entries | 2026-10-05 | 1. Running right now, 2. Done this session, MECH-A closed (Virus Identity table), Datasets the user supplied: scoped (3 parallel agents), EXPL-HPV16, done: the HPV16 positive control passes |
| HANDOFF_2026-10-03.md | 5 entries | 2026-10-05 | Closed in this pass, Advanced, still `[~]`, Running, Waiting on user decisions, Gotchas found |
| conventions.md | 2 sections | 2026-07-06 | scRNA-seq associations: control depth AND %mito, and define labels depth-independently, Cross-validation: fit feature selection inside the split |
| decisions.md | 81 entries (large — read selectively) | 2026-10-25 | 3.0 default-selection design settled by user grill, CMP-06 compares implementations on one reference, DEF-03 contradiction = resolved viral set vs manifest viral set, `--strand` stays opt-in; new manifest options omitted when unset, Parallel implementers with fixed file ownership |
| last-session.md | 24 entries | 2026-07-15 | 2026-07-15 (post-compaction) — T5, T6, SH2.6 closed; EVE job 25237061 running, Pending (as of session end), 2026-07-15 — EVE (Endogenous Viral Element) characterisation analysis, 2026-07-15 — aifi-scrna-pipeline skill pack installed, Pending (as of session end) |
| learnings.md | 77 entries (large — read selectively) | 2026-10-26 | cDNA-only host reference causes false-positive viral signal from GRCh38 non-coding reads, samtools view exits 1 on duplicate BAM header entry (NC_002076.2), covid_viralscan/results/ is gitignored — SURVEY_SUMMARY.md not tracked, summarize_survey.py --cellranger-outs skipped: script expects one barcode set for all samples, 64 GB references/ was untracked but NOT gitignored |
| log/ | 90 sessions | 2026-10-26 | viralscan (89), translated-proteins (1) |
| decisions.md | 81 entries (large — read selectively) | 2026-10-27 | 3.0 default-selection design settled by user grill, CMP-06 compares implementations on one reference, DEF-03 contradiction = resolved viral set vs manifest viral set, `--strand` stays opt-in; new manifest options omitted when unset, Parallel implementers with fixed file ownership |
| last-session.md | 24 entries | 2026-07-15 | 2026-07-15 (post-compaction) — T5, T6, SH2.6 closed; EVE job 25237061 running, Pending (as of session end), 2026-07-15 — EVE (Endogenous Viral Element) characterisation analysis, 2026-07-15 — aifi-scrna-pipeline skill pack installed, Pending (as of session end) |
| learnings.md | 76 entries (large — read selectively) | 2026-10-27 | cDNA-only host reference causes false-positive viral signal from GRCh38 non-coding reads, samtools view exits 1 on duplicate BAM header entry (NC_002076.2), covid_viralscan/results/ is gitignored — SURVEY_SUMMARY.md not tracked, summarize_survey.py --cellranger-outs skipped: script expects one barcode set for all samples, 64 GB references/ was untracked but NOT gitignored |
| log/ | 85 sessions | 2026-10-27 | viralscan (84), translated-proteins (1) |
| findings/ | 4 findings across 31 topics | 2026-10-25 | spikein-flank-index-retired, hhv6b-flank-reads-missed-by-cds-only-index, hsv1-lat-exon-recut, twostep-hpv77-call-is-host-reads-star-missed, evidence-replay-used-primary-ec-numbering, +26 more |
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-10-26 (heuristic)
Last summarized: 2026-10-27 (heuristic)

## Tag clusters

- **performance** (7 entries) — L-72, L-73, L-74, L-75, L-76
- **scrna-seq** (5 entries) — L-6, L-7, L-8, L-9, L-14
- **anellovirus** (4 entries) — L-1, L-18, L-20, L-21
- **covid** (4 entries) — L-16, L-18, L-20, L-21
- **verify-by-artifact** (4 entries) — L-11, L-12, L-14, L-16
- **anndata** (3 entries) — L-22, L-68, L-74

## Most recent (10)

- [2026-10-26] L-77: st_dev is a per-host mount id; use inode alone to prove a rename
- [2026-10-21] D-81: Detection-limit spike-in design approved (user: "approve all"; Q1-Q4, Q6-Q8)
- [2026-10-21] D-80: STAR-03: star-index layout and manifest
- [2026-10-21] D-79: STAR-02: star_hostvirus hardening choices
- [2026-10-21] D-78: Discovery fingerprint and resume
- [2026-10-21] D-77: Discovery candidates: observation state only
- [2026-10-21] D-76: Discovery polyA/T trim before Kraken2
- [2026-10-21] D-75: Discovery per-UMI rule: deepest-consistent
- [2026-10-21] D-74: Discovery Kraken2 defaults: Cellsnake 0.05 / 4
- [2026-10-18] L-76: Perf wave 2 C00: two checkouts at different commits cannot be byte-compared without excusing the checkout sha
- [2026-10-18] L-75: Perf wave 2 B18: a fused reducer can be a wash on its own, and tracemalloc invents CPU regressions

## By tag

- `performance`: L-69, L-70, L-72, L-73, L-74, L-75, L-76
- `scrna-seq`: L-6, L-7, L-8, L-9, L-14
- `anellovirus`: L-1, L-18, L-20, L-21
- `covid`: L-16, L-18, L-20, L-21
- `verify-by-artifact`: L-11, L-12, L-14, L-16
- `anndata`: L-22, L-68, L-74
- `barcodes`: L-4, L-14, L-16
- `bioinformatics`: L-5, L-11, L-12
- `confounding`: L-6, L-7, L-8
- `eve`: L-18, L-20, L-21
- `gotcha`: L-11, L-13, L-14
- `host-homology`: L-1, L-18, L-20
- `testing`: L-15, L-76, L-77
- `bustools`: L-12, L-14
- `byte-identical`: L-72, L-73
- `causal-inference`: L-7, L-8
- `cellranger`: L-4, L-16
- `ci`: L-22, L-68
- `covid_viralscan`: L-3, L-4
- `git`: L-3, L-5
- `gitignore`: L-3, L-5
- `governance`: L-69, L-70
- `hpc`: L-71, L-77
- `io`: L-74, L-75
- `kallisto`: L-11, L-12
- `kb-python`: L-11, L-12
- `mycelium`: L-10, L-19
- `mypy`: L-69, L-73
- `review`: L-6, L-9
- `snakemake`: L-12, L-16
- `summarize_survey`: L-3, L-4
- `testing`: L-15, L-76
- `whitelist`: L-12, L-14
- `10x`: L-14
- `5-prime`: L-14
- `CI-drift`: L-73
- `NC_002076.2`: L-2
- `aav2`: L-21
- `aifi`: L-19
- `annotation`: L-73
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
- `convention-pack`: L-19
- `conventions`: L-10
- `count-data`: L-8
- `covariate`: L-6
- `cpu-limit`: L-17
- `cross-validation`: L-9
- `dependencies`: L-22
- `depth-confound`: L-15
- `documentation-drift`: L-70
- `duplicate-contig`: L-2
- `emcv`: L-21
- `empty-droplets`: L-14
- `emptydrops`: L-16
- `endogenous-viral-elements`: L-18
- `environment-yml`: L-22
- `evalue`: L-15
- `evidence`: L-72
- `exactness`: L-76
- `false-positive`: L-1
- `feature-selection`: L-9
- `fixture-blindspot`: L-72
- `gem-x`: L-14
- `gene-selection`: L-7
- `grchr38-homology`: L-21
- `gse128078`: L-13
- `gzip`: L-12
- `hostresponse`: L-15
- `hpc`: L-71
- `hvg`: L-9
- `index-space`: L-72
- `ingest`: L-5
- `kb-ref`: L-11
- `label-definition`: L-8
- `large-files`: L-5
- `leakage`: L-9
- `lint`: L-73
- `macos`: L-68
- `medium-partition`: L-17
- `memory`: L-75
- `migration`: L-77
- `mito`: L-6
- `mitochondrial`: L-7
- `multimapping`: L-72
- `ngs_tools`: L-11
- `no-deps`: L-22
- `numpy-stubs`: L-73
- `packaging`: L-22
- `panel-screen`: L-21
- `partition`: L-17
- `paths`: L-13
- `phase-a-b`: L-21
- `plasma-cells`: L-20
- `process-safety`: L-71
- `profiling`: L-75
- `provenance`: L-76
- `pyproject`: L-22
- `qc`: L-7
- `qos`: L-17
- `quant`: L-16
- `rebase`: L-72
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
- `shared-account`: L-71
- `silent-correctness`: L-13
- `skill-installation`: L-19
- `slurm`: L-17
- `slurm-exit-code`: L-2
- `specificity`: L-1
- `star`: L-1
- `starsolo`: L-4
- `statistics`: L-6
- `synthetic-data`: L-15
- `thresholding`: L-8
- `timing-noise`: L-69
- `tooling`: L-10
- `transposable-elements`: L-18
- `validation`: L-74
- `viralscan`: L-16
- `worktrees`: L-70
- `yaml`: L-10

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 "$(cat .mycelium/plugin-root)/skills/core/scripts/recall_lessons.py" --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
