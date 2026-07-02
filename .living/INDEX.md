<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-02

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 2 sections | 2026-07-02 | scRNA-seq associations: control depth AND %mito, and define labels depth-independently, Cross-validation: fit feature selection inside the split |
| decisions.md | 0 entries | 2026-07-02 | — |
| learnings.md | 7 entries | 2026-07-02 | 64 GB references/ was untracked but NOT gitignored, Controlling for a composite covariate that CONTAINS the tested feature is circular, Mitochondrial genes need a %mito control before claiming them as biology, Raw-count positivity thresholds silently confound with sequencing depth, Feature selection (HVG) before the CV split leaks into held-out metrics |
| log/ | 4 sessions | 2026-07-01 | viralscan (4) |
| findings/ | 4 findings across 3 topics | 2026-07-02 | host-transcriptome-encodes-viral-infection-state, raw-count-thresholds-confound-with-sequencing-depth, unsupervised-feature-selection-leakage-is-often-negligible |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-02 (heuristic)

## Tag clusters

- **scrna-seq** (4 entries) — L-2, L-3, L-4, L-5
- **confounding** (3 entries) — L-2, L-3, L-4
- **bioinformatics** (2 entries) — L-1, L-7
- **causal-inference** (2 entries) — L-3, L-4
- **review** (2 entries) — L-2, L-5

## Most recent (10)

- [2026-07-02] L-2: Controlling for a composite covariate that CONTAINS the tested feature is circular
- [2026-07-02] L-3: Mitochondrial genes need a %mito control before claiming them as biology
- [2026-07-02] L-7: kb ref / kallisto index have two silent FASTA-vs-GTF contracts
- [2026-07-01] L-1: 64 GB references/ was untracked but NOT gitignored
- [2026-07-01] L-4: Raw-count positivity thresholds silently confound with sequencing depth
- [2026-07-01] L-5: Feature selection (HVG) before the CV split leaks into held-out metrics
- [2026-07-01] L-6: ACTIVE_CONVENTIONS.yaml is malformed after install_convention.py

## By tag

- `scrna-seq`: L-2, L-3, L-4, L-5
- `confounding`: L-2, L-3, L-4
- `bioinformatics`: L-1, L-7
- `causal-inference`: L-3, L-4
- `review`: L-2, L-5
- `bug`: L-6
- `circularity`: L-2
- `classifier`: L-5
- `conventions`: L-6
- `count-data`: L-4
- `covariate`: L-2
- `cross-validation`: L-5
- `feature-selection`: L-5
- `gene-selection`: L-3
- `git`: L-1
- `gitignore`: L-1
- `gotcha`: L-7
- `hvg`: L-5
- `ingest`: L-1
- `kallisto`: L-7
- `kb-python`: L-7
- `kb-ref`: L-7
- `label-definition`: L-4
- `large-files`: L-1
- `leakage`: L-5
- `mito`: L-2
- `mitochondrial`: L-3
- `mycelium`: L-6
- `ngs_tools`: L-7
- `qc`: L-3
- `reference-build`: L-7
- `references`: L-1
- `sequencing-depth`: L-4
- `statistics`: L-2
- `thresholding`: L-4
- `tooling`: L-6
- `verify-by-artifact`: L-7
- `yaml`: L-6

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
