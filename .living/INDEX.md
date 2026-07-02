<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-07-02

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 2 sections | 2026-07-02 | scRNA-seq associations: control depth AND %mito, and define labels depth-independently, Cross-validation: fit feature selection inside the split |
| decisions.md | 0 entries | 2026-07-02 | — |
| learnings.md | 6 entries | 2026-07-02 | 64 GB references/ was untracked but NOT gitignored, Mitochondrial genes need a %mito control before claiming them as biology, Raw-count positivity thresholds silently confound with sequencing depth, Feature selection (HVG) before the CV split leaks into held-out metrics, ACTIVE_CONVENTIONS.yaml is malformed after install_convention.py |
| log/ | 4 sessions | 2026-07-01 | viralscan (4) |
| findings/ | 4 findings across 3 topics | 2026-07-02 | host-transcriptome-encodes-viral-infection-state, raw-count-thresholds-confound-with-sequencing-depth, unsupervised-feature-selection-leakage-is-often-negligible |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-07-02 (heuristic)

## Tag clusters

- **scrna-seq** (3 entries) — L-2, L-3, L-4
- **bioinformatics** (2 entries) — L-1, L-6
- **causal-inference** (2 entries) — L-2, L-3
- **confounding** (2 entries) — L-2, L-3

## Most recent (10)

- [2026-07-02] L-2: Mitochondrial genes need a %mito control before claiming them as biology
- [2026-07-02] L-6: kb ref / kallisto index have two silent FASTA-vs-GTF contracts
- [2026-07-01] L-1: 64 GB references/ was untracked but NOT gitignored
- [2026-07-01] L-3: Raw-count positivity thresholds silently confound with sequencing depth
- [2026-07-01] L-4: Feature selection (HVG) before the CV split leaks into held-out metrics
- [2026-07-01] L-5: ACTIVE_CONVENTIONS.yaml is malformed after install_convention.py

## By tag

- `scrna-seq`: L-2, L-3, L-4
- `bioinformatics`: L-1, L-6
- `causal-inference`: L-2, L-3
- `confounding`: L-2, L-3
- `bug`: L-5
- `classifier`: L-4
- `conventions`: L-5
- `count-data`: L-3
- `cross-validation`: L-4
- `feature-selection`: L-4
- `gene-selection`: L-2
- `git`: L-1
- `gitignore`: L-1
- `gotcha`: L-6
- `hvg`: L-4
- `ingest`: L-1
- `kallisto`: L-6
- `kb-python`: L-6
- `kb-ref`: L-6
- `label-definition`: L-3
- `large-files`: L-1
- `leakage`: L-4
- `mitochondrial`: L-2
- `mycelium`: L-5
- `ngs_tools`: L-6
- `qc`: L-2
- `reference-build`: L-6
- `references`: L-1
- `review`: L-4
- `sequencing-depth`: L-3
- `thresholding`: L-3
- `tooling`: L-5
- `verify-by-artifact`: L-6
- `yaml`: L-5

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
