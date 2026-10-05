# ViralScan 3.0 Correctness, Validation, and Publication Plan

Status: authoritative v3 roadmap, accepted 2026-07-22. Earlier PLAN, PROGRESS,
ROADMAP, publication-readiness, and benchmark documents describe pre-v3 work and
are historical unless this document explicitly incorporates them.

Execution status, dependencies, acceptance commands, and the single next action
are tracked in [`../../PLAN.md`](../../PLAN.md). This document defines the locked
contracts; `PLAN.md` is the operational checklist.

## Release objective and locked contracts

ViralScan 3.0 is a correctness-first major release. Pre-v3 corrected or combined
counts are scientifically incompatible and must be regenerated from retained
corrected BUS/reference inputs or FASTQ; they must never be numerically migrated.

- Production scope is human paired-end droplet scRNA-seq.
- The count unit is one resolved corrected CB-UMI molecule, never BUS record/read
  multiplicity. EC transcripts are projected to distinct genes before ambiguity
  classification. Multiple ECs for one CB-UMI are resolved explicitly; disjoint
  collisions remain unresolved.
- `adata.X` is the complete selected-method molecule matrix. Non-overlapping
  `counts_unique` and `counts_ambiguous_allocated` layers sum to `X`.
- `host-conservative` is primary and preserves mixed host-virus molecule mass by
  allocating it only among compatible host genes. Other methods are `equal`,
  `unique-weighted`, `em-global`, and hierarchical `em-cell`; plain `em` is removed.
- Host homology, sibling-virus ambiguity, unresolved mass, and allocation model
  scope remain visible in outputs and evidence.
- Cell calling is `auto`: an external list when supplied, otherwise EmptyDrops.
  Knee calling is explicit and never a fallback.
- A nonzero molecule is candidate evidence. Probable/strong tiers require
  calibrated read-level specificity and QC.
- Default reference/workflow is curated combined human-plus-virus without the
  expanded anellovirus panel. Anellovirus is opt-in and screening-only without
  full-host-genome competition.
- Run reuse is safe: default refuses non-empty output, `--resume` requires an
  identical fingerprint, and `--overwrite` is mutually exclusive and explicit.
- Pip supplies Python/reporting/validation/CLI parsing. Locked conda and digest-
  pinned containers guarantee the full workflow.

## Ordered execution

1. Freeze unsafe releases and legacy claims; inventory artifacts; enforce data
   governance; create a machine-readable claim registry.
2. Specify algorithms, versioned schemas, invariants, and a hand-verifiable golden
   fixture before extending the counting engine.
3. Replace the multimapping engine with streaming corrected CB-UMI resolution,
   distinct-gene EC projection, exact conservation, deterministic allocation,
   explicit unresolved evidence, diagnostics, and bounded sparse accumulation.
4. Rebuild downstream consumers so their selected layer, method, fingerprints,
   cell calling, and evidence rules are explicit; make reruns create a fresh,
   fully self-consistent result tree.
5. Freeze reproducible reference profiles, fail closed on fetches, make expanded
   anellovirus opt-in, record accession/sequence/licence provenance, detect
   duplicates, and add full-genome homology controls.
6. Keep combined quantification recommended. Reimplement two-step filtering with
   exact fragment identity, paired-read synchronization, loss diagnostics, and no
   CB-UMI-wide deletion.
7. Make `viralscan evidence` accession-specific and end-to-end: exact lineage,
   competitive host-virus alignment/BLAST, deterministic sampling, defined
   deduplication, raw/deduplicated BAM and coverage, IGV assets, per-virus/per-cell
   QC, and diagnostic interpretation flags.
8. Split CLI/service architecture; add `doctor` and `validate-run`; preflight all
   inputs/tools/space; stage outputs atomically; fingerprint sentinels and capture
   structured logs, exit codes, and tool versions.
9. Complete wheel/sdist, locked conda, Docker and Apptainer parity; pin assets and
   toolchains; gate exact-SHA releases; generate SBOM/licence reports; finish
   governance and package metadata.
10. Pre-register and execute a deterministic truth panel spanning abundance,
    infected-cell fraction, homology, sibling viruses, complexity, PCR duplication,
    UMI collisions, ambient contamination, and supported chemistries. Harmonize
    STARsolo, traditional alignment, Venus, and Viral-Track comparisons by inputs,
    references, barcodes, features, denominators, and count layers.
11. Rebuild documentation and executable vignettes from actual v3 contracts,
    including negative and ambiguous cases and balanced workflow comparisons.
12. Regenerate every manuscript number/figure from frozen v3 artifacts; reassess
    host response; keep TTV as a homology case study absent orthogonal validation;
    complete ethics/authorship/citation and independent claim audits.
13. Publish rc1, obtain three-laboratory external testing, resolve correctness and
    installation failures, freeze scientific artifacts, tag stable only from the
    exact green SHA, archive, submit, and maintain scientific-incident triage.

## Acceptance and gates

Core tests cover molecule conservation, PCR duplicates, multiple transcripts and
ECs, collisions, host-virus and sibling-virus ambiguity, zero support, convergence,
determinism, input ordering, and chunk size. CLI safety, real integrations,
references, documentation, packaging, containers, and the 103-million-record EBV
performance baseline each have explicit gates.

No release may proceed without molecule correctness, safe output semantics,
reference provenance, end-to-end evidence, complete environments, real integration
tests, and truthful docs. Stable 3.0 additionally requires rc testing, truth-panel
holdout, harmonized comparison, anellovirus safeguards, and reproducible artifacts.
The manuscript additionally requires Venus and Viral-Track, uncertainty analysis,
public provenance, ethics/authorship, and a claim-to-artifact audit.

Future non-blockers are nonhuman hosts; bulk, spatial, long-read and single-nucleus
models; cluster-informed EM; variant/de-novo discovery; cloud/distributed execution;
orthogonal wet-lab integration; and curator-reviewed panel update automation.
