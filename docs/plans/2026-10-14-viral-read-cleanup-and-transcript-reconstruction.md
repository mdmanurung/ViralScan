# Viral read cleanup and transcript reconstruction: reviewed proposal

Date: 2026-10-14
Last updated: 2026-10-14
Project: ViralScan
Reviewed checkout: `codex/drop-py39`, HEAD `5e551c9`; working tree has existing changes.
Status: exploratory implementation in progress; handoff requested before the corrected data run.

## Objective and scope

Evaluate whether conservative cDNA cleanup reduces reagent/tail artefacts while retaining
genuine viral molecules. Separately assess reference-based viral transcript reconstruction and
targeted host-virus junction reporting. These are proposed additions, not replacements for the
current detection, reference-masking, ambiguity, or evidence contracts.

The initial review authorized no execution. The user subsequently said "run it" on 2026-10-14,
authorizing this bounded exploratory prototype and pilot. Publication and changes to scientific
gates remain outside scope; frozen DSR rounds are read-only.
`PLAN.md` remains authoritative. Proposed work enters its queue only at implementation scoping;
the existing dirty PLAN and other user changes are preserved.

## Review findings and resolutions

| Finding | Resolution |
|---|---|
| The draft treats generic paired-end trimming as applicable to single-cell reads. | Restrict the first pilot to verified R1 barcode/UMI + R2 cDNA layouts. Keep retained R1 records unchanged. Reject unsupported layouts rather than guess. |
| Cleanup before host subtraction conflicts with today's routing. | Start with a separate preprocessing pilot. Do not reorder existing `--read-filter` modes. Package integration must explicitly wire the new input pair into both host alignment and counting. |
| Poly-A presence is not itself an artefact. | Retain the templated body, restrict trimming to configured read ends, and test legitimate terminal/internal A/T-rich viral sequences. Tail origin cannot always be distinguished from sequence alone. |
| A positive and presumed negative cannot establish general sensitivity/specificity. | Use exact tail-aware planted truth plus unplanted twins for mechanism assessment; use real positives/negatives as contextual checks with explicit truth tiers. Broader qualification remains separate. |
| The existing SENS-07 planter does not model poly-A read ends. | Extend or supplement the existing planting workflow with tail/orientation cases before using it to evaluate trimming. Its current internal-fragment results cannot qualify tail trimming. |
| Existing evidence BAMs are unsuitable for general spliced reconstruction. | Add a separate splice-aware alignment for assembly. Preserve the existing `minimap2 -ax sr` evidence route and its outputs. |
| Detection-selected assembly is circular and incomplete. | Label it candidate reconstruction; keep host/sibling ambiguity and selection metadata. A later discovery mode needs broader raw-read recruitment and its own validation. |
| Pooling libraries or deduplicating by start position can distort transcripts. | Assemble per biological sample and accession/strain hypothesis. Preserve library namespaces and distinct read evidence; quantify independent CB/UMI support separately. |
| RNA junctions do not establish DNA integration. | Report candidate host-virus RNA junctions, require molecule-level support and artefact review, and reserve genomic-integration claims for independent evidence. |
| Upstream thresholds and preprocessing descriptions are not interchangeable with code. | Pin source/version/commands before implementation. Do not copy published detection cutoffs or claim default Trim Galore removes every homopolymer. |

All findings above are addressed in the design. Numerical promotion margins, pilot inputs, and
the supported chemistry/strand recipes remain unresolved and must be fixed before outcome runs.

## Current code and reuse boundaries

- `src/viralscan/scripts/read_filter.py`: `artefact` drops selected pairs; `tso-trim` cuts a
  leading TSO from R2 and drops pairs with fewer than 31 remaining bases. Both are opt-in.
  This is not general adapter/poly-A trimming. Preserve their current semantics.
- `src/viralscan/Snakefile`, `runconfig.py`: the current order is raw reads, optional host
  subtraction, optional read filter, then `kb_count`. `kb_r1`/`kb_r2` also drive evidence consumers.
- `src/viralscan/scripts/host_filter.py`: exact fragment IDs and pair-validation helpers exist.
  Reuse them and the lineage conventions; never subtract all fragments sharing a CB/UMI.
- `src/viralscan/evidence.py`: exact BUS-flag replay, extracted-read lineage, competitive evidence,
  coverage, UMI-aware support, and browser BAMs exist. Its viral alignment uses `minimap2 -ax sr`.
  Assembly needs a separate spliced alignment and the full candidate read set, not a BLAST sample.
- `scripts/sens07_planted_sensitivity.py`: the current planter uses internal, forward-strand
  fragments and explicitly excludes poly-A end realism. Extend existing truth tooling rather
  than create a second general benchmark framework.
- `environment.yml` and its lock files contain STAR and samtools but do not declare Cutadapt or
  StringTie. Optional prototype tools need a pinned isolated runtime. No dependency was added.

## Progress summary

- [x] REVIEW: inspect current routing, evidence, truth limitations, and upstream source.
- [ ] VR-01: freeze a bounded pilot contract and source receipts — in progress; v3 pilot recipe frozen.
- [ ] VR-02: prepare and verify a standalone cleanup prototype — in progress; nine tests pass, lint remains open.
- [ ] VR-03: evaluate cleanup and decide package integration — in progress; v2 failed, corrected v3 dry-run passed.
- [ ] VR-04: assess optional viral transcript reconstruction — in progress; native fixture passed, real samples pending.
- [ ] VR-05: scope candidate host-virus junction reporting — not started; deferred.

## Implementation steps

### VR-01: freeze the pilot contract

- **Status:** not started.
- **Outcome:** one written, versioned recipe and evaluation manifest before outcome inspection.
- **Actions:** choose checksum-pinned 3-prime and 5-prime input layouts if both are to be
  supported; resolve cDNA orientation from chemistry/library metadata; choose adapter/tail
  sequences, minimum overlap/error handling, and minimum usable body length. For kallisto arms,
  use the recorded index k rather than assume every index uses 31. Treat poly-G handling as a
  separate instrument-specific decision. Do not add general trimming of all A/C/G/T runs.
- **Comparison arms:** baseline with cleanup off; existing `tso-trim`; one explicit cDNA-cleanup
  recipe. Within each chemistry use the same input fragments, references, strand, cell universe,
  allocation method, and detection settings. Record the changed preprocessing order. The pilot
  cleanup recipe must specify TSO handling so effects are not attributed to tail trimming alone.
- **Dependencies:** recover current guards and execution eligibility; select inputs and truth
  tiers; predeclare acceptable viral-molecule loss and required artefact reduction per target and
  chemistry. No universal margin is inferred from REF-06 or upstream defaults.
- **Affected areas:** future `PLAN.md` task entry, existing planting/protocol records, isolated
  runtime receipt. Pin Venus/Viral-Track revisions and trimming/assembly versions.
- **Validation:** inputs and recipes resolve; every comparison has matching denominators;
  tuning uses training-eligible data only. Held-out data cannot determine the trimming recipe.
- **Evidence:** pending input manifest, source/runtime hashes, acceptance margins, execution status.

### VR-02: prepare a standalone cleanup prototype

- **Status:** not started.
- **Outcome:** an auditable paired FASTQ transformation, separate from package defaults.
- **Actions:** use explicit Cutadapt end-specific settings on R2 only for supported layouts;
  retain corresponding R1 records unchanged; preserve IDs/order and sequence-quality lengths.
  Remove both mates together only when the post-trim cDNA fails the declared length rule.
  Run chemistry preflight on originals. Keep originals immutable and use a separate output root.
  Do not invoke generic paired-end poly-A defaults on barcode reads. Fail on unsupported/custom
  geometry, duplicate fragment IDs without an established namespace, malformed FASTQ, or
  missing tools. Reuse paired-ID checks and lineage writers where their contracts match.
- **Dependencies:** VR-01; future implementation authorization and pinned optional runtime.
- **Affected areas:** smallest analysis prototype under the existing project workflow, existing
  fragment-lineage convention and truth fixtures. No production CLI or dependency change yet.
- **Validation:** meaningful fixtures cover body+tail, reverse orientation, partial TSO,
  reagent-only and too-short reads, genuine A/T-rich bodies, internal homopolymers, mismatched
  mates, truncation, and duplicate IDs. Compare original and retained R1 record hashes. Verify
  `input pairs = retained pairs + removed pairs`, exact mate synchronization, and trimming bounds.
- **Evidence:** pending commands, reports, output hashes, per-fragment retained/removed decisions
  and reasons, plus trim coordinates/lengths in a companion audit keyed by library + fragment ID.

### VR-03: evaluate cleanup and decide integration

- **Status:** not started.
- **Outcome:** bounded evidence for whether the optional recipe helps, with failures retained.
- **Actions:** add end/tail/orientation truth cases to existing planting tools, including genuine
  complex viral bodies adjoining poly-A/poly-T, legitimate viral homopolymers, TSO junctions,
  and host/reagent-only controls. Run all declared arms on identical original pair budgets;
  compare true-molecule retention, false assignments, host-virus ambiguity, short-read losses,
  cell support, runtime, and resource use. Track loss through preprocessing and host subtraction
  separately. Use library-namespaced molecule keys and an outcome-independent cell anchor.
- **Dependencies:** VR-01/02 and actual project execution eligibility. This document opens no
  scored execution gate. An exploratory result stays exploratory and does not qualify a default.
- **Affected areas:** existing truth generator/scorers and receipts. Real public positives and
  artefact-heavy libraries supply contextual checks; presumed negatives are not absolute truth.
- **Validation:** meet the predeclared target/chemistry margins; show per-target results rather
  than only an average. Internal-fragment controls alone do not pass tail-trimming validation.
  Include clean no-tail reads, legitimate low-complexity viral reads, and independent molecules.
- **Evidence:** pending per-arm comparison, loss ledger, uncertainty/replicate counts, and verdict.
- **Integration decision:** only after a useful pilot, add one opt-in preprocessing stage before
  host alignment and counting. Scope all callers: CLI parsers, defaults, `RunConfig`, menu,
  `Snakefile`, host-filter input/dependency routing, `kb_r1`/`kb_r2`, exact evidence/host-capture
  replay, fingerprint/resume handling, runtime locks, output docs, and ship-scope entries if needed.
  Give its outputs distinct paths; do not silently move existing filters.
  Focused tests must cover DAG order, transformed-input replay, changed recipe/input invalidating
  resume, old/default invocation compatibility, and failures leaving no completion sentinel.
  Refresh generated docs/governance only after the authorized final source changes.

### VR-04: assess optional transcript reconstruction

- **Status:** not started.
- **Outcome:** supported viral transcript fragments and junctions for selected samples/targets.
- **Actions:** use exact candidate extraction and lineage, preserve qualities where available,
  and create a separate coordinate-sorted STAR RNA alignment with appropriate strand tags.
  Include host and relevant viral competitors during candidate adjudication; do not treat
  viral-only mapping as uniqueness. Run StringTie with the matching viral annotation as guidance
  where available; retain novel candidate models separately. Keep existing evidence BAMs intact.
  Start per biological sample, namespace every constituent library, and separate strain/accession
  hypotheses. Avoid pooling unrelated donors. Preserve reads carrying distinct exon evidence;
  report independent molecule support without blindly collapsing all reads to one per CB/UMI.
- **Dependencies:** VR-01 plus a separately scoped runtime/execution decision. VR-03 is a priority
  decision point, not a claim that trimming must succeed before assembly is technically possible.
- **Affected areas:** future optional evidence-stage wrapper, existing extraction/competitor
  outputs, STAR/StringTie receipts, per-sample GTF and junction-support tables. No automatic
  reference/catalogue rewrite or programme-state changes.
- **Validation:** a fixture with a known spliced viral transcript and a no-splice control;
  host/sibling-ambiguous and PCR-duplicate cases; no cross-library CB/UMI merging. Check GTF
  contig names, coordinates, strand tags, supported junctions, and real-data browser alignments.
  Flag inadequate coverage instead of interpreting failed assembly as viral absence.
- **Evidence:** pending GTF, junction table, selection/mapping provenance and coverage report.
  Report candidate fragments. Three-prime/five-prime capture cannot establish complete transcripts
  or full viral genomes. Reconstruction from the same selected reads is not independent truth.
  Novel/divergent-virus discovery from broader raw reads is a separate future scope.

### VR-05: scope candidate host-virus junction reporting

- **Status:** not started; deferred until the cleanup/assembly decisions are recorded.
- **Outcome:** a narrow design for selected biologically relevant targets, not an integration claim.
- **Actions:** borrow Venus's combined host-virus chimeric-alignment idea; preserve original
  candidates including host-compatible reads, not only the survivors of host subtraction.
  Specify anchor length, uniqueness, repeat/microhomology checks, independent molecule support,
  strand and virus-specific terminal/splice context. Reuse browser/evidence artefacts.
- **Dependencies:** a separate scoped design decision; positive junction truth and reagent
  chimera controls. Do not copy Venus thresholds or biological classes universally.
- **Affected areas:** proposed optional junction table and browser tracks; avoid altering core
  detection calls or declaring DNA integration from RNA evidence.
- **Validation:** known host-virus junction versus TSO/template-switch, repeat/homology and
  PCR-duplicate controls. Record false junctions and inadequate evidence explicitly.
- **Evidence:** pending bounded design. Genomic integration needs independent supporting evidence.

## Source receipts and interpretation

Sources were read during this session; upstream revision pinning remains VR-01 work.

- [Venus paper](https://doi.org/10.1371/journal.pcbi.1010636): subtractive detection and targeted
  hybrid-genome chimeric alignment. Its preprocessing prose is not an exact executable recipe.
- [Venus integration source](https://github.com/aicb-ZhangLabs/Venus/blob/main/src/module-integration/module-integration.py):
  explicitly sets quality 5 and length 20, without an explicit poly-A option. Borrow the analysis
  concept, not these unqualified defaults or the wrapper implementation.
- [Trim Galore 0.6.7 source](https://github.com/FelixKrueger/TrimGalore/blob/0.6.7/trim_galore):
  default adapter trimming and specialized poly-A handling are separate paths.
- [Cutadapt guide](https://cutadapt.readthedocs.io/en/stable/guide.html#poly-a-poly-t-trimming):
  generic paired-end `--poly-a` acts on R1 poly-A and R2 poly-T; single-cell geometry needs an
  explicit cDNA-only recipe and a pinned version.
- [Viral-Track assembly source](https://github.com/PierreBSC/Viral-Track/blob/master/Viral_Track_transcript_assembly.R):
  extracts aligned viral reads, merges BAMs across input libraries, and runs StringTie with
  `-f 0.01 -g 5`. This is reference-based transcript reconstruction, not de novo genome assembly.
  Its cross-library pooling and settings are not adopted automatically.
- [StringTie manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual): sorted RNA
  alignments, spliced-read strand tags, annotation guidance and junction-support settings.
- [Viral-Track README](https://github.com/PierreBSC/Viral-Track#transcriptome-assembly): warns that
  conventional scRNA-seq coverage bias makes transcript assembly challenging.

## Final validation and promotion boundary

Before any package integration, assign each touched file to one owner and recover the then-current
checkout/dirty state. Run focused boundary/DAG tests and a real native-tool smoke check on that
implementation. Broader package/governance checks belong to the integrated code stage, not this
Markdown review. Existing historical receipts do not qualify future code or scientific outcomes.

Promotion requires recorded margins, adequate biological replication for the intended claim,
training/holdout separation, and the existing project acceptance process. Off remains the default
unless evidence and a later explicit decision support changing it. Assembly and junction outputs
remain optional evidence products with independent recipes and limitations.

## Decision and evidence log

- 2026-10-14: reviewed against HEAD `5e551c9` plus current source files and existing plan notes.
  Saved the revised proposal. No implementation, installation, outcome run, or submission occurred.
- Unresolved before VR-02/03: supported chemistry/strand recipes, exact pilot inputs and truth
  tiers, numerical acceptance margins, source/runtime pins, and execution eligibility.
- 2026-10-14 execution update: the authorized prototype lives in `analysis/viral_read_cleanup/`;
  tests in `tests/test_viral_read_cleanup_pilot.py`. Three first-1M-pair subsets were prepared.
  Output attempts are under `.tmp/viral-read-cleanup-20261014*`; originals and frozen DSR outputs
  were not changed. Configs record runtime/reference hashes and provisional technical margins.
- Cutadapt 5.2 and StringTie 2.2.3 were installed/downloaded into isolated `/tmp` locations.
  Native STAR/StringTie reconstruction recovered the known two-exon and single-exon fixture:
  `.tmp/viral-read-cleanup-20261014-v2/fixture/assembly/assembly_summary.json`.
- The v2 real-data attempt failed before quantification because sequence-search-derived trim
  coordinates were ambiguous in a repeat. The prototype now uses explicit TSO offsets and
  Cutadapt's per-read match coordinates. Nine focused tests pass after the fix; a regression
  using the reproduced real repeat and cleanup of seven lint findings are still required.
- The corrected v3 config and byte-identical copied subset receipts are frozen in
  `.tmp/viral-read-cleanup-20261014-v3/`; its 31-job dry-run passed with existing subsets reused
  through `--rerun-triggers mtime`. No v3 execution or real assembly has started.
  Chemistry receipts validate CB/UMI geometry but their `end` is null: the pilot's sense-tail
  orientation is a declared provisional recipe, not verified 3-prime/5-prime library metadata.
  Full upstream revision pins and biological qualification remain open. No default promotion.
