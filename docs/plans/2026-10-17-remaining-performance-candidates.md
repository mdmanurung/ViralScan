# Remaining performance candidates: implementation and qualification plan

Date: 2026-10-17. Project: ViralScan. Status: **planning complete; implementation not started**.

This document turns the 27 unimplemented survey candidates into bounded work cards.
It also records preservation requirements for the three completed candidates. A card
is ready to start at its first unchecked step; that does **not** mean its proposed
optimization has already passed an equivalence or benefit gate. Some cards may end
with a measured decision to retain the current implementation.

The current request authorizes drafting this plan. It does not authorize executing
its implementations, changing scientific methods or resource budgets, submitting
SLURM jobs, downloading data, deleting retained data, committing, pushing or publishing.
Future execution must use the authorization then in force. Existing retained inputs
are read-only. No held-out or protocol-frozen analysis is a performance fixture.

## Evidence, current source and review

Read these documents before executing any card:

- [Computational survey](../../analysis/computational_bottlenecks_2026-10-17/REPORT.md)
  and its linked scripts, source excerpts, profiles and measurement receipts.
- [First-wave plan](2026-10-17-performance-triage-byte-identical-outputs.md), including
  the complete artifact contract and failed sparse-reduction probes.
- [First-wave execution report](../../analysis/computational_bottlenecks_2026-10-17/execution/REPORT.md),
  [validation receipt](../../analysis/computational_bottlenecks_2026-10-17/execution/validation.json),
  [follow-up pilot](../../analysis/computational_bottlenecks_2026-10-17/execution/followup_assessment.json),
  and `tests/fixtures/performance_original.py`.
- Applicable repository instructions, `PLAN.md`, configuration, serializers and the
  actual callers/tests named by the selected card. Re-read affected bodies after drift.

At **2026-10-17 10:46:40 UTC**, HEAD was
`5e551c96831c0c152da9fa450af6a75b1266e610`, branch `codex/drop-py39`.
All **76** live package/resource files in the captured combined manifest matched
their recorded SHA-256. HEAD alone is insufficient: B05/B07/B20 remain uncommitted
alongside unrelated analysis/tracker changes. Do not reset or stash this checkout.

| Evidence anchor | SHA-256 at planning review |
|---|---|
| `execution/variants/combined/variant.json` | `fa3d1bda07c9bd86489b04e88525a788947a6605e0dd2da2e20f984db92ef11f` |
| `execution/source_baseline.json` | `41005b320e89d16e5a4e88fcc46bd1cd590106a7467bef1fff635cbcaffc3f02` |
| `execution/validation.json` | `54f8f23e4e2a533562646570145c1b669c113cc4261f226b8b8498ad84a27c19` |
| `execution/followup_assessment.json` | `f65aee3ca9f5cac666e442eadbfd21383cbb0591360dce1c51e218eeb70322f7` |

Paths in this table are relative to
`analysis/computational_bottlenecks_2026-10-17/`.
The sorted manifest of `relative_path + TAB + sha256 + NEWLINE` for the 76 live
files has SHA-256
`26de210659189473de4b210090de3ff8641b1037a5b5bf4c0c2fb403eade93a8`.
The combined manifest contains each individual file hash. Recheck it before edits.
The first-wave validation receipt also binds the earlier `PLAN.md`; the additive
link to this new plan is subsequent documentation drift, not new test evidence.
Do not overwrite that receipt to conceal this distinction.

Actual bodies reviewed for this plan include allocation/EM/assignment emission,
detection reducers/denominators, UMAP entry, programme call/pruning, FASTQ filtering,
replay/capture, alignment/QC/BLAST, verdict reduction, anellovirus accumulators,
host-response feature construction/fits, reference parsing/gates/native invocation,
cache validation, chemistry matching and provenance hashing. The survey's original
recommendations are historical: direct sparse sums and the broad CSC experiment are
not approved replacements merely because their prototype timings were favorable.

| Review finding | Resolution |
|---|---|
| Three optimizations are already implemented | Keep them; baseline this wave from the captured combined source, while preserving the pre-first-wave reference independently. |
| B08's early timing looked favorable | Use the later pilot: paired reporting median worsened; require renewed post-wave evidence before an edit. |
| A fast sparse sum changes numerical bits | B04 retains individual dense column reductions; B06 tests bounded dense row blocks with identical reduction width/layout. |
| Disabled UMAP validates malformed inputs | B10 begins with a failure-contract feasibility experiment; a completion marker alone cannot justify skipping checks. |
| Native replay/thread behavior can be nondeterministic | Require original A/A raw-artifact measurements first; resource changes remain a separate authorization gate. |
| Some apparent skips alter logs, warnings or partial files | Capture these interfaces and retain original error ordering; use fallback or defer when they cannot be preserved. |
| Memory gains and component speedups are different outcomes | Measure separately; no multiplication into a claimed whole-pipeline gain. |

## Order, dependencies and progress

Priority reflects likely benefit, edit effort and exactness risk, not an established
production runtime ranking. Effort estimates are engineering time after fixtures and
tools exist; qualification of native paths may take longer. No numerical speedup is
promised for an unimplemented card.

| Wave | Cards | First executable outcome | Rough effort per card |
|---|---|---|---|
| 0 | G00–G02 | Freeze current baseline, variability contract and stage profile | 0.5–1 day |
| 1 | B24, B23, B04, B06 | Local counter/access/pruning or bounded-memory prototype | 0.5–2 days |
| 1, conditional | B08, B10 | Renewed benefit or validation-feasibility decision | 2–8 hours before any edit |
| 2 | B25, B27, B18, B29, B17 | Bounded text/state or shared decoding | 1–3 days |
| 3 | B13, B14, B15, B16, B19 | Filtering/replay/native-I/O prototype | 2–5 days; B15 can be longer |
| 4 | B09, B01, B02, B03, B21, B22 | Reader/state/model workspace work | 2–5 days; numeric refactors can be longer |
| 5 | B26, B28, B12 | Integrity/cache or R-stage feasibility decision | 1–3 days |
| Separate gate | B11, B30 | Resource design and fixed-worker qualification | 1–3 days after authorization |

Do not implement an entire wave at once. Qualify one isolated card against the same
frozen wave baseline before combining accepted changes. G02 may reorder independent
cards if measured stage shares justify it; record why. B17/B18/B29 are independent
consumers of alignment data and should not be folded into one broad rewrite.

- [x] P00 — Read source, reconcile all 30 candidates, review and save this tracker.
- [ ] G00 — Capture the current optimized baseline. **Status: not started.**
- [ ] G01 — Freeze new originals and original A/A behavior. **Status: not started.**
- [ ] G02 — Measure relevant stage shares. **Status: not started.**
- [ ] C00 — Qualify combined accepted changes. **Status: not started.**

Each B-card below has three separately checkable steps: establish feasibility,
implement the smallest qualifying change, qualify it. Every unchecked B-step has
status **not started**. Do not mark implementation complete when only a prototype
or assessment exists. If an assessment rejects implementation, leave its edit step
unchecked and record the disposition as **deferred** with the failed gate.

### G00 — Capture the baseline without disturbing shared work

Outcome: an independently executable, immutable copy of the current optimized source.
Dependency: P00; update `PLAN.md` to show the selected card before implementation.

Capture `git status --short`, branch/HEAD, all runtime/resource/schema/template/test
hashes, relevant untracked package files, interpreter/dependency/native binary hashes,
native defaults, locale, timezone, RNG seeds, configuration, effective thread counts,
CPU model/affinity, host load and input digests. Confirm imported module paths in
every run; an installed stale package is not the working-tree candidate.

Use a new unused directory under the existing survey packet, for example
`analysis/computational_bottlenecks_2026-10-17/followups/<UTC>-<card>/`. Create it only
when executing. Snapshot baseline and isolated candidate separately; no destructive
Git operations or hardlinks to files that will be edited. A frozen-input copy or
read-only file handle is required for reuse based on immutability. Size/mtime alone
does not establish that an external input cannot change.

Recheck the 76-file manifest and re-audit affected claims if it differs. Preserve the
historical `execution/original/`, `runs/original-a/`, portable oracles and receipts.
Historical verifiers intentionally constrain edits to the first three functions;
they must not be weakened to accommodate a new candidate. Create a new wave verifier
with its own explicit scope, while verifying old saved evidence independently.
Evidence: `baseline.json`, source snapshot, environment/input manifests, command log.

### G01 — Establish exactness before production edits

Outcome: frozen outputs produced by G00's baseline and a tested strict comparator.
Dependency: G00. Before edits run baseline A/A in isolated processes on identical
inputs, paths and test-only rendering settings. Capture the candidate-specific
fixture matrix below, plus applicable retained development data when bounded.
Do not claim that changing six method labels on fixed layers is six allocator refits.

Compare relative file sets, lengths and actual raw bytes, recording the first
difference as well as SHA-256. Capture exact types/dtypes/shapes/numeric bits,
dictionary/table/sparse-buffer order, warnings, exception class/message/cause,
exit status, stdout/stderr/log events, completion markers, partial outputs and
preservation of pre-existing outputs. Hash inputs before/after; compare RNG state
where consumed. Validate the comparator with intentional one-byte, missing-file,
type/dtype and failure-status mutations. A failed mutation test blocks qualification.

Default fixtures: integer/fractional/empty/negative-result/boundary values; dense,
CSR and CSC where supported; explicit zeros, unsorted indices, duplicate coordinates,
missing/aliased layers, unknown/duplicate genes, repeated names and existing invalid
inputs. Preserve unsupported-input behavior instead of silently expanding support.

Do not use `allclose`, rounding, sorting, decompression, metadata removal or rewritten
serialization to pass a raw-byte gate. Decoded comparisons are diagnostic only.
Stable compressed/BAM/H5AD/plot bytes require raw comparison. If original A/A differs,
investigate and document the exact variability; do not turn nondeterminism into a
blanket artifact exclusion or assume semantic equality proves raw identity.

Keep truthful runtime/source/tool provenance. Predeclare precise unavoidable
operational fields from observed A/A evidence, with full raw hashes still reported
as different. The first-wave four 64-character provenance spans do not authorize
new exclusions in native paths. No whole-directory identity claim if any bytes differ.
New unexplained variability blocks acceptance for that artifact until resolved.
Evidence: `original-a/`, `original-b/`, immutable golden manifest, artifact contract,
typed/RNG/error records and comparator self-check receipt.

### G02 — Measure benefit on the affected workload

Outcome: a stage profile and fair repeated baseline/candidate comparison.
Dependency: G00–G01; profiling may precede the candidate edit. Use disposable local
fixtures and already retained approved development inputs; no new job/data request
is implied. Do not flush global filesystem caches on this shared host.

Keep inputs, software/native versions, scientific configuration, iteration budgets,
thread counts and numerical environment fixed. Use a warm-up for both variants and
at least seven paired trials with alternating A/B order; repeat in a second timing
session. Separate profiling, timing and memory instrumentation. Report every trial,
median and spread, warm/cold-cache scope, child wall/CPU time and bytes read/written.
Whole-stage measurements must include cache preparation, conversion and serialization.

Default benefit gate: at least 10% median affected-stage wall-time reduction in both
sessions, or at least 25% isolated peak-memory reduction with no more than 5% stage
wall-time regression. These are proposed engineering acceptance thresholds, not
scientific parameters. A departure needs an explicit recorded rationale before
acceptance; it cannot excuse exactness failure. A noisy result stays inconclusive.
Measure optional-feature applicability and report component, stage and whole-run
results separately; if whole-run execution is unavailable, say so.

Measure sampled process RSS across the function/stage and native child tree in fresh
processes, alongside traced allocations and fixture-construction peaks. Report
sampling interval and transient-peak limitations. Lifetime `ru_maxrss` already reached
during setup is not a measurement of additional function RSS. Include steady memory,
transients, cache/spool disk usage and worst-case high-diversity inputs.
Evidence: `profile/`, `timing.json`, `memory.json`, native command/provenance receipts.

## Candidate work cards

All code paths below are relative to `src/viralscan/`; test paths are relative to
`tests/`. Every card depends on G00–G02 and the shared acceptance/rollback rules.
Additional dependencies and acceptance details appear beside the relevant card.
Evidence for each card goes to its new packet as `feasibility.md`, `diff.patch`,
`exact.json`, test log/JUnit, timing/memory receipts and `decision.md`; native cards
also retain real native command, binary hash, exit and artifact receipts. An absent
receipt is pending evidence, never a pass.

### B24 — Reuse reference low-complexity counters

**Priority:** first build-focused candidate; moderate effort, low arithmetic risk.
**Reach:** `scripts/build_reference.py`: `low_complexity_kmer_counts`,
`low_complexity_report`, `validate_reference_records`, `_enforce_low_complexity_gate`,
both builders and `prepare_reference_inputs`; reference FASTA/GTF/manifests/logs.
The classifier currently runs three times per record in the gate. Near 3× is an
upper opportunity for repeated classifier work, not a measured gate/build speedup.

- [ ] B24.1 — Freeze valid/invalid gate calls and count invocations. Include empty
  FASTA/record, duplicate IDs and sequences, duplicate IDs with distinct sequences,
  lowercase/N/IUPAC, short records, all classifier categories, exact fraction and
  homopolymer limits, tied worst records and more than ten tied offenders. Capture
  the current report's last-value/first-key behavior for duplicate identifiers.
- [ ] B24.2 — Compute per-record integer counters once within one gate invocation;
  derive the existing report, pure total and threshold verdict from them. Reuse
  structural/verdict helpers privately while retaining public signatures. Keep
  initial scan logging before structural failure, structural `ValueError` outside
  the low-complexity `RuntimeError` wrapper, stable offender sorts and causes. No
  process-global cache, classifier rewrite, mask change or threshold change.
- [ ] B24.3 — Run `test_build_reference.py`, `test_reference_pipeline_aud15.py`,
  `test_reference_manifest_contract.py`, `test_bundled_panel_mask.py`; add strict
  diagnostic/order oracle checks. Compare prepared FASTA/GTF/report bytes and all
  error/log records. Time whole gates on clean random and N-rich/masked panels at
  several lengths, then offline reference preparation. Native build claims require
  real identical-tool builds; retain fallback if shared helper behavior drifts.

**Stop:** a repeated scan cannot be removed while preserving duplicate/error/log
semantics, or preparation-inclusive benefit fails G02. Treat targeted masking's
per-window classifier optimization as a later separate change, not this card.

### B23 — Avoid programme dictionaries that the writer discards

**Priority:** early optional-feature candidate; benefit unmeasured.
**Reach:** `gene_programs.py:call_cell_programme`, `_nonzero_columns`, `_decide`;
`scripts/gene_programs.py:_marker_matrix`, `run_one`, `main`; summary/log serializers.

- [ ] B23.1 — Profile marker-negative cohorts and freeze the exact prune predicate:
  `state == indeterminate` and zero productive/latent breadth on both layers. A
  zero matrix alone does not establish that a call is prunable. Include selected-only
  evidence, unresolved markers, absent models, absent layers, overlapping groups,
  min-UMI boundaries, non-assessable absence and latent/reactivating flags.
- [ ] B23.2 — Keep the public helper's all-cell return contract. In the script path,
  first try bounded row blocks with the unchanged helper and immediate pruning;
  accumulate state counts for **all** calls, including discarded rows, in original
  first-seen order. Do not log per block. Reconstruct the existing pandas value-count
  tie order for the final log or retain that path. If only memory improves, qualify
  that result before attempting an index-skipping helper with a broader API change.
- [ ] B23.3 — Run `test_gene_programs.py`, `test_gene_programs_audit4.py`,
  `test_gene_programs_unresolved.py`, `test_program_evidence_boundary.py`,
  `test_gene_programs_capturability.py`; compare cells/summary/JSON/manifest bytes,
  marker-resolution metadata and logs. Benchmark sparse-negative and dense-positive
  called-cell cohorts, including helper work and serialization; measure peak RSS.

**Stop:** summaries/log tie order or helper validation timing changes; keep all-cell
calls rather than dropping scientifically meaningful zero-evidence states.

### B04 — Cache a narrow local detection view, keep dense column reductions

**Priority:** early measured hotspot; no approved direct sparse/batched sum.
**Reach:** `scripts/detection.py:_gene_counts_from_matrix`, `detect_genes`,
`preprocessing`; `found_genes.txt`, group detection and downstream reports.

- [ ] B04.1 — Profile post-first-wave detection, recording requested viral-column
  fraction. Freeze duplicate-name validation, repeated accession traversal, Python
  group-sum order, scalar result dtypes and nonfinite error lists. Reuse adversarial
  sparse-sum counterexamples from the first-wave probe as mandatory regression cases.
- [ ] B04.2 — For reusable sequences and a proven supported layout, select requested
  columns once into a bounded **local** CSR view; map each original index to it.
  Continue the same individual `toarray().sum()` with its original shape/dtype and
  order. Preserve input matrices and the duplicate-name check before access. Fall
  back for generators, unsupported/noncanonical representations, all-column cases
  and over-budget views, as needed. Do not import B05's cache assumptions untested.
- [ ] B04.3 — Run `test_detection.py`, `test_var_names_safety.py`,
  `test_infected_rule_and_gate.py`, `test_detection_audit_s2.py` and stage regressions.
  Cover dense/CSR/CSC, all numeric fixture classes, repeated/missing genes and group
  threshold boundaries. Compare typed totals plus actual detection-stage bytes;
  benchmark complete preprocessing with cache construction and high selected fractions.

**Stop:** any total bit/dtype, group decision or error differs, or cache memory defeats
the measured benefit. Never accept a sparse reduction because `allclose` passes.

### B06 — Bound dense row-reduction memory

**Priority:** early memory candidate; direct sparse sum has failed exactness probes.
**Reach:** `utils.py:virus_molecules_per_cell` and its detection/enrichment/UMAP/plot
callers; separately `compute_stats` accession-breadth and ambiguity-total operations.

- [ ] B06.1 — Characterize existing dense layout/strides and row-reduction dispatch
  for each supported sparse representation/dtype and selected gene count. Freeze
  cancellation/large-positive mixtures, signed zeros, float32-to-float64 conversion,
  repeated/missing genes, zero rows/columns and infected-threshold neighbors.
- [ ] B06.2 — First optimize only the shared per-cell helper: dense **row blocks**,
  with all selected columns in the same order, the same float64 conversion, row
  width and array layout, followed by the same `sum(axis=1)`. Never split a row's
  reduction across column blocks or combine partial sums. Choose a private byte
  budget, guard oversized single rows, and retain the old path for unproven layouts.
  Leave flattened scalar ambiguity sums unchanged; they are a different reduction.
- [ ] B06.3 — Run detection, infected-rule, enrichment/UMAP and performance-exactness
  tests; exercise block sizes 1, 2, non-divisors and larger than row count across
  widths around NumPy reduction boundaries. Require exact arrays, warnings and all
  caller artifacts. Benchmark wide-virus and narrow-virus cases; show reduced
  isolated RSS including temporary conversion buffers. Treat breadth-only streaming
  as a separately isolated subchange after helper acceptance.

**Stop:** blocking changes row reduction bits/layout-sensitive behavior, or warning
category/message/multiplicity cannot be preserved for allowed nonfinite inputs.
Nonfinite inputs can use the original path only if the guard itself preserves errors.

### B08 — Reuse the comparable-cell denominator lazily

**Priority:** conditional; later pilot did not improve paired reporting.
**Reach:** `scripts/detection.py:compute_stats`, `_comparable_called_cells`,
`_in_barcode_set`; host-set and combined-reference denominators.

- [ ] B08.1 — Reprofile after B07/B05 on many-virus combined references and external
  host sets. Repeat the stats-only **and** paired reporting measurement in
  `execution/assess_followup.py` using a new packet. The saved 20→1 helper count and
  stats-only gain do not qualify the paired 1.7338→1.9760 s result as a speedup.
- [ ] B08.2 — Only after renewed benefit, initialize a local sentinel and calculate
  the identical denominator at the first original use inside a valid-virus iteration;
  reuse thereafter. Distinguish uncomputed from computed `None`. Do not hoist above
  earlier failures or evaluate it for zero viruses/all-invalid gene lists.
- [ ] B08.3 — Run detection/sensitivity/ambient tests and stage byte comparisons.
  Cover missing/empty host dictionaries, called/comparable `None` and empty sets,
  no host genes/no called cells, malformed masks, candidate-negative viruses,
  threshold boundaries and original caller warning order; assert input immutability.
  Include conversion/preparation/serialization in the benefit measurement.

**Stop:** renewed paired reporting gain fails G02, any lazy branch/error changes, or
custom/stateful supported inputs invalidate reuse. A smaller helper count alone is
not an acceptance criterion.

### B10 — Determine whether disabled UMAP can avoid eager loading exactly

**Priority:** feasibility first; old missing/corrupt/schema failures are mandatory.
**Reach:** `scripts/umap.py:main`, `run`; `Snakefile` UMAP rule; H5AD and found-gene
parsing, run identity, stdout and `log/umap.done`.

- [ ] B10.1 — Freeze enabled/disabled cases with missing/invalid/truncated H5AD,
  malformed X/layers/obs/var/uns, wrong schema, unreadable/malformed found-gene text,
  empty files, missing identity, unavailable Plotly and two simultaneous faults.
  Establish exactly which failures eager `read_h5ad` exposes and in what order.
- [ ] B10.2 — Prototype a reader-only bounded/lazy validation path for disabled
  plots, retaining every observed validation/error and found-gene conversion. Inspect
  backed AnnData behavior against actual installed code; `backed='r'` or reading
  `uns` alone is not automatically equivalent. Always close handles. Retain eager
  loading for unproven encodings. A prior detection marker does not prove present
  bytes are intact; no blind early return or trust-by-mtime shortcut.
- [ ] B10.3 — Run `test_umap.py`, `test_snakemake_script_modules.py` and relevant
  `test_snakefile_dag.py` cases. Require identical done-marker content/presence,
  stdout, errors and untouched input hashes. Measure complete disabled rule wall/RSS
  on retained H5AD encodings; enabled plotting remains on the original path.

**Stop:** equivalent validation still needs the entire read or a malformed case is
silently accepted. Record **no compliant fast path** if so. Relaxing disabled-mode
validation requires a separate explicit behavior decision; this plan does not grant it.

### B25 — Stream reference text without changing byte or failure contracts

**Priority:** memory opportunity, start with one adapter.
**Reach:** `menu.py:_materialize_reference_input`; `scripts/build_reference.py`:
`_fasta_records`, `prepare_reference_inputs`, header/GTF readers; NCBI merge is B28.

- [ ] B25.1 — Measure overlapping raw-byte/text/list/record allocations and choose
  the dominant adapter. Freeze CR/LF/CRLF and trailing-newline runs, blank inputs,
  descriptions, lowercase, Unicode text errors, malformed GTF and duplicate IDs.
  Note that `read_text().splitlines()` recognizes separators beyond a newline iterator.
- [ ] B25.2 — First replace whole-byte concatenation with bounded copying that
  buffers/discards only trailing CR/LF, emitting exactly one LF after each nonempty
  stripped input. Then separately assess streaming headers/GTF with original text
  splitting rules. Preserve uppercase/wrapping/record order and validation-before-
  write boundaries; spool before publication where required. Keep public list-return
  APIs until their callers have a separately qualified streaming implementation.
- [ ] B25.3 — Run build-reference, reference-manifest, reference-pipeline-aud15 and
  CLI tests; compare raw prepared/merged FASTA/GTF and manifest bytes, exception and
  partial-output behavior. Benchmark large synthetic host text and small panels,
  reporting isolated peak RSS, spool footprint and total preparation time.

**Stop:** newline semantics or failure side effects differ, or a sequence-storage
refactor is needed to achieve the proposed gain. Do not bundle parser redesign here.

### B27 — Reuse bounded whitelist matching preparation

**Priority:** low, evidence first; ordinary read sampling is already bounded.
**Reach:** `chemistry.py:sample_r1`, `match_rates`, `detect`, `classify`, `resolve`;
chemistry preflight callers. Keep the separate full-list preflight API unchanged.

- [ ] B27.1 — Profile many samples and candidate lists, counting repeated barcode
  lengths and preparation cost versus actual whitelist I/O. Freeze duplicate entries,
  `-1` suffix stripping, blank/short first lines, gzipped lists, empty/short reads,
  tied rates, label order and unreadable paths.
- [ ] B27.2 — Try per-invocation reuse of sampled barcode Counters by inferred
  barcode length, retaining separate original per-label whitelist scans and ordered
  integer hit accumulation. Bound cache memory; build lazily at the original access
  point. Do not deduplicate whitelist lines: duplicates currently contribute more
  than once. Do not retain million-entry Python whitelist sets or cache across calls.
- [ ] B27.3 — Run `test_chemistry.py`, `test_chemistry_check.py`,
  `test_whitelist_preflight.py`; compare exact rate floats, detections, warnings,
  requested-chemistry errors and output text. Measure complete preflight, worst-case
  distinct lengths and memory; if I/O dominates, record a negative result.

**Stop:** benefit is below G02 or cache memory exceeds the agreed baseline budget.
Already cached native-tool discovery is not new work to optimize.

### B18 — Reduce per-read verdict state while preserving first-seen order

**Priority:** bounded-memory evidence work; no sort-by-name shortcut.
**Reach:** `molecule_verdict.py:best_scores`, `molecule_rows`,
`write_molecule_verdicts`, `write_from_bam`; evidence/dossier consumers.

- [ ] B18.1 — Profile read versus molecule cardinality. Freeze order based on the
  first **usable scored alignment**, not first raw QNAME occurrence. Include secondary
  hits, supplementary/unmapped-only reads, malformed names/AS tags, score ties,
  interleaved molecules and counted-read mappings.
- [ ] B18.2 — Prototype a fused reducer: retain a seen-usable-QNAME structure for
  unique-read counts and per-molecule best-side state; update maxima in stream order.
  Preserve molecule insertion order and the public `best_scores` API. The old code
  parses the entire SAM before parsing molecule names; adversarial mixed faults must
  retain that exception precedence. If not feasible, keep a two-phase disk-backed
  ordered state representation or reject the refactor. Stream output only after
  matching the original parse-before-write boundary; keep gzip writer/header policy.
- [ ] B18.3 — Run `test_molecule_verdict.py`, `test_dsr02_verdicts.py` and evidence
  tests; compare verdict gzip raw bytes, summary, errors and partial files. Benchmark
  high reads-per-molecule and all-unique reads, including disk/state overhead.

**Stop:** first-seen order, unique-read count, validation precedence or compression
bytes change. This may reduce per-read dictionaries without making memory constant.

### B29 — Decode anellovirus BAM once for both existing statistics

**Priority:** optional-stage I/O candidate; independent of B17.
**Reach:** `scripts/anello_align.py:main`, `_alignments`;
`anello_align.py:accession_metrics`, `virus_molecules`, `virus_summary` and writers.

- [ ] B29.1 — Capture both consumers on the same decoded stream, including duplicate
  QNAMEs per accession, same-virus multi-accession reads, cross-virus ambiguity,
  missing/invalid CB/UB, circular wraps, absent lengths and empty/malformed BAM.
- [ ] B29.2 — Share one decode through explicit accumulator update/finalize paths;
  retain public helper behavior. Update each consumer in the original record order,
  keep `weighted += 1/NH` before accession read dedup, and preserve last-valid CB/UB
  assignment in the molecule consumer. Do not use `itertools.tee` with a whole-pass
  lag, which merely buffers the entire stream. Preserve first-pass versus second-
  pass error precedence, or fall back when failures cannot be matched.
- [ ] B29.3 — Run `test_anello_align.py`, `test_anello_detection_merge.py` and, when
  tools exist, `integration/test_anello_align_star.py` and anellovirus-chain checks.
  Compare accession/virus TSVs, metadata, merged detection bytes and completion
  failures. Measure complete summarization with real `samtools` on a fixed retained
  or tiny native BAM; quantify decode count and peak state memory separately.

**Stop:** numerical accumulation, unique-virus molecule semantics, output order or
original failure behavior changes. No STAR flags, filtering or index changes.

### B17 — Share evidence QC inputs with unchanged reducers

**Priority:** repeated BAM work; moderate implementation and native qualification.
**Reach:** `evidence.py:alignment_qc_table`, `per_cell_alignment_qc`,
`coverage_depth_points`, `plot_coverage_comparison`; `scripts/evidence_run.py` raw/
deduplicated QC path, CSV and PNG writers.

- [ ] B17.1 — Instrument exact `samtools` command counts and record order. Freeze
  primary/secondary policies, depth flags, missing NM/AS/CB tags, empty BAM, zero
  reference length, odd/even depth medians and raw-versus-deduplicated outputs.
- [ ] B17.2 — First decode SAM once per distinct BAM and feed the existing text
  reducers; reuse identical-policy depth text between QC and plotting. Scope reuse
  to one immutable BAM at the caller; do not combine raw/dedup state. Keep native
  `samtools coverage` separate where its policy differs. Bound shared buffers/spool;
  do not retain both BAMs' full SAM simultaneously. Streaming accumulator extraction
  is a later isolated change if shared text alone qualifies.
- [ ] B17.3 — Run evidence, molecule-representative and read-start-profile tests plus
  real tiny evidence-chain integration. Compare all QC/coverage TSVs, sampled arrays,
  PNG bytes under identical renderer settings, command errors and output presence.
  Time complete evidence QC including plotting; measure process/child RSS and I/O.

**Stop:** parsing order, float sums/entropy, depth inclusion or diagnostics change.
If an eliminated second command has distinct observable failure behavior, define
and test the same failure surface before accepting reuse.

### B13 — Optimize serial read classification before considering workers

**Priority:** optional artefact path; serial Python/gzip cost is measured.
**Reach:** `scripts/read_filter.py:classify`, `filter_pairs`, `trim_tso`;
`anello_align.py` reagent/body helpers; both FASTQs, lineage and audit writers.

- [ ] B13.1 — Reprofile classifier versus compression on diverse, reagent-rich and
  duplicate-rich reads. Edge tables and pigeonhole rejection already exist. Capture
  first reason, both-orientation rescue, CB/UMI span, TSO trim boundaries, quality
  preservation, mate mismatch/truncation and the exact partially written prefix.
- [ ] B13.2 — First try a bounded invocation-local cache of successful classifier
  results for repeated `(R1 CB/UMI span, R2)` pairs, with serial validation/writing
  unchanged. Include lookup cost and bounded key memory. If ineffective, select one
  measured pure classifier loop for an exact serial optimization; do not add a new
  compiled dependency. Keep gzip backend/level/write sequence and mode defaults.
  Ordered worker chunks are a separate resource change requiring authorization.
- [ ] B13.3 — Run `test_read_filter.py`, `test_anello_align.py`, host-filter/DAG
  interaction tests and tiny real count integration. Compare raw compressed mates,
  lineage, audit, failures and partial files with test-only fixed gzip clocks. Measure
  complete filter wall/RSS on low-repetition as well as cache-friendly fixtures.

**Stop:** high-diversity regressions exceed G02 or keys consume excessive RAM. Never
replace artefact mode with the faster but analytically different `tso-trim` mode.

### B14 — Share host-filter validation and lineage scans safely

**Priority:** optional full-read I/O; native STAR cost not yet isolated.
**Reach:** `scripts/host_filter.py:_starsolo_filter`, `_gzip_file`,
`validate_paired_fastq_ids`, `_write_filter_audit`, paired-ID iterator and callers.

- [ ] B14.1 — Profile two raw and two retained Python pair scans plus compression.
  Freeze exact canonical-ID, duplicate-ID, missing/extra-retained, out-of-order,
  mate-count, late truncation and STAR-more-pairs behavior. Record audit publication
  timing and partial outputs for each failure.
- [ ] B14.2 — First share retained mate validation/counts/ID collection in one
  streaming pass, passing privately collected state into the audit helper while
  keeping its direct-call validation behavior. Preserve aggregate audit publication
  before later lineage failures. Do not assume STAR output preserves input order.
  Original raw prevalidation stays before STAR; original raw lineage scan remains
  unless a separately qualified ordered spool can replace it. Compression fusion
  is another isolated experiment, with unchanged gzip byte/write policy.
- [ ] B14.3 — Run `test_host_filter.py`, `test_read_filter.py`, DAG tests and
  `integration/test_host_filter_star.py` when available. Require unchanged compressed
  mates, lineage/audit bytes and failure artifacts. Measure whole Python filter and
  actual STAR-inclusive stage separately, including retained-ID RAM and spool I/O.

**Stop:** skipped validation or changed failure/output order; saving a scan does not
justify changing exact fragment-level filtering or STAR arguments.

### B15 — Share one exact replay across multiple targets

**Priority:** potentially large multi-target gain; high effort and strict native gate.
**Reach:** `reads.py:extract_virus_reads`, `scripts/evidence_run.py:run_evidence`,
`evidence.py:_replay_quant_bus`, capture/replay/extraction helpers;
`scripts/auto_evidence.py` and `scripts/host_capture.py` entry paths.

- [ ] B15.1 — Measure one/two/many-target replay shares with real tiny tools. Run
  original A/A first: multithreaded kallisto can assign different EC IDs. Inventory
  all target workdir intermediates/manifests, slug collision handling, extraction
  order and target-1-success/target-2-failure artifacts before designing reuse.
- [ ] B15.2 — First share quantification only within one multi-target invocation;
  retain target-specific capture, EC maps and output paths. Bind reuse to FASTQ/index/
  transcript/whitelist bytes, chemistry geometry, strand, correction/sort policy,
  native binaries/versions and effective workers. Always use the replay's own EC map.
  Materialize required per-target intermediates without changing their contract.
  Batch FASTQ dispatch is a second isolated subchange, preserving per-target row,
  CB/UMI/read-number order and all exact lineage. Keep single-target fallback.
- [ ] B15.3 — Run reads, exact-lineage, evidence, auto-evidence, host-capture and
  strand-consistency tests; run `integration/test_exact_lineage.py`,
  `integration/test_strand_replay_kallisto.py`, evidence-chain and truth-panel checks
  with installed tools. Include gzipped/custom/no whitelist, stranded chemistries,
  overlapping/empty targets, cache mismatch/corruption and interrupted replay.
  Require all per-target raw artifacts, exits and provenance contracts to qualify.
  Benchmark 1/2/5 targets with unchanged workers, including staging/hash costs.

**Stop:** baseline EC/native bytes are not repeatable under the required contract,
sharing changes observed intermediates, or cross-target failure order cannot match.
No cross-run cache or pooled competitive reference in the first implementation.

### B16 — Bound the minimap2-to-samtools SAM handoff

**Priority:** native memory candidate; piping can change headers and failure order.
**Reach:** `evidence.py:_run`, `align_reads_to_viral`; `run_evidence` competitive
alignment, BAM/index and downstream evidence/provenance artifacts.

- [ ] B16.1 — Establish real original A/A BAM/BAI bytes and command logs with fixed
  tools/paths/workers. Freeze minimap2 failure before sort, sort failure, index
  failure, absent tools, stderr content, old-output preservation and partial outputs.
- [ ] B16.2 — First try a sequential disk-backed SAM handoff: minimap2 must finish
  successfully before invoking sort, with the same argv/stdin interface. Preserve
  full stderr and RuntimeError formatting. This bounds Python SAM memory while
  retaining failure precedence. A concurrent pipe is a later experiment requiring
  explicit two-child waiting, deadlock-safe stderr handling and matched failure
  precedence; do not introduce it merely to claim streaming.
- [ ] B16.3 — Run evidence/subcommand tests, real evidence-chain/competitor integration
  and fault-injection process tests. Compare **raw** BAM/BAI/header bytes and all QC
  artifacts; no removing `@PG` records or comparing decoded SAM as a pass. Include
  spool write failure, full stderr and producer/consumer failure. Measure full
  align/sort/index wall, child-tree RSS and temporary disk volume.

**Stop:** command/header/partial-file bytes or exception behavior differ. A memory
gain may qualify under G02 even without a speedup, but disk cost must be reported.

### B19 — Reuse an immutable BLAST database without changing search policy

**Priority:** optional evidence/build reuse; benefit unmeasured.
**Reach:** `evidence.py:blast_identity`, `competitive_blast_identity`,
`write_competitor_fasta`, sampling/parsing helpers; evidence manifests/raw-hit files.

- [ ] B19.1 — Profile database construction separately from sampling/search. Freeze
  single/multi-class paths, empty samples, complete subject counts, identical-score
  subjects in different order, tie delta, low complexity and build/search failures.
  Determine whether generated DB files are part of the required raw file set.
- [ ] B19.2 — Prototype invocation-local reuse for **identical ordered FASTA bytes**,
  binding makeblastdb/BLAST versions and arguments. Preserve per-target database
  artifacts and paths where required; verify database integrity. Keep deterministic
  sample/seed, subject order/count and `max_target_seqs` policy. Defer persistent
  caches until corruption, ownership and publication behavior are specified.
- [ ] B19.3 — Run `test_evidence.py`, `test_evidence_competitors.py`, subcommand tests
  and real competitor/evidence integration. Compare sampling manifests, raw hits,
  classifications, failure files and all DB bytes where stable. Benchmark cold build
  and repeated reuse, including hashing/materialization. Test invalidated/corrupt DBs.

**Stop:** cache hits suppress required validation/failures, native DB artifacts are
not stable, or fixed query/subject search results change. No cheaper incomplete search.

### B09 — Read only necessary layers; preserve the persisted count schema

**Priority:** large-memory opportunity; reader work only in this plan.
**Reach:** `scripts/multimap.py:final_results` defines aliases; consumers include
detection, UMAP, programmes, host response, validation and `menu._swap_multimap_layer`.

- [ ] B09.1 — Profile one consumer and enumerate every X/layer/obs/var/uns read,
  schema validation and fallback alias it needs. Freeze current H5AD bytes and
  alias identities after read. Start with programme/host-response readers, not the
  validator or method switcher which intentionally inspect broad state.
- [ ] B09.2 — Prototype consumer-local backed/selective reads or called-cell views,
  retaining required metadata/layers and the existing stored file. Close handles;
  preserve slicing dtype/order/copy and input immutability. Audit malformed omitted
  layers: if eager loading previously failed there, preserve the check or keep the
  eager fallback. Do not remove any of the 17 layer names or change writers.
- [ ] B09.3 — Run affected consumer tests, `test_validation.py`, `test_mech_a_consumers.py`,
  rerun-multimap and method-switch checks. Compare full downstream byte sets and input
  H5AD hashes, across aliases/missing layers/legacy and v3 encodings. Measure load,
  slice, complete consumer wall/RSS and file-descriptor usage.

**Stop:** omitted data hide errors or unsupported backed operations change behavior.
Persistent deduplication is **outside strict H5AD identity** and needs a separate
schema/migration decision; this card cannot label it byte-preserving.

### B01 — Compact allocation state while retaining entry order

**Priority:** substantial allocator work; seek memory reduction first.
**Reach:** `multimapping.py:build_multimap_layers`, `_csr_from_entries`, molecule
iteration/audit; `scripts/multimap.py` storage, method switching and lineage consumers.

- [ ] B01.1 — Profile live sizes of ambiguous records, coordinate lists, EC metadata
  and sparse construction on all six actual allocation methods. Freeze every layer's
  data/indices/indptr/dtype, conservation/audit/diagnostics and assignment bytes.
- [ ] B01.2 — First replace one high-volume Python coordinate list family with
  bounded typed append buffers preserving every entry and insertion order, followed
  by the **same single** sparse construction from the full ordered entry sequence.
  Include final concatenation transients. Preserve integer range/overflow behavior;
  retain fallback for unsupported ranges. Compact ambiguous record storage only in
  a second isolated step with identical iteration. No EC/cell aggregation, chunk
  sparse summation, diagnostics removal or new accumulation order.
- [ ] B01.3 — Run multimap/counting-contract/property/tolerance/audit tests and actual
  allocation refits for each method, with CB–UMI intersections, collisions, unresolved
  molecules, repeated coordinates, sibling groups and host-virus boundaries. Compare
  raw H5AD/count-audit/assignment outputs and typed intermediate buffers. Measure
  peak RSS during sparse construction, not only list growth, plus time and disk.

**Stop:** typed buffers or consolidation change sparse reduction order, diagnostic
aliases or errors; a more ambitious bounded allocator needs a new exactness argument.

### B02 — Reuse full-width cell-EM workspace before shrinking it

**Priority:** numeric risk; preserve fitted trajectories and convergence exactly.
**Reach:** `multimapping.py:em_cell_abundances`, `em_gene_abundances`, em-cell caller
in `build_multimap_layers`; fitted compact retained state, audit/assignment consumers.

- [ ] B02.1 — Capture each iteration's theta, update, denominator, convergence ratio
  and diagnostics on small independent original oracles. Cover no compatible unique
  support, zero prior/global sum, large gene count with few compatible genes,
  max-iteration exits, near-tolerance values, skewed priors and invalid inputs.
- [ ] B02.2 — First precompute repeated gene-index arrays and reuse **full-length**
  scratch buffers, resetting via operations proven bit-identical to `unique + prior`.
  Keep global-prior normalization, insertion order, full-width sums and difference
  reduction unchanged; returned arrays must not alias future mutations. Retain the
  early all-zero fallback. A compatible-only iteration is a separate derivation:
  it cannot replace full-width sums by summing only active entries without proof.
- [ ] B02.3 — Run all multimap properties/tolerance/counting-contract tests and exact
  per-iteration oracle cases; execute real em-cell/global refits and assignment
  reconstruction. Compare iteration counts, diagnostic reasons, numerical bits and
  final H5AD/TSV bytes. Benchmark many sparse-support cells and broad-support cells;
  measure scratch lifetime/peak RSS separately from already compact retained state.

**Stop:** convergence changes by even one iteration or one fitted bit. Existing
EMC-01 retained-state compaction is complete and must not be claimed as a new fix.

### B03 — Reduce assignment replay overhead without dropping lineage

**Priority:** applies only when molecule assignments are requested.
**Reach:** `multimapping.py:build_multimap_layers` local `bus_records`,
`_write_molecule_assignments`, `iter_cb_umi_molecules`; rerun-multimap safeguards.

- [ ] B03.1 — Profile path-backed versus DataFrame-backed replay, repeated stable
  sort, resolution and gzip writing. Freeze every method's weights, EC-list order,
  unresolved rows, row counts, metadata headers, failed input prefixes and input
  immutability. Confirm actual production input type before choosing a fix.
- [ ] B03.2 — First reuse the same lazily prepared stable DataFrame ordering across
  its two traversals, if that path is relevant; do not add memory on path-backed
  production runs. For the path-backed bottleneck, assess a bounded on-disk resolved
  molecule spool retaining every assignment field and exact order, then emit with
  the unchanged gzip/CSV writer after fitting. Retain original replay on unproven
  paths. Validate how reread I/O errors and externally changed BUS files behave;
  reuse requires an actual immutable input contract.
- [ ] B03.3 — Run multimap property/counting-contract, rerun-multimap and BUS-ingestion
  tests. Require assignment raw gzip bytes, status totals and all count layers to
  match for each real allocation method. Benchmark enabled/disabled lineage, spool
  disk/RSS and complete allocation-plus-writing; disabled output incurs no overhead.

**Stop:** gains require retaining all molecules in extra RAM, second-pass failures
are lost, or compressed bytes change. Never disable the requested assignment file.

### B21 — Bound dense host-feature construction without changing the model

**Priority:** major memory opportunity; solver changes excluded.
**Reach:** `scripts/hostresponse.py:_select_features`, `_run_cohort`, normalization,
HVG selection and cell/group/differential/enrichment consumers; entrypoint contract.

- [ ] B21.1 — Profile float64 dense intermediate versus final float32 storage and
  ineligible-cohort overhead. Freeze raw/log-normalized sparse/dense inputs, duplicate
  coordinates, normalization mutations, feature order and the current full-matrix
  warning. Eligibility-first reordering currently also moves warning/error/network
  annotation behavior: inventory those effects before proposing a skip.
- [ ] B21.2 — First fill a final C-order float32 array from bounded dense row blocks
  using the same **dense-then-cast** operation, avoiding a full float64 dense copy.
  Keep feature selection and solver input dtype/layout/order identical. For
  noncanonical sparse input retain the original path unless duplicate accumulation
  is proven identical. Do not cast sparse entries before duplicate summation, change
  to sparse solvers, move HVG learning outside training folds or reduce features.
  A lazy no-eligible-cohort path is a later separately qualified subchange.
- [ ] B21.3 — Run hostresponse, groups/strata/baseline/null/symbol and entrypoint tests.
  Compare feature array bits, predictions/weights, folds, statuses, warnings and all
  CSV/JSON artifacts. Exercise depth matching, panel-in-fold and multiple strata.
  Measure conversion-stage and full-cohort peak RSS, setup transients and wall time.

**Stop:** feature layout/dtype/numerical bits or observable preprocessing behavior
changes. This reduces transient memory, not the final dense matrix's unavoidable size.

### B22 — Reuse deterministic fit preparation; gate parallel execution separately

**Priority:** potentially expensive optional modeling; exact numerical environment.
**Reach:** `scripts/hostresponse.py:evaluate_group_cv`, `_run_stability_selection`,
`_run_l2_regression`, `_panel_depth_fold_prob`, `group_cv_permutation_null`,
`run_structured_null` and table writers.

- [ ] B22.1 — Measure actual fit counts and preparation cost by seed/fold/stratum/
  replicate. Freeze full RNG draws/state, canonical cell/fold ordering, failure rows,
  convergence warnings, stable-panel selections, per-fit outputs and final reductions.
  Null refits exclude baselines/panel arguments; do not multiply the wrong workload.
- [ ] B22.2 — First share only immutable, label-independent preparation proven equal
  within a fixed fit context, leaving seed/fold/iteration loops serial and draws in
  their existing order. Never carry observed label-dependent cohorts/HVG/folds into
  null replicates. Draft any parallel design separately: preplan ordered task seeds,
  preserve per-fit thread settings, ordered warnings/errors/results and identical
  final reduction. Worker/resource changes need separate authorization under the
  prior thread-budget constraint; do not silently enable a new default.
- [ ] B22.3 — Run hostresponse groups/baseline/null/strata tests and independent small
  full-procedure refits. Compare probabilities, coefficients, stability counts, null
  rows, empirical p-values and output bytes; injected fold/fit failures must retain
  the same first error. Measure full cohort time, total CPU and combined worker RSS
  if a later authorized parallel prototype is attempted.

**Stop:** RNG use, label isolation, warnings or floating reductions differ. Fewer
seeds, folds, stability iterations or requested null replicates are never this fix.

### B26 — Reuse digests only within a genuinely immutable operation

**Priority:** low until measured; integrity checking is required work.
**Reach:** `run_safety.py:sha256_file`, `sha256_path`, `build_run_manifest`,
`build_run_complete`; `validation.py:validate_run`; evidence/reference manifests.

- [ ] B26.1 — Count duplicate byte reads within one operation versus deliberately
  separate resume/completion/validation boundaries. Freeze changed bytes with same
  size/mtime, symlink retargeting, directory member changes, corrupt outputs, missing
  files and permission errors. Preserve sorted path-plus-hash directory identity.
- [ ] B26.2 — Choose one private call chain that already owns immutable snapshots;
  compute its digest once and pass the verified value with that captured payload.
  If no immutable boundary exists, retain rehashing. No global path cache, no reuse
  across completion or later `validate-run`, and no mtime/size-only fingerprint.
  Keep actual candidate software hashes truthful and preserve public helper behavior.
- [ ] B26.3 — Run run-safety, run-complete-marker, software-identity, reference-
  provenance, validation and evidence tests. Mutation tests must still detect every
  changed payload and stale completion marker. Compare fingerprints, manifests,
  validation issue ordering and errors exactly; benchmark complete validation/startup
  and bytes read, including snapshot/copy cost.

**Stop:** the only way to save I/O is weakening revalidation or missing concurrent
changes the original detects. Record no safe reuse rather than adding stale hashes.

### B28 — Reuse verified cached-accession data within an offline merge

**Priority:** low; deliberate network throttles remain unchanged.
**Reach:** `scripts/ncbi_fetch.py:_cache_valid`, `_gtf_cache_valid`, `_write_cached`,
`cached_accession_provenance`, `_fetch_one`, `fetch_reference`; `data_fetch.py`
integrity checks only if separately shown redundant. Coordinate with B25/B26.

- [ ] B28.1 — Run fully offline cache fixtures, counting hash/parse calls per
  accession. Freeze legacy/recorded retrieval status, GTF format stamps, accession
  dedup order, corrupt/missing sidecars, FASTA/GenBank mismatch, trailing-newline
  behavior, credentials redaction and failure/partial merged files.
- [ ] B28.2 — Pass verified digests/parsed metadata only for captured immutable bytes
  within the same merge. Keep each new invocation's cache integrity checks. Streaming
  merge is a separate B25-style subchange, staging all chunks before publication so
  late errors preserve original outputs. Do not change source snapshots, timestamps,
  remote request order, retry/backoff or throttle, and do not add network parallelism.
- [ ] B28.3 — Run `test_ncbi_fetch.py`, `test_data_fetch.py`,
  `test_reference_manifest_contract.py` and relevant CLI tests with network disabled.
  Compare merged FASTA/GTF/provenance/sidecar bytes, logged requests and failures.
  Benchmark cached end-to-end preparation with digest reuse overhead and peak RSS;
  mocked requests demonstrate contracts only, not remote-service performance.

**Stop:** reuse permits stale/corrupt cache data or changes failure publication.
No live requests or large downloads are authorized by this card.

### B12 — Profile and preserve the required emptyDrops computation

**Priority:** scientific-contract protection; measurement first, no current algorithm fix.
**Reach:** `scripts/cellcalling.py:_run_emptydrops`, `emptydrops_cells`, host-matrix
selection; `scripts/emptydrops.R`, `read_emptydrops_mask`, called-cell artifacts.

- [ ] B12.1 — Verify availability/version of the approved R runtime and DropletUtils;
  do not infer it from an old skipped or passing receipt. Capture original A/A TSV,
  mask, R warnings/stderr, RNG settings and conversion timing on tiny nondegenerate
  kb/STARsolo fixtures. Preserve exact `matrix.mtx` preference and square-matrix
  orientation behavior. Profile parsing/conversion/emptyDrops/barcodeRanks separately.
- [ ] B12.2 — Only if adapter/conversion cost is material, prototype one storage-
  lifetime change preserving identical `gxc` sparse buffers and the exact sequence
  `set.seed` → `emptyDrops` → `barcodeRanks` → output construction. Keep caller,
  FDR/lower/niters/seed and numerical libraries unchanged. Cache/reuse across runs
  requires its own identity/failure/provenance design and is not the first fix.
- [ ] B12.3 — Run `test_cellcalling.py`, infected-rule/sensitivity tests and actual R
  fixtures. Compare raw TSV, exact FDR/total/knee/inflection bits, called masks and
  downstream detection bytes. Include gzip/plain triplets, dimension errors and
  threshold-boundary cells. Measure total detection plus R child peak RSS/CPU.

**Stop:** R dependencies unavailable, original A/A not repeatable, or Monte Carlo
dominates with no safe adapter gain. Record measured unavoidable cost; never reduce
niters, substitute knee calling, or claim a skipped R test qualifies this path.

### B11 — Reconcile counting resource declarations without changing workers yet

**Priority:** resource correctness; separate from strict-byte optimizations.
**Reach:** `Snakefile:kb_count`, `multimap`; `scripts/multimap.py` correction/sort
commands; `RunConfig.cores`, CLI/menu configuration and command/provenance tests.

- [ ] B11.1 — Record actual installed kb/kallisto/bustools defaults and effective
  workers against Snakemake declarations at requested cores 1/6/8/32. Use command
  inspection/dry runs first. The historical kb default of eight is a pinned older
  observation, not a timeless assumption. Record oversubscription and concurrency.
- [ ] B11.2 — Draft explicit semantics for total workflow budget versus native job
  threads, including both whitelist branches and Snakemake's resource scaling.
  Test a fixed-effective-worker version first. **Do not implement changed thread
  counts under the earlier no-thread-change instruction.** Obtain separate scope
  authorization before wiring requested workers or altering scheduling defaults.
- [ ] B11.3 — After authorization, run DAG/config/CLI/multimap-tool tests plus real
  tiny counting/strand/lineage integration. Compare native BUS/EC/H5AD artifacts,
  failures and reports with original A/A. Measure worker counts, CPU utilization,
  total stage wall/RSS and concurrent budget adherence. Any thread-dependent raw
  bytes require an explicit separate compatibility decision, never normalization.

**Stop:** authorization absent or native identity cannot be established. Resource
correctness can be valuable without being a byte-identical performance improvement;
keep those claims separate.

### B30 — Make reference worker controls explicit only behind the resource gate

**Priority:** resource correctness; reference-build-only.
**Reach:** `scripts/build_reference.py:_run_cdhit_est`, `measure_host_homology`,
`_finish_reference_build`, `_run_kb_ref`; `menu.py:_build_kb_ref`, build CLI/config.

- [ ] B30.1 — Inventory `cd-hit-est -T 0 -M 8000`, minimap2 defaults and both kb-ref
  entrypoints on the installed binaries. Profile stage CPU/RSS and captured PAF
  volume without a large reference build. Freeze clustering representative/tie
  order, homology records and prepared/native reference bytes.
- [ ] B30.2 — Draft one explicit resource policy covering all entrypoints and
  provenance. Preserve clustering/search/scientific options and current effective
  workers during byte qualification. Changed thread counts require separate
  authorization as in B11. Streaming the captured PAF is a separate memory subchange
  with unchanged parse/order/error contract; do not combine it with worker rewiring.
- [ ] B30.3 — After authorization, run build-reference/manifest/pipeline tests,
  host-homology-loci tests and real tiny CD-HIT/minimap2/kb builds when tools exist.
  Compare representative FASTA/cluster files, PAF, homology TSV, index/t2g/native
  artifacts and downstream count bytes. Report actual workers and child peaks;
  retain the existing `max_over_all_child_processes` RSS scope label where used.

**Stop:** worker changes alter cluster/tie/native bytes or lack authorization.
No raised CD-HIT memory cap, changed identity cutoff, tool version or large build
is implied by this plan.

## Completed candidates to preserve

These are implemented, not new open optimization cards. Their qualification scope
and limits remain those of E01–E05; no universal or complete native-run claim follows.

| Candidate | Status and current implementation | Mandatory checks on every later combined change |
|---|---|---|
| B05 | Complete: bounded canonical CSR subset cache inside `summarize_multimap_evidence`; old paths retained elsewhere | Same dense column arithmetic/order, missing-layer aliases, canonical/noncanonical/CSC/dense/generator/over-budget fallbacks, untouched input buffers and peak memory. Broad CSC conversion was not adopted. |
| B07 | Complete: guarded scalar count formatting in `_count_value` | Rounding/tolerance, exact return types/bits, nonfinite/errors/NumPy error policies, huge-value fallback and serialized tables. |
| B20 | Complete: binary bool/integer ndarray enrichment fast path | Original neighbor query/positive set, random draws/permutation order, final reductions, exact bits/Generator state and unsupported-input fallback. |

Run `test_performance_exactness.py` against its independent frozen original helpers;
do not regenerate that oracle from modified code. Reuse the old 689-file corpus and
12 retained-slice comparisons where applicable, but extend fixtures for every newly
touched path. The original stage comparison had four truthful provenance differences;
report them as differences, not whole-directory byte identity.

## C00 — Combined qualification and completion criteria

Dependency: every selected card independently passed its exactness and G02 benefit
gate. Deferred cards need documented failed/pending gates, not checked edit steps.

1. Freeze all accepted source hashes before final runs. Compare each isolated card
   and the combined source against G00's same immutable baseline; compare the
   pre-first-wave original where its artifact contract applies. Run integration at
   overlap boundaries: B04/B06/B08/B05 detection; B09/B21 host features;
   B15/B16/B17/B18 evidence; B24/B25/B28 reference preparation. One change must not
   conceal another's drift.
2. Run the selected cards' focused tests, first-wave exactness suite, then the
   standard package suite using the configured network/integration/research policy.
   Run required native/R integration separately; label unavailable checks as pending.
   Use targeted Ruff/format, mypy when Python changes warrant it, and governance/
   ship-scope checks for affected files. Do not repair unrelated failures silently.
3. Repeat paired timing and memory trials on the frozen combined snapshot. Measure
   any authorized representative whole run separately; no component ratio products,
   hidden skipped diagnostics, lower iteration budgets or changed threads.
4. Save final source/input/tool/harness hashes, exact comparisons, A/A variability,
   truthful provenance spans, test commands/exits, individual trials, memory method,
   native/R gaps and rollback decisions. Check current live hashes still match the
   measured candidate. A late edit invalidates affected qualification.
5. Update this tracker and `PLAN.md` with passed, deferred and pending work. Report
   analytical exactness separately from operational bytes, resource correctness,
   scientific validity and release readiness. No commit/push/publication is implied.

Rollback: restore only the selected card's isolated source diff after checking for
concurrent edits, or leave the previous implementation/fallback active. Never reset
the shared tree, delete retained outputs, rewrite old goldens or loosen a comparator.
An exactness failure blocks that candidate even if its speedup is large.

## Suggested commands and execution handoff

Run from the repository root in the already installed full test environment. The
example B24 command is focused and offline; replace the test list with the selected
card's list. These commands are instructions for execution, not receipts of runs
performed while drafting this plan.

```bash
PERF_PYTHON=/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/viralscan_test_full/bin/python
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8 TZ=UTC
export MPLCONFIGDIR=/tmp/viralscan-perf-followups-mpl
export NUMBA_CACHE_DIR=/tmp/viralscan-perf-followups-numba
"$PERF_PYTHON" -m pytest tests/test_build_reference.py tests/test_reference_pipeline_aud15.py tests/test_reference_manifest_contract.py tests/test_bundled_panel_mask.py
```

The single-thread library settings reproduce the synthetic test environment; they
are not permission to change a production run's native or library worker counts.
Snapshot-specific harnesses must set `PYTHONPATH` to the relevant immutable snapshot
and assert module paths/hashes before computation. Use unique receipt labels. The
old `execution/qualify.py` covers the first wave only, not all cards in this document.
Do not run survey scripts against drifted live source by disabling their hash guards.

Execution prompt:

> Read this plan and its linked evidence completely. Start G00–G02, then action one
> eligible card at a time in priority order. Inspect actual affected bodies/callers/
> serializers/tests and capture baseline A/A goldens before production edits. Keep
> B05/B07/B20 intact. Make an isolated minimal change, require exact outputs and
> failure behavior, measure complete stage time and isolated peak memory, and defer
> any card that fails its gates. Preserve shared work and historical receipts.
> Resource/scientific-contract changes remain gated; no SLURM submission, large
> download, retained-data deletion, commit, push or publication. Finish with C00 and
> durable evidence identifying accepted changes and every remaining limitation.

## Decision and evidence log

- 2026-10-17 — Planning review complete. Actual source checked against all 76
  captured combined files. All 30 IDs reconciled: 27 unchecked implementation
  cards and three completed preservation entries. Added explicit gates for B08's
  negative paired timing, B10's validation behavior, numerical reduction order,
  native variability and resource authorization. No production edits or performance
  experiments performed during this planning task.
- Pending execution — Append dated entries per card with packet path, baseline/
  candidate hashes, gate results, changed scope, failures, decisions and limitations.
  Reorder priorities only with recorded evidence. G00/G01/G02/C00 remain not started.

### Planning validation receipt

2026-10-17 — Local document/source checks passed using the full test interpreter:
27 unique open B-cards, 81 unique ordered B-substeps, all 30 survey IDs accounted
for, every open card assigned a wave, all Markdown link targets present, and all
33 literal test filenames found under `tests/`. Embedded Bash passed `bash -n`.
The document is below the 500 KiB repository limit, with LF endings, final newline,
no tabs or trailing whitespace. The current 76 package/resource hashes still match
the combined snapshot. Ship-scope helpers confirm this plan is outside the shipped
allowlists; no package/schema/export configuration was changed.

`PYTHONDONTWRITEBYTECODE=1 <full-test-python> scripts/check_data_governance.py`
returned exit 0, `data-governance check passed`. `git diff --check -- PLAN.md`
passed. These are planning/hygiene checks only. No package regression suite, new
candidate benchmark, native integration or R qualification was executed for this
documentation-only task. The first-wave receipts remain unchanged.
