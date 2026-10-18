# Remaining performance candidates: implementation and qualification plan

Date: 2026-10-17; status reconciled 2026-10-18 (`PERF-DOC-01`). Project: ViralScan.
Status: **wave 1 (B24, B23, B04, B06) implemented and qualified at component level on
`codex/perf-wave1-b24-b23-b04-b06`; B08 deferred on its benefit gate; B10 and wave 2
(B25, B27, B18) in flight in the `20261018-112322-wave2-g00` packet; the remaining cards
are not started. G01, G02 and C00 stay partial.**

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
- [x] G00 — Capture the current optimized baseline. **Status: done (2026-10-17).** Frozen,
  read-only `baseline_source/` (= 5e551c9 + the uncommitted first-wave edits), `baseline.json`
  and `capture_baseline.py` under `analysis/computational_bottlenecks_2026-10-17/followups/20261017-122240-G00/` (untracked evidence packet).
- [~] G01 — Freeze new originals and original A/A behavior. **Status: partial.** Exact-mode
  harness sections b24/b23/b04/b06 ran baseline-vs-candidate with identical `file_manifest`
  per card and in C00 (baseline A/A twice, byte-identical); frozen oracles are in each card's
  `extra_tests`. The comparator mutation self-check receipt was not produced in packet 1; the
  wave-2 packet produced one for its own baseline (2026-10-18,
  `followups/20261018-112322-wave2-g00/g01_comparator_selfcheck.json`, verdict pass: one-byte,
  missing-file, dtype, one-ulp and error-message mutations all rejected; an identical failure
  reproducing and warning records captured are accepted, as intended). The wave-1 goldens still
  lack their own mutation receipt.
  **Wave-2 update (2026-10-18):** re-frozen for the wave-2 baseline (wave-1 tip `faf5667`)
  in packet `...followups/20261018-112322-wave2-g00/`: 80-file `baseline_file_manifest.tsv`,
  A/A identical (`g01_baseline_aa.json`) and the comparator mutation self-check now exists
  (5/5 injected mutations rejected, `g01_comparator_selfcheck.json`). Still single-session.
- [~] G02 — Measure relevant stage shares. **Status: partial.** Paired timing/memory exists
  per card and in C00 on synthetic component fixtures (single session of 3 baseline + 3
  candidate processes). No stage profile of a real run, no second timing session, no
  representative whole-run measurement. **Wave-2 update (2026-10-18):** the B08 replot is
  9 alternating paired trials per fixture in one session at host load 17.1 → 9.4; the second
  session requirement is still outstanding.
- [~] C00 — Qualify combined accepted changes. **Status: wave-1 (B24, B23, B04, B06)
  qualified at component level on 2026-10-18; native/R/whole-run checks pending.** Evidence:
  `analysis/computational_bottlenecks_2026-10-17/followups/20261017-122240-G00/c00/C00_REPORT.md`. Next: remaining cards, a second timing session on a quieter host.
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

- [x] B24.1 — Freeze valid/invalid gate calls and count invocations. Include empty
  FASTA/record, duplicate IDs and sequences, duplicate IDs with distinct sequences,
  lowercase/N/IUPAC, short records, all classifier categories, exact fraction and
  homopolymer limits, tied worst records and more than ten tied offenders. Capture
  the current report's last-value/first-key behavior for duplicate identifiers.
  **Evidence (2026-10-18):** gate fixtures (empty/header-only, duplicate ids and sequences, N/IUPAC, short records, ties, >10 tied offenders, nextafter limits, invocation counts) frozen in `tests/test_perf_wave1_b24.py` against a verbatim copy of the 5e551c9 module.
- [x] B24.2 — Compute per-record integer counters once within one gate invocation;
  derive the existing report, pure total and threshold verdict from them. Reuse
  structural/verdict helpers privately while retaining public signatures. Keep
  initial scan logging before structural failure, structural `ValueError` outside
  the low-complexity `RuntimeError` wrapper, stable offender sorts and causes. No
  process-global cache, classifier rewrite, mask change or threshold change.
  **Evidence (2026-10-18):** one classifier call per record per gate (3N -> N; 2N -> N on structural failure), FASTA read once; public signatures unchanged; verdict body reused verbatim.
- [~] B24.3 — Run `test_build_reference.py`, `test_reference_pipeline_aud15.py`,
  `test_reference_manifest_contract.py`, `test_bundled_panel_mask.py`; add strict
  diagnostic/order oracle checks. Compare prepared FASTA/GTF/report bytes and all
  error/log records. Time whole gates on clean random and N-rich/masked panels at
  several lengths, then offline reference preparation. Native build claims require
  real identical-tool builds; retain fallback if shared helper behavior drifts.
  **Evidence (2026-10-18):** four named test files (124 passed, 3 dustmasker skips) plus 65 oracle tests; exact harness `file_manifest` identical (43 artifacts) in C00 incl. error/log records; gate timing C/B 0.31-0.39 over 11 jobs (B spread 6-18% on seven clean jobs, 78-117% on four jobs with a baseline outlier process; medians unaffected); preparation-inclusive C/B 0.27-0.33 with B spread 4-7% in one session, but one candidate process (started at load 61) overlapped baseline in one job, so session-2 repeat is still advisable. Real dustmasker/`kb ref` native builds not run (binary absent); masking exercised with a faked `_run_dustmasker`. Memory unchanged (traced peak -1 to -6%).

**Stop:** a repeated scan cannot be removed while preserving duplicate/error/log
semantics, or preparation-inclusive benefit fails G02. Treat targeted masking's
per-window classifier optimization as a later separate change, not this card.

### B23 — Avoid programme dictionaries that the writer discards

**Priority:** early optional-feature candidate; benefit unmeasured.
**Reach:** `gene_programs.py:call_cell_programme`, `_nonzero_columns`, `_decide`;
`scripts/gene_programs.py:_marker_matrix`, `run_one`, `main`; summary/log serializers.

- [~] B23.1 — Profile marker-negative cohorts and freeze the exact prune predicate:
  `state == indeterminate` and zero productive/latent breadth on both layers. A
  zero matrix alone does not establish that a call is prunable. Include selected-only
  evidence, unresolved markers, absent models, absent layers, overlapping groups,
  min-UMI boundaries, non-assessable absence and latent/reactivating flags.
  **Evidence (2026-10-18):** prune predicate (indeterminate and zero breadth on both layers) frozen by 129 oracle tests incl. tie order and zero-cell validation; no separate profile of marker-negative cohorts beyond the timing/memory fixtures.
- [x] B23.2 — Keep the public helper's all-cell return contract. In the script path,
  first try bounded row blocks with the unchanged helper and immediate pruning;
  accumulate state counts for **all** calls, including discarded rows, in original
  first-seen order. Do not log per block. Reconstruct the existing pandas value-count
  tie order for the final log or retain that path. If only memory improves, qualify
  that result before attempting an index-skipping helper with a broader API change.
  **Evidence (2026-10-18):** `run_one` calls `call_cell_programme` in row blocks of 8192, prunes immediately, counts all states in first-seen order, logs once; public helper untouched. C00 adds a `None`-layer pass-through (`unique if unique is None else unique[start:stop]`) so mypy stays clean and a `None` layer is validated by the helper exactly as before.
- [~] B23.3 — Run `test_gene_programs.py`, `test_gene_programs_audit4.py`,
  `test_gene_programs_unresolved.py`, `test_program_evidence_boundary.py`,
  `test_gene_programs_capturability.py`; compare cells/summary/JSON/manifest bytes,
  marker-resolution metadata and logs. Benchmark sparse-negative and dense-positive
  called-cell cohorts, including helper work and serialization; measure peak RSS.
  **Evidence (2026-10-18):** focused tests 157 passed; exact harness identical (96 files, 2184 cases) with block sizes None/1/2/3/7. Memory (traced, 150k cells): sparse-negative 0.20x, dense-positive 0.82x of baseline; wall time at 30k cells C/B 0.52 for both cohorts but baseline-vs-baseline spread 44-48% (one baseline process was as fast as the candidates), so no speedup is claimed. Process RSS is shared with other sections and is not an isolated measure.

**Stop:** summaries/log tie order or helper validation timing changes; keep all-cell
calls rather than dropping scientifically meaningful zero-evidence states.

### B04 — Cache a narrow local detection view, keep dense column reductions

**Priority:** early measured hotspot; no approved direct sparse/batched sum.
**Reach:** `scripts/detection.py:_gene_counts_from_matrix`, `detect_genes`,
`preprocessing`; `found_genes.txt`, group detection and downstream reports.

- [~] B04.1 — Profile post-first-wave detection, recording requested viral-column
  fraction. Freeze duplicate-name validation, repeated accession traversal, Python
  group-sum order, scalar result dtypes and nonfinite error lists. Reuse adversarial
  sparse-sum counterexamples from the first-wave probe as mandatory regression cases.
  **Evidence (2026-10-18):** mandatory sparse-sum counterexamples and duplicate/nonfinite/dtype fixtures frozen (600 oracle tests); the requested viral-column fraction was not profiled on real detection data.
- [x] B04.2 — For reusable sequences and a proven supported layout, select requested
  columns once into a bounded **local** CSR view; map each original index to it.
  Continue the same individual `toarray().sum()` with its original shape/dtype and
  order. Preserve input matrices and the duplicate-name check before access. Fall
  back for generators, unsupported/noncanonical representations, all-column cases
  and over-budget views, as needed. Do not import B05's cache assumptions untested.
  **Evidence (2026-10-18):** `_requested_column_view` selects the requested columns once into a local CSR copy for plain CSR input; per-column `toarray().sum()` unchanged; fallbacks retained. Cost: traced peak +4 MB at 100k x 2000 x 300 viral columns (1.36 -> 5.35 MB).
- [~] B04.3 — Run `test_detection.py`, `test_var_names_safety.py`,
  `test_infected_rule_and_gate.py`, `test_detection_audit_s2.py` and stage regressions.
  Cover dense/CSR/CSC, all numeric fixture classes, repeated/missing genes and group
  threshold boundaries. Compare typed totals plus actual detection-stage bytes;
  benchmark complete preprocessing with cache construction and high selected fractions.
  **Evidence (2026-10-18):** 239 focused tests + 600 oracle tests pass; exact harness identical (57 cases); B04+B06 overlap run (detect_genes -> compute_stats with forced 1/few-row blocks, 110 layout/dtype/budget cases) identical to baseline. Timing C/B 0.34 (grouped) and 0.56 (per gene) but B spread 56-79%, so benefit inconclusive under the noise rule. Full `preprocessing()` stage bytes not run.

**Stop:** any total bit/dtype, group decision or error differs, or cache memory defeats
the measured benefit. Never accept a sparse reduction because `allclose` passes.

### B06 — Bound dense row-reduction memory

**Priority:** early memory candidate; direct sparse sum has failed exactness probes.
**Reach:** `utils.py:virus_molecules_per_cell` and its detection/enrichment/UMAP/plot
callers; separately `compute_stats` accession-breadth and ambiguity-total operations.

- [x] B06.1 — Characterize existing dense layout/strides and row-reduction dispatch
  for each supported sparse representation/dtype and selected gene count. Freeze
  cancellation/large-positive mixtures, signed zeros, float32-to-float64 conversion,
  repeated/missing genes, zero rows/columns and infected-threshold neighbors.
  **Evidence (2026-10-18):** bitwise oracle across layouts/dtypes (bool/int/uint/float32/64, complex, longdouble, float16 dense), signed zeros, cancellation, >2**53 integers; 384 tests.
- [x] B06.2 — First optimize only the shared per-cell helper: dense **row blocks**,
  with all selected columns in the same order, the same float64 conversion, row
  width and array layout, followed by the same `sum(axis=1)`. Never split a row's
  reduction across column blocks or combine partial sums. Choose a private byte
  budget, guard oversized single rows, and retain the old path for unproven layouts.
  Leave flattened scalar ambiguity sums unchanged; they are a different reduction.
  **Evidence (2026-10-18):** row-blocked densify (32 MiB private budget, whole rows, same float64 C-order reduction), original path for CSC/dense/unproven layouts and any nonfinite sum. **Round-2 guard found by the adversarial verifier:** the blocked path applies only to bool/int/uint/float32/float64; complex (per-block ComplexWarning multiplicity), longdouble and object dtypes keep the single-shot path.
- [~] B06.3 — Run detection, infected-rule, enrichment/UMAP and performance-exactness
  tests; exercise block sizes 1, 2, non-divisors and larger than row count across
  widths around NumPy reduction boundaries. Require exact arrays, warnings and all
  caller artifacts. Benchmark wide-virus and narrow-virus cases; show reduced
  isolated RSS including temporary conversion buffers. Treat breadth-only streaming
  as a separately isolated subchange after helper acceptance.
  **Evidence (2026-10-18):** 266 focused tests + 384 oracle tests pass; exact harness identical (96 cases incl. forced-blocked). Memory benefit reproduced in C00: isolated RSS delta 0.42x (100k x 200) and 0.085x (100k x 1000), deterministic (B spread 0-0.3%); time neutral (0.86-1.07, not separated). Not run: byte comparison of full enrichment/UMAP caller artifacts beyond the repo tests; breadth-only streaming stays out of scope.

**Stop:** blocking changes row reduction bits/layout-sensitive behavior, or warning
category/message/multiplicity cannot be preserved for allowed nonfinite inputs.
Nonfinite inputs can use the original path only if the guard itself preserves errors.

### B08 — Reuse the comparable-cell denominator lazily

**Priority:** conditional; later pilot did not improve paired reporting.
**Reach:** `scripts/detection.py:compute_stats`, `_comparable_called_cells`,
`_in_barcode_set`; host-set and combined-reference denominators.

- [x] B08.1 — Reprofile after B07/B05 on many-virus combined references and external
  host sets. Repeat the stats-only **and** paired reporting measurement in
  `execution/assess_followup.py` using a new packet. The saved 20→1 helper count and
  stats-only gain do not qualify the paired 1.7338→1.9760 s result as a speedup.
  **Evidence (2026-10-18):** 9 alternating paired trials per fixture at host load 17.1 → 9.4
  (`...followups/20261018-112322-wave2-g00/cards/b08/reprofile_b081.json`). Paired reporting
  e06 1.1884 → 1.1077 s (+7.3%), many-virus-200x10 14.1581 → 13.6397 s (+3.8%); stats-only
  +6.2% / +5.7%; prototype faster 9/9 trials on both helper-bearing fixtures; external
  host set 0.0% (the helper is never called when the caller supplies `host_cells`); traced
  peaks unchanged (19.8/59.1 MB). Both measurements improve where the helper runs, so the
  old paired regression is not reproduced — but the gain is below G02.
- [ ] B08.2 — Only after renewed benefit, initialize a local sentinel and calculate
  the identical denominator at the first original use inside a valid-virus iteration;
  reuse thereafter. Distinguish uncomputed from computed `None`. Do not hoist above
  earlier failures or evaluate it for zero viruses/all-invalid gene lists.
  **Deferred (2026-10-18): failed gate G02.** Renewed paired-reporting gain is 3.8–7.3%
  median affected-stage reduction (below the 10% default) with no memory reduction, so the
  card's stop condition applies; the saved helper count is explicitly not an acceptance
  criterion. No source edit exists. Evidence: `cards/b08/decision.md`, `feasibility.md`.
- [ ] B08.3 — Run detection/sensitivity/ambient tests and stage byte comparisons.
  Cover missing/empty host dictionaries, called/comparable `None` and empty sets,
  no host genes/no called cells, malformed masks, candidate-negative viruses,
  threshold boundaries and original caller warning order; assert input immutability.
  Include conversion/preparation/serialization in the benefit measurement.
  **Not run (2026-10-18):** qualifies an edit that was deferred above.

**Stop:** renewed paired reporting gain fails G02, any lazy branch/error changes, or
custom/stateful supported inputs invalidate reuse. A smaller helper count alone is
not an acceptance criterion.

### B10 — Determine whether disabled UMAP can avoid eager loading exactly

**Priority:** feasibility first; old missing/corrupt/schema failures are mandatory.
**Reach:** `scripts/umap.py:main`, `run`; `Snakefile` UMAP rule; H5AD and found-gene
parsing, run identity, stdout and `log/umap.done`.

- [x] B10.1 — Freeze enabled/disabled cases with missing/invalid/truncated H5AD,
  malformed X/layers/obs/var/uns, wrong schema, unreadable/malformed found-gene text,
  empty files, missing identity, unavailable Plotly and two simultaneous faults.
  Establish exactly which failures eager `read_h5ad` exposes and in what order.
  **Evidence (2026-10-18):** 37 runs / 27 failures through the real `umap.run`
  (`...followups/20261018-112322-wave2-g00/cards/b10/failure_matrix_b101.json`). Precedence:
  identity → Plotly (enabled only) → H5AD read (h5py open, then per-element anndata validation
  and encoding registry) → multimap schema → found-genes parse → plots. Two-fault cases fix the
  order. Read-tolerated (and so must stay tolerated): deleted X/var/obsm groups, zeroed CSR
  `indptr`, duplicate `obs` index, empty/extra-segment found-genes, missing Plotly when disabled.
- [ ] B10.2 — Prototype a reader-only bounded/lazy validation path for disabled
  plots, retaining every observed validation/error and found-gene conversion. Inspect
  backed AnnData behavior against actual installed code; `backed='r'` or reading
  `uns` alone is not automatically equivalent. Always close handles. Retain eager
  loading for unproven encodings. A prior detection marker does not prove present
  bytes are intact; no blind early return or trust-by-mtime shortcut.
  **Deferred (2026-10-18): no compliant fast path.** Six candidates were driven through the
  unmodified `umap.run` against the freeze (control reproduced it 26/26). `uns_only` silently
  accepts 7 malformed inputs, `backed='r'` 3 plus an error-order change, `hand_min` 5,
  `read_lazy` needs `xarray` (unusable), and even an adversarial hand-rolled validator only
  matches 24/26 by duplicating anndata's messages and private registry. Evidence:
  `cards/b10/probe_b102.json`, `feasibility.md`, `decision.md`.
- [ ] B10.3 — Run `test_umap.py`, `test_snakemake_script_modules.py` and relevant
  `test_snakefile_dag.py` cases. Require identical done-marker content/presence,
  stdout, errors and untouched input hashes. Measure complete disabled rule wall/RSS
  on retained H5AD encodings; enabled plotting remains on the original path.
  **Not run (2026-10-18):** qualifies an edit that does not exist; the frozen matrix
  (`failure_matrix_b101.json`) is retained as the disabled rule's regression net instead.

**Stop:** equivalent validation still needs the entire read or a malformed case is
silently accepted. Record **no compliant fast path** if so. Relaxing disabled-mode
validation requires a separate explicit behavior decision; this plan does not grant it.

### B25 — Stream reference text without changing byte or failure contracts

**Priority:** memory opportunity, start with one adapter.
**Reach:** `menu.py:_materialize_reference_input`; `scripts/build_reference.py`:
`_fasta_records`, `prepare_reference_inputs`, header/GTF readers; NCBI merge is B28.

- [x] B25.1 — Measure overlapping raw-byte/text/list/record allocations and choose
  the dominant adapter. Freeze CR/LF/CRLF and trailing-newline runs, blank inputs,
  descriptions, lowercase, Unicode text errors, malformed GTF and duplicate IDs.
  Note that `read_text().splitlines()` recognizes separators beyond a newline iterator.
  **Evidence (2026-10-18):** `...followups/20261018-112322-wave2-g00/cards/b25/`
  (`allocations_b251.json`, `contract_b251.json`). Size sweep: whole-byte concat peaks at
  3.0× the largest input (48 MB for 16 MB files, unbounded) while a bounded form is flat at
  2.01 MB; `_fasta_records` 1.73×, GTF reader 3.16×, header reader 2.22× of their inputs.
  Dominant adapter = the whole-byte concatenation (only unbounded term, byte-contract-only
  replacement). Frozen: 15 edge inputs × 10 combinations byte-identical, trailing runs 0–8,
  chunk-boundary runs, `IndexError` on a bare `>`, `UnicodeDecodeError`, malformed-GTF
  `ValueError`, duplicate IDs, and the eight `splitlines()`-only separators.
- [x] B25.2 — First replace whole-byte concatenation with bounded copying that
  buffers/discards only trailing CR/LF, emitting exactly one LF after each nonempty
  stripped input. Then separately assess streaming headers/GTF with original text
  splitting rules. Preserve uppercase/wrapping/record order and validation-before-
  write boundaries; spool before publication where required. Keep public list-return
  APIs until their callers have a separately qualified streaming implementation.
  **Evidence (2026-10-18):** implemented in `menu._materialize_reference_input` (per-file
  streaming, trailing CR/LF holdback, `_REFERENCE_COPY_CHUNK`). Headers/GTF assessed and
  **deferred**: a streamed header reader reproduces values exactly (incl. all eight extra
  separators) and cuts peak 29×, but chunked incremental UTF-8 decoding reports
  `UnicodeDecodeError` relative to the failing chunk (chunk-size-dependent messages), so
  keeping the failure contract needs a hand-built error; list-return APIs stay. Evidence:
  `assessment_b252.json`, `cards/b25/feasibility.md`.
- [x] B25.3 — Run build-reference, reference-manifest, reference-pipeline-aud15 and
  CLI tests; compare raw prepared/merged FASTA/GTF and manifest bytes, exception and
  partial-output behavior. Benchmark large synthetic host text and small panels,
  reporting isolated peak RSS, spool footprint and total preparation time.
  **Evidence (2026-10-18):** 284 passed / 3 skipped / 1 deselected
  (`pytest_b253.log`) plus 1324 passed governance+wave-1 (`pytest_b253_governance_wave1.log`;
  no re-pin needed — `menu.py` is not digest-pinned). `exact_b253.json`: 5/5 identical
  (multi-file merges, single-file bypass, empty-argument errors, missing-later-input partial
  output, `prepare_reference_inputs` prepared FASTA + manifest with timestamps excluded).
  `memory_b253.json`: traced peak 2.01 MB flat vs 16/64/128 MB baseline → 7.96×/31.84×/63.68×
  lower, byte-identical, wall ratio ~0.5. Oracle detection power: `mutation_b252.json`
  (4/4 distinguishable mutants caught; one mutant byte-equivalent by construction).

**Stop:** newline semantics or failure side effects differ, or a sequence-storage
refactor is needed to achieve the proposed gain. Do not bundle parser redesign here.

### B27 — Reuse bounded whitelist matching preparation

**Priority:** low, evidence first; ordinary read sampling is already bounded.
**Reach:** `chemistry.py:sample_r1`, `match_rates`, `detect`, `classify`, `resolve`;
chemistry preflight callers. Keep the separate full-list preflight API unchanged.

- [x] B27.1 — Profile many samples and candidate lists, counting repeated barcode
  lengths and preparation cost versus actual whitelist I/O. Freeze duplicate entries,
  `-1` suffix stripping, blank/short first lines, gzipped lists, empty/short reads,
  tied rates, label order and unreadable paths.
  **Evidence (2026-10-18):** `...followups/20261018-112322-wave2-g00/cards/b27/`
  (`profile_b271_s1.json`, `profile_b271_s4.json`, `contract_b271.json`). Real bundled
  on-lists (10xv2 0.74 M, 10xv3 6.79 M, 10xv4 7.37 M lines) plus a `-w` list, 100 k
  reads/sample: preparation is 1.24 % (1 sample) / 2.23 % (4 samples) of the call and
  streaming is 97.8–98.8 %; 10xv2/3/4 all infer `cb_len` 16, so the repeated-length
  count is maximal (one distinct length for four labels). Frozen: duplicates count
  twice (rate 1.5), `-1` stripped, blank first line → `cb_len` 0, short first line
  truncates reads, short reads excluded from counts, gzip opened by suffix, empty
  reads → 0.0, ties keep the first label in insertion order, label order preserved,
  `FileNotFoundError`/`IsADirectoryError` carry the path.
- [ ] B27.2 — Try per-invocation reuse of sampled barcode Counters by inferred
  barcode length, retaining separate original per-label whitelist scans and ordered
  integer hit accumulation. Bound cache memory; build lazily at the original access
  point. Do not deduplicate whitelist lines: duplicates currently contribute more
  than once. Do not retain million-entry Python whitelist sets or cache across calls.
  **Not implemented — negative (2026-10-18).** Tried as a prototype (identical rates,
  measured 0.990 ×/0.985 × the baseline call): the removable preparation is only
  0.91 %/1.66 % of the preflight, below G02's 10 %, because list I/O dominates. The
  cost side also fails the card's stop condition in the worst case: with four distinct
  inferred lengths the cache retains every label's counters at once and traced peak is
  **1.95 ×** the baseline (168.2 MB vs 86.2 MB). No source edit. Evidence:
  `cards/b27/decision.md`, `feasibility.md`, `worstcase_b271.json`.
- [x] B27.3 — Run `test_chemistry.py`, `test_chemistry_check.py`,
  `test_whitelist_preflight.py`; compare exact rate floats, detections, warnings,
  requested-chemistry errors and output text. Measure complete preflight, worst-case
  distinct lengths and memory; if I/O dominates, record a negative result.
  **Evidence (2026-10-18):** 53 passed (`pytest_b273.log`) with `chemistry.py`
  unchanged. Complete preflight 11.73 s (1 sample) / 25.05 s (4 samples), essentially
  all whitelist I/O; worst case for distinct lengths and its memory cost is in
  `worstcase_b271.json`; the change-vs-baseline comparison was made against the
  prototype, and "if I/O dominates, record a negative result" is the outcome.

**Stop:** benefit is below G02 or cache memory exceeds the agreed baseline budget.
Already cached native-tool discovery is not new work to optimize.

### B18 — Reduce per-read verdict state while preserving first-seen order

**Priority:** bounded-memory evidence work; no sort-by-name shortcut.
**Reach:** `molecule_verdict.py:best_scores`, `molecule_rows`,
`write_molecule_verdicts`, `write_from_bam`; evidence/dossier consumers.

- [x] B18.1 — Profile read versus molecule cardinality. Freeze order based on the
  first **usable scored alignment**, not first raw QNAME occurrence. Include secondary
  hits, supplementary/unmapped-only reads, malformed names/AS tags, score ties,
  interleaved molecules and counted-read mappings.
  **Evidence (2026-10-18):** `...followups/20261018-112322-wave2-g00/cards/b18/`
  (`profile_b181.json`, `contract_b181.json`). The frozen two-phase shape holds
  `{qname: {side: score}}` (184 B per unique usable read) across the whole
  molecule-building pass: 300.3 B/read of per-read state at 200 k unique reads, and the
  peak is `max(scores + molecules, molecules + rows)` — 625 B/read at 1 read per
  molecule, 341 B/read at 8. Frozen (20 cases): secondary hits usable, supplementary/
  unmapped-only leave no molecule, no-`AS` records unusable, malformed names create no
  molecule, ties stay `tie`, order = first *usable* alignment, `reads` = unique usable
  read names, `counted` follows `assigned_weight > 0` and is sticky, `ValueError` in
  stream order, `FileNotFoundError` for the lineage before the SAM is touched, `KeyError`
  without `assigned_weight`.
- [x] B18.2 — Prototype a fused reducer: retain a seen-usable-QNAME structure for
  unique-read counts and per-molecule best-side state; update maxima in stream order.
  Preserve molecule insertion order and the public `best_scores` API. The old code
  parses the entire SAM before parsing molecule names; adversarial mixed faults must
  retain that exception precedence. If not feasible, keep a two-phase disk-backed
  ordered state representation or reject the refactor. Stream output only after
  matching the original parse-before-write boundary; keep gzip writer/header policy.
  **Implemented (2026-10-18).** `_parse_molecules` (seen-QNAME set + `_Molecule` map,
  maxima in stream order) with `_row_of`/`_iter_rows` and `molecule_rows` yielding
  through them; the lineage read stays before the parse, and both writers parse before
  opening the gzip and then stream rows from the molecules. Fused parse alone is a wash
  at 1 read per molecule (0.968×, the seen set then overlaps the molecules), so the
  streaming writer is the other half; materializing rows on the fused parse measures
  0.970× there, i.e. that variant is the trade the stop condition forbids. Peak 0.728×
  (1 read per molecule), 0.641×, 0.571×, 0.526× (8), 0.468× at 1 M reads; wall 1.034×
  (2/5 trials faster) and 0.936× (5/5) in alternating best-of-5 trials without
  tracemalloc. No temp file (gzip `FNAME` carries the output name, so spooling would
  change bytes). Evidence: `cards/b18/{feasibility,decision}.md`, `diff.patch`,
  `memory_b182.json`, `mutation_b182.json`.
- [x] B18.3 — Run `test_molecule_verdict.py`, `test_dsr02_verdicts.py` and evidence
  tests; compare verdict gzip raw bytes, summary, errors and partial files. Benchmark
  high reads-per-molecule and all-unique reads, including disk/state overhead.
  **Evidence (2026-10-18):** 246 passed (`pytest_b183.log`); full suite 5,028 passed /
  19 skipped / 138 deselected (`pytest_full_b182.log`, wave-1 C00 recorded 4,948 + this
  wave's 32 B25 and 48 B18 tests). `exact_b182.json`: 21 frozen cases × 6 comparisons
  (rows, gzip payload, raw gzip bytes, summary record, summary file, file existence) and
  8 `write_from_bam` arms through a stub `samtools` (success and non-zero exit) all
  identical, no partial or `.tmp` file. `mutation_b182.json`: 11 non-equivalent mutants,
  10 caught, 1 documented equivalent control. No disk state; the card's disk-backed
  fallback was not needed.

**Stop:** first-seen order, unique-read count, validation precedence or compression
bytes change. This may reduce per-read dictionaries without making memory constant.
Both hold: identity is exact (`exact_b182.json`) and the peak is reduced, not constant.

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
  Reorder priorities only with recorded evidence. **Discharged 2026-10-18:** the wave-1
  entries below carry packet paths, hashes and gate results; G00 is done and G01/G02/C00
  are partial. The earlier wording here ("G00/G01/G02/C00 remain not started") predated
  the wave-1 runs and was corrected in the `PERF-DOC-01` reconciliation.

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

- 2026-10-18 — C00 wave-1 combined qualification (B24, B23, B04, B06) in worktree
  `vs_work/perf_wave1`, branch `codex/perf-wave1-b24-b23-b04-b06`. Combined source vs the G00
  baseline (as a git checkout, so provenance kind matches): 148-artifact `file_manifest`
  identical, baseline A/A identical. Full suite 4,946 passed / 19 skipped / 2 failed; both
  failures are repo governance pins on `build_reference.py` (artifact inventory, claim registry),
  expected after B24 and re-pinned in follow-up commits. Ruff clean; mypy clean after the B23
  `None` pass-through (the card candidate had one `[index]` error). Benefit, in one session on a
  host with load 8-61: B24 gate yes; B24 preparation supported but with one overlapping outlier
  process; B06 memory yes, time neutral; B23 memory yes for sparse-negative, time inconclusive;
  B04 time inconclusive and +4 MB view cost. No ratio products, no whole-pipeline claim.
  Packet: `analysis/computational_bottlenecks_2026-10-17/followups/20261017-122240-G00/c00/C00_REPORT.md`.
  **Pre-rebase tip suite (added at the `PERF-DOC-01` reconciliation):** after the two governance
  re-pins, the suite at the then-tip of the *pre-rebase* branch was **4,948 passed, 19 skipped,
  138 deselected, exit 0**
  (`c00/pytest-full-final.log`, 843.7 s, sha256
  `05477922db7a5cad626ce68f5628c9c3e8e9eac859dc4ba3b55c48e04288eb78`), with the JUnit receipt in
  `pytest-full-final.log.xml`. The `4,946 passed / 2 failed` numbers above describe the
  pre-re-pin combined state; both failures were stale governance pins on `build_reference.py`
  and are resolved by `6c8cb96` and `f852214` (the re-pin commits as rebased onto `main`;
  the pre-rebase SHAs `5ea27cd`/`faf5667` no longer exist in this history). That tip was replaced
  by the rebase onto `main`; the rebased tip's numbers are in the last entry of this log.

- 2026-10-18 (reconciliation pass, `PERF-DOC-01`) — Corrected this tracker against the evidence
  receipts. No code, golden, comparator, oracle or gate verdict was changed, and no card was newly
  claimed as benefit-qualified. Changes: (1) the header status line no longer reads
  "implementation not started"; (2) the "Pending execution" bullet no longer claims G00/G01/G02/C00
  are not started; (3) the certified final tip-suite result (4,948 passed / 19 skipped /
  138 deselected, exit 0) is recorded with its log path and sha256 alongside the pre-re-pin
  4,946/2-failed state; (4) `C00_REPORT.md` §3b, which was an empty placeholder, is filled from
  that log. Scope note: B08 was decided in the separate wave-2 packet — **deferred**, 3.8–7.3 %
  median paired-reporting reduction against the 10 % default gate with no memory substitute and an
  external-host-set path that never calls the helper
  (`followups/20261018-112322-wave2-g00/cards/b08/decision.md`); B10's failure matrix
  (`cards/b10/failure_matrix_b101.json`) and wave 2 (B25, B27, B18) are still running there. This
  pass deliberately does **not** rewrite those card bodies or the B08/B10 sections: the wave-2
  branch owns them and must rebase onto this commit. Unchanged by design: G02 stays partial (one
  loaded-host session, component fixtures, no whole-run measurement), B08.2/B08.3 stay unchecked
  behind the failed gate, and every wave-1 benefit statement stays component-level.

- 2026-10-18 — **Rebase of wave 1 onto `main` (`6910c2d`), and the merge fix it forced.** The
  branch's original base was 188 files behind. The rebase took one real source conflict in
  `build_reference.py` (AUD3-05's near-duplicate screen × B24's split helpers; the resolution keeps
  AUD3-05's signature and screen call and delegates to B24's helpers in AUD3-05's order — structure
  → screen → low-complexity) plus four tracker/governance conflicts. Both governance digests were
  recomputed at the merged bytes, not side-picked: `build_reference.py` →
  `88934e091793fe48a7dffcd8b8f1f3315eaeeaf50bfd52ab840263d42ed66c9a` pinned at `6c8cb96`, and the
  `claims/registry.json` inventory row →
  `3a2a0b1620013e918d8a87087e0715f04a07d470e54337e5f69a6b337b18d39e` pinned at `f852214`. The
  pre-rebase tip is preserved locally as `codex/perf-wave1-prerebase` (`dec9885`).

  The rebase also surfaced a **merge-induced defect**, which is why the branch was not landed as
  rebased. `main`'s Q123 per-virus positive-cell count (`329b3dc`) indexes the layer that
  B05/B07/B20's budget-gated view had already rebound to a narrow column span, and it used an
  absolute position for it. Each feature was correct on its own branch; combined they disagree
  about the index space. It is reachable in production: `detection.py` groups every viral gene
  present and `build_reference.py` concatenates the host GTF first, so viral genes are always a
  trailing block. A 2,040-column matrix (2,000 host + 40 viral) raised
  `IndexError: index (2019) out of range` where `main` returns a complete frame. Every existing
  grid missed it because its fixtures use leading genes, where the view coincides with the
  absolute indices. Fixed in `5d7de93`: resolve the index through `gene_positions` (the mapping the
  sibling per-gene access already uses) and accumulate the mask per integer column through
  `_gene_vector`, which keeps the cache budget's no-submatrix promise and needs no test edit.
  Verified: 90/90 exact against main's implementation over the method × layout × dtype grid, and
  identical on the production-shaped input and on the non-leading layouts that previously raised.

  The frozen evidence oracle was re-based in `b6bec78` (user-approved). The two
  `tests/fixtures/performance_original.py` blocks that `main` changed are re-pinned to `main`'s own
  source; the other six stay byte-identical to `5e551c9`; both re-based blocks match `main`'s
  blocks byte-for-byte; and minimality was proved by reverse-swapping the two blocks and recovering
  the original file byte-for-byte. No assertion, comparator, tolerance or gate verdict changed, and
  nothing was derived from the optimized code.

  **Rebased tip suite: 5,512 passed, 41 skipped, 139 deselected, 6 failed** (541.3 s,
  `c00/pytest-full-rebased.log`, sha256
  `90151e142d1404de69a4a22a630a059ce48d273471573cc24c80b119f15ed22d`). The six failures are
  `tests/test_release_static.py::test_weekly_audit_keeps_exit_code_and_log[*]`, pre-existing and
  environmental: the same six fail on a pure `origin/main` checkout because the test shells out to
  `python`, which is not on `PATH` here (`bash: line 3: python: command not found`). For contrast,
  the rebased tip *before* the two fixes was 98 failed / 5,420 passed — 92 of them the fixture
  vintage recorded below, plus these same six.

  Note for the wave-2 session (B25/B27/B18): `main` changed `molecule_verdict.py` by +108 lines of
  behavioural code (EVID-UMI-01/Q137 and AUD3-04: `matched_bases`, `_read_rank`,
  `deduplicate_competitive_alignments`). B18's frozen oracle is pinned at `faf5667`, which predates
  all of it, and the line numbers that card cites (`:61,96,101,153`) have shifted; re-check the
  oracle target against the post-merge file before trusting it — the same class of issue this entry
  records for B24.

  Tags: PERF-EXEC-02, PERF-DOC-01, rebase, byte-identity, LOCK-01, wave-2
