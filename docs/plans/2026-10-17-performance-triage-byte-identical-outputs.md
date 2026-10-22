# Performance triage and byte-identical output plan

> Path placeholders: `$SITE_ROOT` = site work root (institutional filesystem paths, replaced by placeholders).

Date: 2026-10-17. Project: ViralScan.
Status: first-wave implementation complete; E01–E05 passed; E06 assessed/deferred.
Original source baseline: `5e551c96831c0c152da9fa450af6a75b1266e610`; the 56-file
survey inventory matched before editing. Captured source: `execution/original/` in
the survey packet. Final results: `execution/REPORT.md` and `execution/validation.json`.

## Objective and decision

Implement the smallest high-benefit changes first, preserving the exact analytical
outputs of the current implementation. Start with **B07 count formatting**, followed
by **B20 binary-label UMAP enrichment** and **B05 local column-access optimization**.
Keep all counting policies, thresholds, diagnostics, random draws, native-tool
versions, configured threads, output schemas, row order, and serializers unchanged.

The user authorized local implementation and validation of E01–E05 on 2026-10-17,
followed by assessment of B08/B24. The original triage consisted only of review and
experiments; execution evidence is recorded separately under `execution/` below.
No SLURM submissions, large downloads, retained-data deletion, commits, push, or
publication are authorized.

The survey and receipts are under
`analysis/computational_bottlenecks_2026-10-17/`. Its previous DataFrame/tolerance
checks establish useful evidence but are not a general proof of raw byte identity.
The new `byte_identity_probe.py` specifically compares types, numeric bits, warnings,
errors, and actual serialized output where stated.

## Priority: benefit, effort, and byte risk

Effort estimates include a focused implementation and meaningful regression tests
in a prepared environment. They are planning estimates, not elapsed-time receipts.
Speedups are measured for individual synthetic workloads, not the complete pipeline.
This table retains the pre-implementation triage measurements; accepted production
code measurements and the narrower B05 implementation are recorded in E02–E05 and
the execution report. The experimental CSC ratio is not the implemented B05 speedup.

| Rank | Candidate / applicability | Benefit evidence | Estimated effort | Byte-preservation strategy / decision |
|---:|---|---|---|---|
| 1 | **B07: `_count_value()`**, default reporting | Full `compute_stats`: 2.525 → 0.439 s, **5.75×**, guarded variant | 1–2 h | Preserve arithmetic/tolerance and Python result types; guarded scalar fast path with original fallback. Strict renderer probe passed. First change. |
| 2 | **B20: `viral_neighbor_enrichment()`**, UMAP enabled | 2.972 → 0.130 s, **22.9×** on 2,000 cells / 300 positives | 1–3 h | Vectorize the inner binary-label means only; preserve RNG calls, neighbor queries, positive-cell set, permutation order, and final reductions. 24 exact-bit cases passed. |
| 3 | **B05: `summarize_multimap_evidence()`**, default multimapping report | 1.796 → 0.316 s including CSC conversion, **5.69×** | 2–4 h | Optimize local reads, leaving stored layers and existing per-column computations intact. Exact CSV checks passed on two fixtures; expand adversarial and peak-memory checks before accepting. |
| 4 | **B24: repeated reference classification**, reference builds only | Three scans per record; potential near 3× for classifier work, **not yet a measured rewritten gate** | 3–5 h | Reuse identical integer counters within one invocation; preserve structural validation, failure type/message, offender ordering, logging, and reference bytes. Next build-focused change. |
| 5 | **B08: repeated denominator**, default reporting | 2.635 → 2.541 s, about **3.6%** in that fixture | 0.5–1.5 h | Hoist the same read-only calculation; preserve all host-set/None/empty-set branches. Small, safe cleanup; lower benefit than ranks 1–3. |

B10, skipping the disabled UMAP load, looks like a tiny code edit and could avoid a
large read. It is conditional rather than first in this list: the old code also
checks the schema and surfaces missing/corrupt H5AD or found-gene inputs. An early
return must not silently weaken those failure contracts. No runtime speedup for
this change has been measured. Resolve the validation contract before implementing.

Rank 2 is conditional on UMAP use; for workloads without UMAP, rank 3 follows rank 1.
Ranking is not a claim that UMAP dominates every run. A component speedup `s` affecting
fraction `f` of a run gives approximately `1 / ((1-f) + f/s)` total speedup. For
example, a 5.75× improvement to 10% of runtime improves the whole run only about 1.09×.
Measure stage shares rather than adding the survey's ratios.

The B07 timing above is the fresh guarded-variant measurement in
`receipts/guarded_formatting_timing.json`: one warm-up and three trials per variant,
20,000 cells, 2,000 genes, and 38,187 per-cell rows. Exact statistics, table values,
and serialized per-cell bytes passed. The survey's earlier 5.56× result used the
unguarded experiment; it is not substituted as evidence for the guarded proposal.

## Review findings incorporated into this plan

| Review finding | Resolution |
|---|---|
| Numeric closeness and ordinary DataFrame checks do not establish bytes | Acceptance uses raw-file comparisons plus exact dtype/bit checks; no `allclose` acceptance gate. |
| Sparse summation changes addition order | Demonstrated counterexamples; direct B04/B06 reduction rewrites moved out of the first wave. |
| Unrestricted scalar replacement could change large-value/error behavior | Proposed B07 fast path is guarded; original implementation remains the fallback. Exceptions and warnings are compared. |
| Vectorizing UMAP could change the null calculation or random stream | Retain the existing rule, all neighbor queries, and RNG call sequence; fast path only for supported binary arrays. |
| Faster column access can increase RAM or change sparse duplicate behavior | Require canonical/noncanonical, explicit-zero, dtype, raw-byte, and peak-memory tests; keep changes local to reporting. |
| Disabled UMAP currently validates inputs even when plotting is off | B10 remains conditional until validation and failure behavior are specified. |
| Whole run directories already contain variable and intentionally changing bytes | Inventory baseline-versus-baseline variability; distinguish exact artifacts from explicitly identified operational/provenance fields. Never claim whole-directory identity while excluding fields. |
| Goldens generated by the candidate would be circular evidence | Freeze reference outputs from the captured original source before edits, hash them, and never update them to make a candidate pass. |

The executed artifact contract is fixed in `execution/artifact_contract.json`.
The only allowed detection-stage operational span is the truthful 64-character
build hash in reference provenance; every other byte is compared unchanged.
B10's validation behavior remains unresolved and outside this implementation.

## What byte-identical means

The strict requirement is **same relative file set, same file lengths, and exactly
the same bytes** for each artifact designated byte-stable. SHA-256 is a receipt and
fast screen; the local comparator should also compare the files themselves and
report the first differing byte. Do not sort tables, round numbers, replace NaNs,
drop columns, normalize whitespace, decompress files, or rewrite HDF5 before this
check. Those operations could conceal regressions. Decoded comparisons are useful
diagnostics after a raw comparison fails; they cannot turn that failure into a pass.

Exact equality of floating values alone is also insufficient: distinguish `int` from
`float`, preserve dtypes and `-0.0`, and compare floating bit patterns. Rows, columns,
gene/barcode order, index types, explicit zeros, sparse indices/indptr, missing-value
representations, quoting, newline conventions, and JSON field/array order can affect
bytes and consumers.

### A complete run directory cannot currently be promised byte-identical

The actual source writes runtime into `summary.txt`
(`menu._write_sample_summary`, line 1070), generates a current date for `report.html`
(`detection.generate_html_report`, line 1945), and records source hashes/Git state
in run/evidence/host-response provenance. A successful performance change should
change measured runtime and must honestly change its recorded software identity.
Gzip headers can change with creation time even when payload bytes are identical;
the probe demonstrates this. Plotly-generated IDs, directory enumeration, native
tool headers/ties, and HDF5 serialization need baseline repeatability checks rather
than assumed stability.

Therefore use two explicitly different checks:

1. **Controlled regression fixtures:** same paths, explicit report date, same pinned
   dependencies, seeds, locale, native tools, thread counts, and controlled clocks or
   renderer IDs where needed. These controls are test-only and applied identically
   to original and candidate. Compare the resulting original files byte-for-byte.
2. **Honest execution receipts:** retain actual source hashes, runtime, timestamps,
   commands, and paths. Record any operational differences separately with precise
   field locations and expected values. Never forge production provenance or freeze
   real runtime to make files match.

Proposed operational exceptions must be exact: for example, the one appended
`Runtime: <number> seconds.` line in `summary.txt`, the specific rendered run date,
or `run_manifest.json`'s `software_identity` fields. Preserve every other byte/field,
validate that actual new hashes identify the candidate, and report the raw-file
hashes as different. Do not allow blanket `*.json`, logs, image, H5AD, gzip, or
timestamp-regex exclusions. Enumerate actual pointers/fields after baseline capture;
other artifacts default to strict comparison.

**If the requirement is literally every file in the complete run directory, with
no operational exceptions, no implementation may be labelled compliant with that
requirement under the current output contract.** Test-fixture byte identity and
identical analytical payloads are narrower claims and must be named as such.

## Pre-implementation strict evidence

The receipt `analysis/computational_bottlenecks_2026-10-17/receipts/byte_identity_triage.json`
and its SHA-256-pinned probe record:

- **40,172 scalar cases:** same Python return type, representation, float bits where
  applicable, exception type/message, and warnings. Inputs include signed zero,
  subnormal values, tolerance boundaries and their neighboring floats, fractional
  counts, large integers/floats, NaN/inf, and four `ndigits` values.
- **Actual reporting functions:** guarded formatting produced byte-identical
  `viral_summary.tsv`, `per_cell_viral.tsv`, `host_homology_loci.tsv`, and `report.html`
  on integer and fractional fixtures. The original and candidate render into the
  same temporary path with an explicit report date; this does not establish that
  default timestamps/paths repeat across independent production runs.
- **24 UMAP cases:** the three returned numeric values and return types match exactly
  for bool/int64 labels, two seeds, single/no/all/partial positives, neighbor clamping,
  and up to the full 1,000 permutations. This is fixture evidence, not a universal
  API or statistical-validity proof.
- **B05:** current CSR and CSC summary outputs match with `check_exact=True` and
  identical serialized CSV/TSV bytes on integer and fractional fixtures. Noncanonical
  duplicates, boundary-focused values, and retained production matrices remain gates.
- **Reduction and compression hazards:** positive-count dense and sparse reductions
  produce different bits; equal gzip payloads can produce different compressed
  bytes. Neither tolerance checks nor decompression is an acceptable strict-byte gate.

For example, one row with `[1.0] + [1e-16] * 32` gives
`1.0000000000000029` through the current dense helper and
`1.000000000000003` through the proposed direct sparse sum. All entries are positive;
the disagreement does not depend on negative values or integer overflow.

Finite fixtures do not prove equivalence for every representable input. The source
argument and guarded fallback constrain the domain where a fast path can be used.

## Candidate-specific implementation constraints

### B07: keep the number-formatting policy unchanged

Preserve the original sequence `float(value)`, `round(value)`, tolerance check, and
final `int`/`round(value, ndigits)` return. For the proposed scalar fast path, use the
same asymmetric NumPy default inequality, not `math.isclose` or a new tolerance.
Restrict it to rounded integers exactly representable in float64 (at most `2**53`
in magnitude); call the original NumPy check outside that domain. Nonfinite input
must continue failing at the original conversion/rounding point.

Keep native Python float operations and pinned dependencies; do not introduce
float32, fused arithmetic, approximate comparisons, or an alternate rounding rule.
Tests must cover value/type/bytes and warnings/errors. `_count_value` is used beyond
the timed per-cell rows, so include headline and aggregate serializers as well.

### B20: change the loop mechanics, not the statistic

Keep `np.random.default_rng(random_state)`, the full `.kneighbors(coords)` call,
`k` clamping, self-neighbor removal, original positive-cell indices, and exactly
one `rng.permutation(labels)` call in the same order per permutation. Retain the
outer one-dimensional mean and the original `>= observed` comparison and denominator.

For binary bool/integer labels the inner sums are small exact integers before the
same floating division; a vectorized row mean is a promising way to remove scalar
dispatch. Keep the old path for nonbinary/nonstandard inputs unless equivalence is
established separately. Add duplicated coordinates/ties, zero/negative permutation
budgets as currently handled, empty inputs, multiple integer dtypes, and explicit
Generator-state tests to the existing 24 cases. Do not change the original statistical
rule or query only the positive coordinates as part of this optimization.

### B05: local access changes only in the first version

Convert or cache local reporting views without assigning back to `adata.X` or
`adata.layers`, rewriting H5AD, changing gene iteration order, or batching arithmetic
reductions yet. Keep `_sum_gene`, `_gene_vector`, rounding, tier logic, and output
construction unchanged where possible. Preserve absent-layer fallback aliases.

Test empty matrices, dense/CSR/CSC, float32/float64, explicit zeros, unsorted indices,
duplicate coordinates, missing layers, repeated gene IDs in groups, unknown genes,
near-rounding-boundary values, and all multimap methods. If conversion changes bytes
on an allowed representation, retain the original path for that representation.
Measure peak memory: retaining five full CSC copies can defeat the speed benefit on
large runs. Prefer bounded layer/view caching where it can preserve existing
per-column computation; do not claim that every batched sparse reduction is equivalent.

### B24 and B08: reuse work without reordering decisions

B24 should reuse exactly the existing integer classifier results, not rewrite the
classifier, change masking, reduce thresholds, or change native reference builders.
Retain public validation behavior, structural checks, offender/tie order, log messages,
and error types/messages. Avoid a process-global sequence cache, which can retain
large references or reuse stale state. Test valid and invalid FASTA fixtures and
require identical prepared FASTA/GTF/report bytes before measuring builds.

B08 should hoist only the invariant comparable-cell calculation. Cover no host sets,
provided host sets, `None`/empty comparable sets, no host genes, no called cells,
zero viruses, and candidates failing detection. Do not compute a previously skipped
helper before an earlier validation error, or mutate input masks. Its small benefit
does not justify a broad rewrite of `compute_stats`.

## Full candidate disposition

Every original survey candidate is accounted for below. E02–E05 record the bounded
implementation and qualification of the first wave; other dispositions remain triage.

| ID | Disposition | Reason under the byte-identity requirement |
|---|---|---|
| B01 allocation/state | Later, substantial refactor | Aggregation/chunking can change sparse accumulation and ordering; preserve all method/lineage contracts. |
| B02 cell EM | Later, high numeric risk | Workspace/prior/convergence changes can alter fitted weights; no exact-proof implementation yet. |
| B03 assignment replay | Later | Gzip bytes, row order and lineage must remain exact; do not disable requested output. |
| B04 gene reductions | Defer direct batching | A reduction-only timing is not a completed speedup; changed float addition order is unsafe without a constrained exact path. |
| B05 evidence column access | Implemented, E04 | Narrow CSR subset cache; exact outputs, bounded cache and measured memory tradeoff. |
| B06 dense viral totals | Defer direct sparse sum | Large memory benefit, but an executed positive-count counterexample disproves universal bit identity. Explore row-block dense reduction preserving dtype, selected-column order, and reduction width, then test exact bits across block sizes. |
| B07 count formatting | Implemented, E02 | Guarded scalar path; strict typed/rendered output checks and measured benefit. |
| B08 denominator reuse | Small follow-up, rank 5 | Same read-only arithmetic, modest measured gain. |
| B09 alias storage/read | Split: selective reads later; storage migration deferred | Removing/deduplicating persisted layers necessarily changes H5AD bytes/schema. Reader-only changes can preserve artifacts but need validation coverage. |
| B10 disabled UMAP load | Conditional quick win | Skipping input checks can alter failure behavior; specify and preserve the validation contract first. |
| B11 thread wiring | Separate resource-correctness work | Native execution/EC/tie order can change; speedup unmeasured and exact native artifact bytes not established. |
| B12 emptyDrops | Preserve required analysis | Do not reduce iterations or substitute caller; denominator and RNG behavior are scientific contracts. |
| B13 read artefact filter | Later, moderate/high effort | Classifier/parallelism must retain mate/order checks and identical serialization/compression. |
| B14 host-filter I/O | Later | Share scans while preserving exact validation/lineage/error behavior; native run still needs qualification. |
| B15 shared replay | Later, high impact but high effort | Cache identity, EC maps, correction/strand, target dispatch and record ordering require native evidence. |
| B16 piped alignment | Later | Child failure propagation, headers, tie order and raw BAM bytes must be checked. |
| B17 combined QC passes | Later | Changed accumulation order and raw/dedup policies need exact comparisons. |
| B18 bounded verdict reducer | Later | Group/name sorting can change existing first-seen output order; preserve it explicitly. |
| B19 BLAST cache | Later | Benefits unmeasured; require immutable database identity, subject order and tie/search policy. |
| B20 UMAP enrichment | Implemented, E03 | Binary vectorization; exact bits/RNG checks; original fallback retained elsewhere. |
| B21 host feature storage | Later, major memory work | Sparse solver/dtype/feature changes can change floating trajectories and predictions. |
| B22 fit parallelism | Later | Preserve fold membership, independent RNG state, fit/reduction order and numerical environment. |
| B23 programme empty-row work | Conditional small follow-up | Promising simple pruning, but gain unmeasured; marker-negative state and summary semantics must be exact. |
| B24 reference counter reuse | Next build-focused wave, rank 4 | Integer counter reuse is safer than changing the classifier; potential gain applies to builds only. |
| B25 streaming reference text | Later | Primarily a memory opportunity; preserve newline, uppercase, headers, record order and wrapping exactly. |
| B26 hashing/provenance | Preserve checks; same-operation reuse only later | Never skip integrity checks or replace content hashes with size/mtime; real new software identity must remain honest. |
| B27 whitelist scanning | Low priority | Already bounds reads and streams lists; speed/RAM tradeoff needs evidence. |
| B28 retrieval/cache | Low priority | Throttle is deliberate; duplicate digest reuse may help but native/network benefit is unmeasured. |
| B29 dual BAM decode | Later | Compatible accumulators may share a pass; preserve read dedup and unique-virus molecule rules. |
| B30 native build workers | Separate resource-correctness work | Thread choices can change native ordering/clustering artifacts; not an initial byte-safe performance change. |

## Progress and execution gates

### [x] E00 — Review and exactness triage

- Status: complete.
- Outcome: reviewed ranking, all 30 dispositions, explicit operational-byte limitation,
  executed strict probe, and this tracker.
- Evidence: `byte_identity_triage.json`; original survey source hashes verified.
- Scope at E00: fixture experiments only; production implementation had not started.

### [x] E01 — Freeze independent reference and artifact contract

- Status: passed. Dependency: E00. Original A/A matched all 689 captured files;
  40,192 scalar cases, 103 UMAP cases, 255 evidence cases, and 20 reporting cases.
  Existing focused suite: 675 passed, two existing warnings. Source inventory of
  56 runtime files matched; 76 package/resource files captured. Goldens read-only.
  Evidence: `execution/source_baseline.json`, `execution/artifact_contract.json`,
  `execution/runs/original-{a,b}/exact.json`, and `execution/original_tests.xml`.
- Capture the original source bytes, dependencies/tool binaries, locale/timezone,
  thread variables, seeds/configuration, immutable input hashes, and serializer
  settings. Use isolated copies; do not reset/stash the shared user checkout.
- Run original-versus-original before editing. Inventory the relative file set,
  lengths and hashes; identify every nonrepeatable artifact and exact changing field.
  Do not silently exempt a file that fails this test.
- Freeze analytical goldens from the original source. Keep originals immutable and
  hash them. Recompute neither goldens nor expected hashes using the candidate.
- Define the controlled-rendering test settings and the precise operational-field
  policy. If a strict whole-tree policy is required without exceptions, stop short
  of claiming compliance; current runtime/provenance serialization prevents it.
- Evidence saved: original A/A receipt, golden manifests, environment listing,
  artifact/field contract. No SLURM submissions were made.

### [x] E02 — Implement B07 in isolation

- Status: passed. Dependency: E01 passed. All 689 files byte-identical to original;
  123 focused exactness tests passed. Three-trial `compute_stats` median
  2.4868 → 0.4550 s (5.47×); negative fixture 0.2114 → 0.2024 s.
  Evidence: `execution/runs/b07/{exact,timing}.json`, `execution/b07_tests.xml`.
- Affected areas: `scripts/detection.py`, its actual helper callers/serializers, and
  focused tests for formatting and complete reporting artifacts.
- Apply only the guarded fast path; keep fallback and serializers unchanged.
- Require all strict scalar cases, aggregate/headline cases, and actual deterministic
  reporting bytes to match; compare unchanged H5AD/count artifacts by hash.
- Measure three-trial timings separately from equality tests. Require a repeatable
  benefit on representative formatting-heavy and sparse-negative cases, not a noisy
  single run. Revert the change if exactness or existing failure tests fail.
- Evidence saved: baseline/candidate raw manifests, test receipt and timings.

### [x] E03 — Implement B20 in isolation

- Status: passed. Dependencies: E01–E02 passed. All 689 files match original;
  141 focused tests passed, including RNG/ties/fallbacks. Three-trial enrichment
  median 2.8694 → 0.1408 s (20.39×). Only binary bool/integer arrays use the fast path.
  Evidence: `execution/runs/b20/{exact,timing}.json`, `execution/b20_tests.xml`.
- Affected areas: `scripts/umap.py`, its caller, and focused UMAP tests.
- Preserve the complete statistic and default RNG/neighbor behavior; optimize only
  the proven binary-label path. Keep the old path for unsupported inputs.
- Require exact return types/bits, p-values, labels, and controlled HTML/report bytes;
  also assert identical RNG consumption and unchanged input arrays.
- Evidence saved: expanded edge/seed/tie checks and a separate timing receipt.

### [x] E04 — Implement B05 access optimization in isolation

- Status: passed. Dependencies: E01–E03 passed. A narrower CSR subset cache keeps
  the original per-column CSR calls/storage format. All 689 files match original;
  708 focused tests passed. Three-trial median 1.8594 → 1.0497 s (1.77×).
  Original A/A and B05 match on a retained 1,024-row em-cell slice, 256 genes,
  all six method labels over fixed em-cell layers, including 2,362 fractional
  corrected entries. No allocation refits or full reload.
  Cache limited to 32 MiB CSR buffers, with unchanged original path above that
  budget or for CSC/dense/noncanonical inputs. 100,000-cell memory probe passed.
  Evidence: `execution/runs/b05/{exact,timing,memory}.json`, `execution/b05_tests.xml`,
  and `execution/retained/{original,repeat,b05}/receipt.json`.
- Affected areas: `multimapping.summarize_multimap_evidence`, `_sum_gene`/
  `_gene_vector` callers if necessary, and evidence-summary tests.
- Keep arithmetic and iteration order initially unchanged; cache only local views.
  Do not persist a new storage format or modify the original layer buffers.
- Require exact table types/bits/bytes across the representation/edge fixture matrix,
  unchanged evidence tiers, hashes of every input sparse buffer, and representative
  retained-output comparisons when bounded local execution is feasible.
- Record peak RSS as well as wall time for the function. Do not accept a speed gain
  that causes an unacceptable peak-memory regression on a large retained matrix.
- Evidence saved: strict artifact comparisons, input immutability and memory/time.

### [x] E05 — Qualify the combined first wave

- Status: passed. E02–E04 individually passed. Combined source matches all 689
  original fixture files and 12 retained-slice files. Actual detection entry point
  covers seven success/failure scenarios, 77 files, with original A/A repeatability.
  Four provenance files differ only in the truthful 64-character build hash;
  every other byte matches. 1,055 tests passed; one R/DropletUtils test skipped.
  Independent portable original helpers/schema frozen; 125 exactness tests
  reverified afterward. No production-source changes after the combined snapshot.
  Combined Python workload median 6.9309 → 1.6463 s (4.21×), not a whole-run gain.
  Evidence: `execution/validation.json`, `execution/combined_tests.xml`,
  `execution/runs/combined/`, `execution/detection_stage/combined/receipt.json`.
- Compare each isolated candidate to the same original reference, then compare the
  combined version. One optimization must not hide another's output drift.
- Check file presence, raw bytes, arrays/types/order, errors/warnings, original
  package tests relevant to touched code, and the exact operational-metadata policy.
- Include positive/negative, fractional allocations, all supported multimap methods,
  threshold-boundary cases, host-set branches, and enabled/disabled UMAP paths.
- Report isolated and combined timing improvements without multiplying speedup
  ratios. Document platform/thread scope and remaining native/end-to-end gaps.
- Evidence saved: combined manifest/diff receipts, tests and timings; final summary
  in `execution/REPORT.md`. No release or scientific claim follows from this gate.

### [x] E06 — Assess optional follow-ups after the first wave

- Status: assessed; B08/B24 implementation deferred. Dependency: E05 passed.
  B08's post-wave memoization pilot matched the checked fixture and reduced helper
  calculations 20 → 1, but paired reporting medians 1.7338 → 1.9760 s did not show
  a consistent gain. Full host/empty/error branches remain unqualified for that change.
  B24 is build-only and has no equivalent rewritten gate/diagnostic/reference-byte
  qualification. Evidence: `execution/followup_assessment.json` and execution report.
  The remaining bullets are future work, not implemented changes.
- B24: share existing counters, compare reference bytes and all failure diagnostics.
- B08: hoist unchanged denominator work and compare every branch and output table.
- B10/B23: proceed only after their validation/state contracts and benefit are shown.
- B04/B06 or resource/native refactors require a new exactness argument and their
  own implementation plan; the original survey's tolerance checks are insufficient.

## Artifact comparison matrix

| Artifact / interface | Required comparison |
|---|---|
| Detection/evidence/sensitivity/called-cell TSVs and `log/found_genes.txt` | Raw bytes, file set/order, plus semantic diagnostics on failure; no sorting or rounded comparison to pass. |
| `_count_value`, UMAP statistics, intermediate count arrays | Exact type/dtype/shape and numeric bytes, preserving signed zero and missing-value behavior; compare warnings/errors as well. |
| Input and untouched output H5AD/count files | Require unchanged raw file hashes/bytes. If a stage writes H5AD, first establish A/A serialization stability; compare all sparse buffers and metadata exactly as diagnostics, not as a substitute for a required raw-byte gate. |
| `report.html` and plot artifacts | Raw bytes under declared fixture date/path/ID/plot settings. Inventory real-run variability before claiming strict production-file identity. |
| Gzip lineage/assignment files | Raw compressed bytes required where designated stable. Decompressed-byte equality is a separate diagnostic claim. Do not change compression backend/level/headers implicitly. |
| `summary.txt` | Analytical text exact; the actual appended runtime is intrinsically different. Any runtime-field exception must be precise and reported, not silently stripped. |
| Run/evidence/host-response manifests | Analytical fields and referenced analytical artifact hashes exact; actual new source/tool identities truthful. Enumerate permitted operational pointers; reject other changes. |
| Failed/empty runs | Same exception class/message, warning category/message, exit status, required status artifacts and absence/preservation of previous outputs. No validation removed for speed. |

A candidate fails if an analytical file appears/disappears, a stable byte differs,
an input is mutated, a decision/tier changes, or a failure becomes success. Do not
update a golden, loosen a tolerance, or add an exclusion in response to a failure.
Investigate the difference and keep the original code path if equivalence cannot be
established. Rollback is the isolated source change; never replace retained reference
outputs to make a candidate look compatible.

## Decision/evidence log and reproduction

- 2026-10-17 — Ranked by measured benefit, applicability, edit/test effort and exactness
  risk. No production optimization started. All B01–B30 candidates assigned a disposition.
- 2026-10-17 — Strict formatter/rendering, binary UMAP, and CSR/CSC summary probes
  passed on the saved fixtures. Sparse-sum and gzip-header counterexamples recorded.
- 2026-10-17, execution — E01–E05 passed with independently frozen references,
  isolated/combined checks, retained data and actual detection-stage files.
  E06 assessed; B08/B24 deferred. Source-only changes are the three named functions.
- Remaining limits — Native/full FASTQ-to-report benefit, production peak RSS,
  cross-version/platform equivalence and universal byte identity remain unqualified.
  B10's failure semantics are unchanged and deferred.

Current execution verification (does not regenerate goldens):

```bash
$SITE_ROOT/conda/envs/viralscan_test_full/bin/python analysis/computational_bottlenecks_2026-10-17/execution/verify.py
```

The historical triage probe below expects the pre-edit source inventory. It now
intentionally rejects the optimized live checkout. Do not reset source or suppress
its hash guard; use the execution verifier and captured original instead. Historical
command, for a checkout containing the original source and pinned environment:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_CACHE_DIR=/tmp/viralscan-survey-numba MPLCONFIGDIR=/tmp/viralscan-survey-mpl $SITE_ROOT/conda/envs/viralscan_test_full/bin/python analysis/computational_bottlenecks_2026-10-17/byte_identity_probe.py
```

The probe verifies the package source hash inventory before execution and records
its own SHA-256. It creates only temporary synthetic files, does not modify package
source, and runs no native aligner, R caller, network request, or scheduler job.
