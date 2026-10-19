# Performance code qualification and remaining gates

This record supplements the frozen historical receipts. PRs #73/#74 merged as
`5e05a27`/`6e70eda`. Remaining-wave baseline `d056986` has the same tree as
`6e70eda`, retained by `provenance/perf-wave2-2026-10-18`. The wave-1 provenance tag
also remains required for fresh governance checks; neither tag is a release tag.

No scientific thresholds, model/fold/iteration settings, search policy, effective
workers or retained datasets are changed. Local fixtures do not establish a
production speedup. No SLURM jobs, downloads or retained-data deletion are authorized.

## Additional accepted code

| Card | Change | Evidence at current stage |
|---|---|---|
| B21 | Bounded dense-then-cast CSR host features | 103 strict +435 caller tests; two seven-pair conversion/cohort sessions; traced peak0.450/0.452x; combined checks below |
| B13 | Invocation-local classifier cache, 256 entries/512 key characters | 47 focused tests; duplicate stage0.0527/0.0525x, reagent0.1356/0.1331x; unique1.0232/1.0173x; combined checks below |
| B28 | One parse of immutable GenBank text, preserving grouped emission | 96 focused tests; full offline stale-GTF preparation0.8964/0.8863x; mutable-path integrity reads retained; combined checks below |
| B16 | Sequential anonymous SAM file between unchanged native commands | Native raw BAM/BAI A/A and A/B pass; 180 focused tests; two independent sessions wall1.0081/1.0016x, treeRSS0.5106/0.4999x; combined checks below |

B21 applies to the full-gene `use_hvg=False` canonical float64 CSR path.
AnnData HVG CSR views retain the original path; no default-HVG memory benefit is
claimed. B21 retains original paths for CSC, other dtypes, noncanonical sparse matrices,
nonfinite/out-of-range data and observable underflow policies. No sparse-model or
eligibility shortcut. Real cell/group/strata fitted fixtures produce27/30 raw equal
files; three manifests differ only in executing-source values, checked as literal
value spans. Baseline A/A is raw identical. Full raw hashes and mutation checks are
retained; this is not a whole-directory byte-identity claim.

B13 preserves raw case, exact barcode/UMI span, per-record validation, reason order,
headers/qualities, gzip/output order and failure prefixes. High-diversity traced
allocations rise9.3% within the explicit cache bounds; long reads use the baseline
classifier. Benefits depend on repeated sequences; no production duplicate-rate
assumption is made.

B28 preserves CDS-first output, separate group ordinals, token collisions, interval
order, warnings and exceptions. Every new cache invocation still revalidates disk.
Warm-cache timing controls are1.0243/0.9697x. No remote-service improvement is claimed.

B16 is a memory candidate, not a runtime improvement claim. SAM disk grows with
alignment text size (28,626,864 bytes for100,000 fixture reads). The producer finishes
before sort; errors include full native stderr. Earlier timing runs instrumented
stdout hashing and are excluded from benefit claims. Dedicated plain timing runs
and separate sampled parent/child memory runs replace them.

## Other remaining cards

| Card | Disposition | Evidence / remaining work |
|---|---|---|
| B01 | Deferred, benefit failure | Full builder wall1.0209x/peak0.9241x; one family, seven pairs; stopped before broader qualification |
| B02 | Deferred, benefit failure | 691 tests; seven paired complete allocator trials show6.9% wall gain, peak1.001x; second session/scratch/trajectory qualification stopped |
| B03 | Deferred, failure contract | Fresh Path reread rejects changed malformed BUS; spool reuse succeeds instead |
| B09 | Deferred, validation contract | Backed read accepts X encoding version rejected by eager read |
| B12 | Deferred, no material adapter opportunity on fixture | Real R4_51 A/A TSV identical; source matrix74736bytes vs~770MB process peak;10000 draws retained |
| B14 | Deferred, failure contract | Retained scan can fail after successful count and aggregate-audit publication |
| B15 | Deferred, failure contract | Later target quantification can independently fail; successful earlier replay does not remove that obligation |
| B17 | Deferred, failure contract | Removing repeated view/depth calls suppresses their independent failures |
| B19 | Deferred, no repeated workload | Current multi-target extraction disables BLAST; direct evidence builds at most one database per invocation |
| B22 | Deferred, insufficient affected-stage share |102 real logistic fits; repeated design preparation0.696% of bounded profiled run; no fits changed |
| B26 | Deferred, integrity contract | Same size/restored mtime, symlink retarget and directory changes alter fresh digests |
| B29 | Deferred, failure contract | Second decode failure is independently observable before final statistics/marker |
| B11/B30 | Separate authorization gate | Concrete resource proposal in2026-10-19-performance-resource-policy.md; no worker/environment/resource edit |

B08 (benefit), B10 (validation) and B27 (negative I/O result) retain their earlier
dispositions. A rejected early gate stops further prototype/timing work; its broad
fixture/edit/qualification checkboxes remain unchecked. Deferred does not mean
implemented or universally impossible to improve in a future contract/workload.

## Combined qualification

The source freeze binds326 source/test files. All five frozen original function
bodies match `d056986` both as AST and raw function source. The B21 oracle formatting
change adds only three preamble blank lines. No original function, golden or
comparison tolerance is regenerated or weakened.

- New independent-oracle suites:154 passed. Earlier exactness suites remain in the
  required current full CI matrix.
- Real offline native integration:53 passed, no skips (589.57s), covering reference,
  count, artefact filtering, strand/read-number replay, competitive evidence, BLAST,
  deduplication, QC, visualization, doctor and validate-run. R4_51 is selected explicitly;
  its actual emptyDrops A/A computation is the separate B12 receipt.
- Real dustmasker preparation: baseline A/A and combined A/B prepared FASTA bytes
  match; inputs including GTF unchanged. Manifest raw hashes retained; no whole
  manifest/directory identity claim for timestamp/source metadata.
- Retained development slice:1024 rows, six fixed-layer method labels,12/12 raw
  reducer files match independent pre-performance-main `6910c2d` A/A and combined
  A/B. These are not six allocator refits. The oldest historical12 files have a
  different schema (`locus_id`, `n_unique_positive_cells_called` were added on main)
  and zero raw matches; they remain unchanged. New baseline evidence is a separate
  receipt, never a replacement golden.
- Final combined cell/group/strata refits:27/30 raw files identical, three source-
  provenance manifests differ. The comparison is recorded in
  `c00/combined_downstream.json`, including every raw hash and precise literal
  executing-source spans. Baseline A/A must be raw identical; only the predeclared
  `software.build_sha256`, `software.git_commit`, `software.source_dirty` values
  may differ. This is not a whole-directory identity claim.
- Ruff check/format275files, mypy60sourcefiles and data-governance checks pass.
  Final PR full CI (including cross-platform unit/native/DAG gates) must pass before
  merge; final revision and check results are recorded on that PR and in local receipts.

Failed runner attempts are retained: an initial path-object hashing mistake, followed
by the retained historical-schema assertion above. Successful matching-baseline
receipts supersede neither the old logs nor the historical golden.

## Repeated performance and remaining gates

G01 historical mutation receipt passes:148 saved comparison artifacts unchanged,
nine control/mutation checks. Six legacy timing groups now each have seven alternating
fresh-process pairs in two sessions; all17 component/fixture rows, spreads, source/
input hashes and every trial are retained in `g02/REPEATED_TIMINGS.md` and JSON.
Representative wall ratios in sessions1/2: count formatting0.2034/0.2014,
evidence summary0.6691/0.6621, enrichment0.0497/0.0507, grouped detection0.6702/0.6652,
reference preparation0.6259/0.6239, concatenation0.7170/0.7552. Memory-focused controls
remain within the median5% regression gate; no new10% runtime claim for B06/B23/B18.
The B24 timing uses a fixture dustmasker emitter; real dispatch is qualified separately.
These are shared-host synthetic component/stage results, not a quiet-host claim.

The implementation assessment has13 accepted cards,15 deferred cards and two
separately gated resource cards. Failed early gates stop broader grids; unchecked
rejected edit steps stay unchecked. G02 production stage attribution/whole-run
measurement and B11/B30 resource-policy qualification remain open. No component
ratios are multiplied. The bounded whole-run and resource-policy proposals are
`2026-10-19-performance-whole-run-proposal.md` and
`2026-10-19-performance-resource-policy.md`; no jobs or changed workers are authorized.

Detailed local receipts, raw artifacts, scripts, source/input hashes and failed
attempt logs are retained under the root checkout's
`analysis/computational_bottlenecks_2026-10-17/followups/20261018-remaining-candidates/`.
