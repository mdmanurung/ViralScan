# Explicit legacy workers: bounded qualification

The user approved this local phase on 2026-10-19. Baseline is merged main
`9f0c397d484759a0ddf83430fb1455dcfabb0bb5`; the isolated candidate adds only
`kb count -t 8` in both whitelist branches, `kb ref -t 8` in both reference
entrypoints, and homology `minimap2 -t 3`. These preserve the installed native
worker defaults. Ruff also reformatted two existing expressions without changing
their AST. No dependencies, public options or scheduling declarations were added.

## Native evidence

Immutable local fixture inputs were hashed before/after execution. Each of six
stages ran original A/A and candidate A/B with identical paths, numerical-library
environment, tool binaries, CPU affinity and scientific options. The reference
fixture includes both reference entrypoints; count covers packaged and provided
whitelists. Separate deterministic two-record reference inputs exercise homology
and CD-HIT clustering.

| Stage | Files per arm | Raw equal A/A | Raw equal A/B |
|---|---:|---:|---:|
| Count, packaged whitelist | 16 | 13 | 13 |
| Count, provided gzip whitelist | 15 | 12 | 12 |
| Menu reference build | 7 | 5 | 5 |
| Build-reference entrypoint | 8 | 6 | 6 |
| Homology TSV | 1 | 1 | 1 |
| Cluster FASTA and cluster membership | 2 | 2 | 2 |

All scientific/native artifacts match raw bytes: BUS, EC, matrix, H5AD, index,
t2g, prepared FASTA, reproducibility report, homology TSV and cluster outputs.
Ten operational files per comparison differ. They are retained unmodified,
including their raw SHA-256 values. No complete-directory byte identity is claimed.

The comparator permits only individually validated JSON value spans for timestamps,
resource measurements, random staging paths and the explicit worker argument,
plus the two executing-source identity fields. All JSON literal bytes outside those
spans remain exact. Logs permit only timestamp prefixes and exact random staging
path spans. Every other field, type, file and literal byte is checked. Five controls
reject a changed worker, scientific byte, missing file, metadata whitespace or an
unlisted numeric field. Original fixtures, goldens and comparators are unchanged.

A native minimap2 missing-input check preserves exit code, empty stdout, absent
partial output, exact error messages and all non-timing stderr bytes. Full stderr
is already variable in original A/A because minimap2 reports `seconds*CPU` timing
prefixes; these raw streams remain available. The candidate argv differs only by
`-t 3`. Mocked boundary tests additionally verify full exception/stderr propagation.

## Scheduler and failure contracts

The 16 actual-rule cases span 1/6/8/32 cores, both whitelist paths and external
`--set-threads kb_count=2`. An argv/environment recorder checks the explicit native
request remains eight and inherited OMP/OpenBLAS/MKL/NumExpr requests stay one,
or `min(cores, 2)` under the external override. These are scheduler/argv contract
tests using a fixture producer; they are separate from native byte evidence.

Existing DAG tests cover optional host filtering, read filtering, anello alignment,
auto evidence, competitor evidence and host capture. Existing rerun tests exercise
error, missing/empty result and cancelled-promotion failures for both whitelist
paths, preserving previous counts, user metadata and temporary-file cleanup. CLI
and rerun-multimap regressions preserve the shared Snakemake launch contract.
A separate paired dry-run matrix compares default and all optional rules together
at every core budget, both whitelist paths and the external override. Arbitrary
subsets of optional rules are not exhaustively crossed.

The targeted default regression suites passed 237 tests (35 integration cases
deselected), and rerun-multimap/tool-resolution/reference-pipeline regressions
passed 48 tests. All 34 scheduler/failure cases passed, and all 32 paired DAG
cases (64 dry runs) preserved reservations. Initial native CI passed 87 tests
with 23 skips. Its first full-suite run caught stale registered source hashes;
the source and claim-registry evidence pins were refreshed against immutable
commits, preserving the `implemented_not_calibrated` claim status. A non-release
provenance tag retains these commits after rebase. Runtime source is unchanged.
Ruff check, Ruff format (276 files), mypy (60 source files) and
data governance passed. The final local native count/strand/lineage/end-to-end
suite passed 54 tests with no skips; governance/worker regressions passed 27
(16 integration cases deselected in that run). Final PR CI is the merge gate. The first count sampler missed descendants
launched by Snakemake's worker thread. A supplemental audit traverses every task's
child list; the original incomplete observations are retained, not presented as
measured worker counts. The supplemental audit observed native count bus/sort
requests of eight in every arm, with sampled bus peaks of one or two threads.
Reference index peaks were two to four, homology one to four and CD-HIT eight.
These transient peaks are not equal to requested workers. Tiny fixtures and 10 ms
sampling cannot establish sustained
production CPU utilization or catch every brief worker peak.

## Remaining boundary

This qualifies explicit legacy requests, without claiming a speedup or total-budget
compliance. A fixed eight-worker count request still exceeds one/six admitted cores;
CD-HIT still requests all visible CPUs with `-T 0 -M 8000`. A total-admission policy
must choose between rejection and changed effective workers, both compatibility
changes requiring a later decision. B11/B30 therefore remain partial. Production
G02 attribution and whole-run timing remain unmeasured; no SLURM jobs, downloads,
held-out/protocol-frozen inputs or retained-data deletion were used.

Raw receipts, manifests, failed harness attempts and the comparator are retained at
`analysis/computational_bottlenecks_2026-10-17/followups/20261019-084842-resource-legacy/`
in the shared research checkout. They are not shipped as package assets.
