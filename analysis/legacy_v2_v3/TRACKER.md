# ViralScan legacy v2 versus v3 diagnostic tracker

Status: **active**

Frozen: 2026-07-24

Branch: `codex/viralscan-v3`

Baseline Git commit: `a8f32876e6ec4a5d8cf29b21aaf325ccaa50e64f`

This tracker covers a diagnostic comparison of archived ViralScan 2.2.0
outputs with v3 molecule-aware counting. Historical outcomes have already been
viewed. They are excluded from confirmatory v3 validation and cannot close
truth-panel, comparator, calibration, release, publication, or
package-superiority gates.

## Status and evidence rules

- `[ ]` pending;
- `[~]` active;
- `[x]` complete with evidence;
- `[!]` blocked or failed, with the blocker or failure retained.

A task may be marked `[x]` only when its command, input hashes, output path, and
validation evidence are recorded. Unfavourable, unsupported, and failed rows
remain in the denominator. Technical repeats are evaluated for reproducibility
but are not independent biological samples.

## Next action

Attempt-2 arrays `25331035` and `25331037` are terminal. All ten rows failed and
every outcome is retained; see the attempt-2 terminal outcome section below.
Repair the three causes, re-validate the four v2 rows that actually succeeded
without rewriting attempt-2 evidence, and freeze a v3-only attempt-3 packet
**without submitting it**. This diagnostic is outcome-ineligible and yields
priority to `SCI-03`; the matched-FASTQ comparison at `LVC-13`–`LVC-14` stays
unstarted and still has no comparison tooling.

## Task status

| ID | Status | Deliverable or exit gate |
|---|:---:|---|
| LVC-00 | `[x]` | Freeze execution packet, tracker, and diagnostic-only protocol. |
| LVC-01 | `[x]` | Inventory and hash the exact 44-row/42-input frozen cohort. |
| LVC-02 | `[x]` | Freeze v2, v3, tool, source, and reference identities. |
| LVC-03 | `[x]` | Implement manifest-driven run, summarize, and validate behavior. |
| LVC-04 | `[x]` | Pass focused, existing, lint, and full-suite verification. |
| LVC-05 | `[x]` | Pass deterministic low-count skin pilot without source mutation. |
| LVC-06 | `[x]` | Pass the SRR12682296 scale pilot and resource gate. |
| LVC-07 | `[x]` | Execute all 44 BUS rows through bounded SLURM arrays. |
| LVC-08 | `[x]` | Validate every successful row and retain every failure. |
| LVC-09 | `[x]` | Generate sanitized BUS comparison tables. |
| LVC-10 | `[x]` | Freeze identical five-control FASTQs after network approval. |
| LVC-11 | `[!]` | Rerun and audit all five controls with v2.2.0. Blocked: `SRR6825024` out-of-memory at the frozen tier ceiling. |
| LVC-12 | `[!]` | Rerun and validate all five controls with frozen v3. Blocked: the pinned viral-panel Zenodo DOI is unregistered, so the cache cannot be fetched or pinned. See `REF-11`. |
| LVC-13 | `[~]` | Interpret validated results within the frozen claim boundary. |
| LVC-14 | `[~]` | Re-audit, freeze hashes, and leave a restart handoff. |

## LVC-00 evidence

Status: `[x]` completed 2026-07-24.

Commands:

```text
git status --short --branch
git rev-parse HEAD
sha256sum PLAN.md docs/review-clear-execute-plan.md docs/review-clear-execute-tasks.md docs/handoffs/2026-07-24-legacy-v2-v3-execution.md docs/handoffs/handoff-CFtSVg.md
python3 -c "import pathlib, yaml; ..."
python3 -c "import pathlib; forbidden=('/'+'exports/','/'+'home/','hg-'+'funcgenom','para-'+'lipg'); ..."
git diff --check -- analysis/legacy_v2_v3
```

Inputs:

| Relative path | SHA-256 |
|---|---|
| `PLAN.md` | `46c0064139a6f7675d063043f4be66f3cf5c9117dac61f5d0f9d4288d9873818` |
| `docs/review-clear-execute-plan.md` | `27e8832e1eb7ff118d771583cb260e6a5a076e7fc386fde4d8333876ae5064dd` |
| `docs/review-clear-execute-tasks.md` | `166204b9956f4f5b505496adebd1fa769302cdf2f38b7f2d96398ff641dd72e7` |
| `docs/handoffs/2026-07-24-legacy-v2-v3-execution.md` | `ed6215a0148876cf7c1d702d9075efc8a306100730cdbfa2b681aefe2d50ec17` |
| `docs/handoffs/handoff-CFtSVg.md` | `d0c5df6267f1a7fa656a7a63a1528e6736211c7fbfe1990a3f4188aa2da1a15b` |

Committed 2026-10-02:
- `review-clear-execute-plan.md` still matches its frozen hash.
- `review-clear-execute-tasks.md` is now `5e9b80328fa0cecf…`: its checkboxes were ticked after the freeze. No frozen copy survives, so the exact diff can't be reconstructed. The ticks lag this tracker: LVC-03's "no institutional absolute paths" item is unchecked, and the LVC-11/12 items predate attempt 2.

Outputs:

- `analysis/legacy_v2_v3/TRACKER.md`;
- `analysis/legacy_v2_v3/protocol.yaml`;
- `analysis/legacy_v2_v3/README.md`.

Preserved pre-existing worktree changes:

- modified `.living/INDEX.md`;
- modified `.living/log/2026-07-21-001-viralscan.md`;
- modified `.living/log/LOG_REGISTRY.md`;
- untracked `.living/log/2026-07-21-002-viralscan.md`;
- modified execution-packet plan and task files;
- untracked `docs/handoffs/` artifacts.

No archived result, installed v2 environment, script, test, or unrelated
worktree file was modified by LVC-00.

Validation evidence:

- protocol YAML parsed successfully;
- frozen cohort totals reconciled to 44 technical rows and 42 logical inputs;
- the diagnostic claim boundary was present and confirmatory eligibility was
  false;
- prohibited EM and unique-weighted outcome arms were present;
- the tracked path-sanitization scan returned no institutional absolute paths;
- `git diff --check -- analysis/legacy_v2_v3` passed.

## LVC-01 evidence

Status: `[x]` completed 2026-07-25.

- Raw manifest: ignored `benchmark_runs/legacy_v2_v3/inventory/raw_manifest.tsv`,
  SHA-256 `1669409ad889d33a66fb3efbfa2a1be2cf6269aa5120c0a5ded23e1fd7aa6ae8`.
- Sanitized manifest: `analysis/legacy_v2_v3/cohort_manifest.tsv`, SHA-256
  `37c9b6df89e0d8260446cd5ea15048e9b0faefaa6023547685a88b98cd3c408c`.
- Verified 44 technical rows, 42 logical inputs, six 10x v2 and 38 10x v3
  whitelists.
- Verified two-member repeat groups for SRR6825025 and
  KCL10525740 S3 L002-input.
- Every required config, summary, BUS, EC, transcript list, whitelist,
  multimap H5AD, base H5AD, barcode list, gene IDs, gene names, run metadata,
  viral-feature list, index, t2g, and transcriptome exists and has a size and
  SHA-256 record.
- The tracked manifest contains no institutional absolute path.

## Active pilot-safety gate

No pilot job has been submitted. Independent review found that the first
benchmark draft allowed `prepare_resolved_bus` to default its corrected BUS
beside the archived raw BUS. A second review found that `run-row` created its
output before containment validation and did not hash every file actually
consumed. Both defects now have regression-tested fixes: corrected BUS writes
go to scratch, containment is checked before output creation, and the refreshed
manifest covers every local file consumed. The archived skin-pilot BUS and H5AD
still match the frozen manifest, and no generated corrected/resolved BUS exists
in its source directory.

The frozen full-profile environment and runtime source/environment hash checks
subsequently passed. The earlier BUS-only snapshot was not used for the pilot.

## LVC-02 evidence

Status: `[x]` completed 2026-07-25.

- Emma's read-only environment: ViralScan 2.2.0, kb-python 0.29.5,
  kallisto 0.51.1, bustools 0.45.1.
- Frozen v3 wheel SHA-256:
  `72ce891e448ac27873922034804aa57dfd8cf8ec4d5d738103f17bef00de4b03`.
- Frozen comparison script SHA-256:
  `374757159f536b1cb8758d9599c14b9e1c113c725816ae172aea5832289eae17`.
- Frozen BUS benchmark script SHA-256:
  `fbe0c64ab7b9e4c7db207046cc888afc78eea8fdde810c137e5f38d1db0f435b`.
- Full environment differs from the tracked `environment.yml` only because
  Bioconda no longer exposes DropletUtils 1.24.0; the available R 4.4 build,
  DropletUtils 1.26.0, was frozen before pilot submission.
- Resolved environment spec SHA-256:
  `d7a36278fb8f163c28c2b5619d4f76da40e5f06d3820f231b2866f858aa521ec`.
- Explicit Conda export SHA-256:
  `53871cc7e85f8e9966b8746f8bd70a90e63a5751319c44fc23fcf6b635bbc39f`.
- Full-profile doctor passed; archived JSON SHA-256:
  `c8132219806b94c7d701504c8ceb3addfdc105cdaacacfdac526ac3e2334cf83`.
- Execution provenance SHA-256:
  `380f7eba49c574abeeb6914955c01014c570d4793ca2ecaae835c548152f9bac`.

## Verification and pilot state

- Focused comparison tests: 15 passed.
- Focused benchmark tests: 10 passed.
- Existing BUS preparation tests: 2 passed, 12 deselected.
- Repository suite: 732 passed, 36 deselected.
- Ruff, diff check, data governance, and v3 protocol validation passed.
- Skin pilot attempt 1, SLURM job `25329961`, failed safely with exit code 1
  after 48 seconds and 585,192 KiB peak RSS. Legacy reconstruction completed
  and retained the expected EBNA-2 total of 2.0.
- The failure occurred before BUS preparation or allocation. The frozen
  benchmark script assumed package modules were present below
  `source_dist/src/`, but the execution snapshot contains the installed v3
  wheel. The missing fingerprint target was
  `src/viralscan/scripts/multimap.py`.
- All 16 recorded source, reference, BUS, H5AD, and configuration hashes
  matched before and after the failed attempt. The failure/status records,
  scheduler logs, and resource record remain under
  `benchmark_runs/legacy_v2_v3/pilots/skin_pilot1/`.
- Pilot 2, the scale pilot, and all full arrays remain unsubmitted.
- The installed-wheel fingerprint fix passed its focused regression and
  real frozen-environment import checks. Generated Python bytecode is excluded
  from the execution identity and disabled in the retry job. Refreshed
  benchmark SHA-256:
  `b115f5236b8f61a7a5322560783d72f03eeb09fc092009cd6791050d9d33cc53`;
  refreshed execution-provenance SHA-256:
  `256899e5ffed27d0ad8b6dfdfc8b7a2097e76af5f12282ea3dc3dd31dda9491e`.
- Skin pilot attempt 2 completed successfully as SLURM job `25329976`
  (exit 0, elapsed 4:08, batch MaxRSS 2,095,288 KiB). It processed 617,076
  cells in one allocator pass; the v3 scientific-parameter hash is
  `0b6ac0da303ad2d13a6ac4c2e8662fd1a7bfff9f77587b5895bdfed6beb33ef2`.
  Input fingerprints before and after, plus the outer archived-source hash
  checks, are identical. Scratch was cleaned and not retained.
- The independent deterministic repeat was submitted to a new output root as
  SLURM job `25329977` and completed successfully (exit 0, elapsed 4:13,
  batch MaxRSS 2,028,152 KiB).
- The deterministic comparison report at
  `benchmark_runs/legacy_v2_v3/pilots/skin_pilot_determinism.json` is a match
  across 11 gates: scientific parameters, input hashes, installed code,
  virus grouping, count audit, per-virus and per-cell tables, H5AD obs/var,
  X, and all layers. Focused comparator tests passed 3/3.
- LVC-05 is closed. Both valid attempts used the archived 10x v3 whitelist,
  completed in one allocator pass, satisfied matrix and molecule invariants,
  cleaned scratch, and left every recorded archived-source hash unchanged.
- LVC-06 preflight verified the SRR12682296 archived 10x v2 whitelist at
  SHA-256
  `b0dda4b114d8fc7ea8def94fb472e6205c1831a72e3ac245ed13cc332d9f20d1`
  and a 9,390,038,355-byte scratch requirement. Scale-pilot job `25329992`
  was submitted with 8 CPUs, 32 GiB, and a 6-hour limit.
- The scale analysis completed successfully in one allocator pass: 53,709,866
  input molecules, 50,952,132 unique, 1,388,191 ambiguous, 1,369,543
  unresolved, and 52,340,323 selected matrix mass. All 99 legacy feature
  totals matched, inner and outer source/input hashes were unchanged, and
  scratch was cleaned.
- SLURM job `25329992` ended `FAILED` only because the final shell-only status
  print contained invalid f-string escaping after all scientific and hash
  gates had passed. The corrected post-check returns zero against the retained
  `success` status. The engineering failure is retained and the expensive
  scientific computation was not rerun.
- GNU time recorded 6:25.66 wall time, 4,231,476 KiB peak RSS, and
  5,370,306,338 scratch bytes before cleanup. The below-4-GiB resource tier is
  retained. The descriptive baseline comparison and its no-whitelist,
  non-exchangeability limitation are recorded in `EBV_SCALE_PILOT.md`.
- The existing streamed allocator row-order/buffer-size invariance test passed.
  LVC-06 is closed. Full arrays remain unsubmitted pending frozen-packet
  hardening.
- Frozen-packet hardening passed 7/7 focused tests. The full default repository
  suite passed 738 tests with 46 deselected; Ruff, diff check, data governance,
  and v3 protocol validation also passed.
- The full-array execution freeze now includes the hardened helper. Execution
  provenance SHA-256:
  `3100a53ef69c627818224fb223d1eeab80476509ce559b9f2ed1a87ed79ad67e`.
  The frozen packet has 44 rows, 41 below-4-GiB and three at-least-4-GiB,
  maximum concurrency four, and no pre-existing row outputs.
- The first 41-row below-4-GiB array was submitted as SLURM job `25330030`,
  then canceled after a pre-validation review found that the packet emitted
  direct v3 directories rather than the full wrapper/legacy/v3 layout required
  by the frozen `summarize` and `validate` commands. Tasks 0, 1, 2, and 6 were
  canceled while active; task 5 completed; the remaining tasks were canceled
  before execution. All partial outputs and logs are retained under
  `slurm_packet_frozen/`.
- The mismatch was caught before broad execution. A replacement packet will
  invoke the frozen full `run-row` comparison harness in a new run directory.
  The corrected helper passed 9/9 focused tests and was frozen at execution
  provenance SHA-256
  `51cb4d16e5a2e1e80d9dd5c41523064fcb57459c363adf79e135f212122ba2d1`.
- Replacement below-4-GiB array job `25330045` was submitted at maximum
  concurrency four. Its first tasks produced the required wrapper layout
  (`hash_check_before.json`, `legacy/`, and `v3/`) and entered BUS allocation.
- Job `25330045` completed 41/41 rows successfully. Interim complete-root
  validation reported exactly 41 successes and three missing planned rows,
  with no error in a completed row. The missing rows were precisely the three
  prespecified at-least-4-GiB tasks.
- The three-row 128-GiB large-tier array was then submitted as SLURM job
  `25330512`. All three tasks completed successfully with exit code zero
  (elapsed 14:11, 4:47, and 4:36).
- The replacement full-wrapper run root contains 44/44 `success` status
  records: 41 from job `25330045` and three from job `25330512`. No scientific
  parameter was changed between tiers, maximum concurrency remained four, and
  the canceled direct-v3 attempt remains retained separately.
- LVC-07 is closed. Aggregate generation retained all 44 planned rows, emitted
  4,664 run-virus rows over the frozen 106-virus universe, and produced a
  header-only failure table. Full-root validation is still the LVC-08 gate.
- The frozen full-root validator then passed with `status: valid`, 44 planned
  rows, 44 `success` rows, and an empty error list. This closes LVC-08.
- Sanitized copies of the validated aggregate and validation evidence were
  promoted under `analysis/legacy_v2_v3/`:

  | Artifact | SHA-256 |
  |---|---|
  | `legacy_reproduction.tsv` | `edad8e74af9409da196d1856fd15694d64345004b0589170e51e526c937bbf88` |
  | `run_metrics.tsv` | `c7e7c5c2177b1083d1449ff3db28128c25369660c9787e80e6907b62ed4655a3` |
  | `virus_metrics.tsv` | `b7b5aa0f134ab2b57435a60dffb3a3fd9237b93e57287eddcc444b4cd40c6135` |
  | `failures.tsv` | `a25406a8c914d43603b541ae2303a8524e613c537a56796e8a56a269d47130d8` |
  | `validation.json` | `e6ca98fa1ae960713a687d9a68e24686c1a7cd66771ab4668c738957b74956e4` |

- The cell-concordance summarizer passed its regression test and generated 44
  rows. Every legacy/v3 barcode universe matches exactly. A separate
  technical-repeat audit records zero L1 distance for unique, equal, and
  host-conservative virus vectors in both declared pairs.
- `cell_concordance.tsv` SHA-256:
  `b6f17d8063fb36b394bbea8bf12f81c10566b5d36f4cdb0f0ae7ffb9e092cd47`.
  `technical_repeat_audit.tsv` SHA-256:
  `18af6da692ab3dbe6a395e04777d94a691e24e2fceac6b8526e37e7d154fe34f`.
  This closes LVC-09.
- `REPORT.md` was drafted only from validated and audited tables. The report
  does not promote these diagnostic findings to truth, performance, release,
  or publication claims.
- `candidate_persistence.tsv` classifies all 101 union run-virus entries in
  skin as persisted, lost, newly nonzero, or absent separately for every v3
  endpoint. Its SHA-256 is
  `12a5421d7732b1c2c5583a0184cd1b059a375ab126b36b2d7f480f39b3d91d4e`.
  Seven same-sample host-conservative lane recurrences are present, all for
  Cercopithecine herpesvirus. LVC-13 remains active only because the fresh
  matched-FASTQ comparison has not yet run.
- The three archived EBV FASTQ pairs referenced by the historical configs are
  absent. The two HIV-control pairs are present but total approximately
  307 GB. The streaming pair-integrity, MD5, and SHA-256 auditor passed 2/2
  focused tests and was submitted as bounded two-task SLURM array `25330549`.
  LVC-10 remains active; no network download was attempted.
- Array `25330549` completed both tasks successfully. SRR6825024 contains
  498,768,252 valid paired records and SRR6825025 contains 284,283,729; both
  pairs passed structure, paired-identifier, byte-count, MD5, and SHA-256
  auditing.
- Authoritative ENA metadata, paired layout, URLs, expected sizes, and MD5
  values were retrieved for SRR12682296, SRR12682297, and SRR12682298. Their
  six compressed FASTQs total 37,451,192,938 bytes.
- Download array `25330802` failed before transferring data because the
  compute-node curl lacks `--retry-all-errors`. The failure and empty logs are
  retained. The compatibility option was removed without changing inputs or
  checksum gates, and corrected bounded array `25330814` was submitted with
  two concurrent tasks; tasks 0 and 1 entered `RUNNING`.
- Corrected array `25330814` completed all six EBV mates with exit code zero.
  Every stored byte count and ENA MD5 matched the frozen metadata; ignored
  status rows additionally record SHA-256.
- A gzip-aware pair auditor now validates decompressed FASTQ structure,
  paired identifiers, and the declared 10x-v2 16-base-CB plus 10-base-UMI
  minimum without materializing uncompressed scratch. Three-sample EBV array
  `25330878` and two-sample retained-HIV chemistry re-audit `25330881` were
  submitted with bounded concurrency. LVC-10 remains active until all five
  terminal audit records agree and the shared manifest is frozen.
- Final verification after adding the tracked cell and candidate summaries and
  FASTQ auditor: 742 default tests passed, 48 were deselected, and the 59
  warnings were the existing host-response convergence warnings. The focused
  diagnostic set passed 22 tests with 27 deselected. Ruff, data governance,
  tracked path sanitization, and `git diff --check` passed.
- Verification after the gzip/chemistry auditor and manifest freezer:
  745 default tests passed, 48 were deselected, and the same 59 pre-existing
  host-response convergence warnings remained. Data governance, v3 protocol
  validation, tracked path sanitization, Ruff, and `git diff --check` passed.
- Arrays `25330878` and `25330881` completed all five full-stream pair audits
  with exit code zero. Every pair has matching identifiers, exact 10x-v2
  geometry (26-base R1; 98-base R2), and frozen storage/content digests.
  `analysis/legacy_v2_v3/control_inputs.tsv` contains the five sanitized rows;
  its SHA-256 is
  `1424504d3670a52243adce473345f09e66ed95ea8012278d02d4c278995b8433`.
  This closes LVC-10.
- The fresh runner refuses prior output/status/log replacement, carries the
  audited input sizes and hashes into every task, and verifies the immutable
  execution packet before running. The post-freeze full suite passed 758 tests
  with 48 deselected and the same 59 existing convergence warnings.
- Fresh attempt 1 arrays `25331024` and `25331026` failed before ViralScan
  execution because Slurm copied `run_task.sh` to its spool and the script
  derived the wrong packet root. All ten scheduler failures and logs are
  retained; no workflow status or output was created.
- Attempt 2 passes an explicit immutable packet root, with a regression test
  that executes the runner from a separate simulated Slurm spool. Small array
  `25331035` (six EBV tasks, maximum four concurrent) and dependency-gated
  large array `25331037` (four HIV tasks) were submitted.

## Attempt-2 terminal outcome (2026-07-26)

Both arrays are terminal. All ten planned rows are retained; none was excluded.
Evidence: `benchmark_runs/legacy_v2_v3/fresh_control_packet_attempt2/status/*.json`,
`sacct -j 25331035,25331037`, the packet log directory, and the on-disk output
trees under `benchmark_runs/legacy_v2_v3/fresh_control_runs/`.

| Row | Array task | Recorded | Terminal disposition |
| --- | --- | --- | --- |
| `v2__SRR12682296` | `25331035_0` | `failed`, exit 65 | Workflow succeeded (`workflow_exit_code` 0); validator defect |
| `v2__SRR12682297` | `25331035_2` | `failed`, exit 65 | Workflow succeeded (`workflow_exit_code` 0); validator defect |
| `v2__SRR12682298` | `25331035_4` | `failed`, exit 65 | Workflow succeeded (`workflow_exit_code` 0); validator defect |
| `v2__SRR6825025` | `25331037_2` | `failed`, exit 65 | Workflow succeeded (`workflow_exit_code` 0); validator defect |
| `v2__SRR6825024` | `25331037_0` | `failed`, exit 1 | Genuine out-of-memory, retained failure |
| `v3__SRR12682296` | `25331035_1` | `failed`, exit 1 | Missing viral-annotation cache |
| `v3__SRR12682297` | `25331035_3` | `failed`, exit 1 | Missing viral-annotation cache |
| `v3__SRR12682298` | `25331035_5` | `failed`, exit 1 | Missing viral-annotation cache |
| `v3__SRR6825024` | `25331037_1` | `failed`, exit 1 | Missing viral-annotation cache |
| `v3__SRR6825025` | `25331037_3` | `failed`, exit 1 | Missing viral-annotation cache |

Three independent causes:

1. **Nested-output validator defect.** `_v2_artifact_errors` in
   `scripts/run_fresh_control.py` checked `config.yaml`, `summary.txt`, and the
   `kb-python` artifacts directly beneath the task output directory, but the
   legacy 2.2.0 CLI writes them one level deeper under the sample identifier.
   The check found nothing, reported six phantom missing artifacts, and
   manufactured exit 65 on four rows whose workflow had exited zero. All six
   required artifacts exist and are non-empty for each of those four rows.
2. **Unpopulated viral-annotation cache.** Every v3 row raised
   `ViralScanDataError` during config creation. The frozen packet never ran
   `viralscan data fetch`, so the Zenodo-backed annotation cache
   (DOI `10.5281/zenodo.20112332`) did not exist in the execution environment.
   No v3 row reached quantification.
3. **Genuine out-of-memory.** `v2__SRR6825024` was killed in the legacy
   `multimap.py` at roughly 121.4 GiB peak resident set against the frozen
   128 GiB tier ceiling, after 5 h 52 m. This is a real ViralScan 2.2.0
   resource limitation on the large HIV control and is retained as a failed row
   under `sensitivity-failure-description`.

### Attempt-3 cache pin is blocked

The approved remedy for cause 2 was to pre-fetch the viral annotation panel on a
network-capable node, pin its `manifest.json` SHA-256 into the packet, and
export `VIRALSCAN_CACHE` for every task. The fetch cannot be performed: the
Zenodo record pinned in `src/viralscan/data_fetch.py`
(`VIRAL_DATA_DOI = "10.5281/zenodo.20112332"`) is not registered.

Verified from `res-hpc-exe029`, which has working outbound HTTPS:

- `https://zenodo.org/api/records/20112332` returns
  `{"status": 404, "message": "The persistent identifier is not registered."}`
- `https://doi.org/10.5281/zenodo.20112332` returns 404
- an unrelated third-party Zenodo DOI referenced in
  `scripts/ttv_public_datasets.json` returns 200 from the same host, so this is
  not a network or egress restriction

The 195 panel GTFs are still present in `src/viralscan/data/` in the source
tree, but the installed package in `env_full` carries only
`anellovirus_accessions.tsv` and `__init__.py`. Tracked as `REF-11` in
`PLAN.md`. Until it is resolved, no cache pin can claim the DOI provenance that
`cache_valid` requires, so the attempt-3 packet is not frozen. All runner,
manifest, and wrapper wiring for the cache is complete and tested.

Two contract notes recorded here because they were not documented previously:

- The status payload written by `run_control` satisfies none of `stage`,
  `attempt_id`, or `scientific_parameter_hash` from `protocol.yaml`
  `required_failure_fields`, so all ten attempt-2 records are non-compliant as
  failure records. `scientific_parameter_hash` is satisfied by the existing
  `command_sha256`, matching the convention used by the identical-BUS arm.
- Tier assignment for fresh runs uses `LARGE_INPUT_BYTES` in
  `scripts/prepare_fresh_controls.py` — a raw-FASTQ byte proxy — because no BUS
  file exists before a fresh run against which to evaluate the protocol's
  literal `raw-bus-size-at-least-4-gib` condition. This proxy placed
  `SRR6825024` in the 128 GiB tier and is unchanged for attempt 3.
