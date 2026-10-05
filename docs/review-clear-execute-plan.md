# ViralScan v2.2.0 versus v3 legacy diagnostic execution plan

## Objective

Build and execute a reproducible diagnostic comparison between Emma Vonk's
archived ViralScan 2.2.0 outputs and the current `codex/viralscan-v3` molecule
counting implementation.

The comparison has two arms:

1. Reprocess all 44 archived technical result trees from their retained BUS,
   EC, transcript, whitelist, t2g, matrix, and summary artifacts.
2. Rerun the three EBV and two HIV public controls end to end from identical
   FASTQs and the original combined reference.

The 44 technical trees represent 42 logical inputs after grouping the second
SRR6825025 run and the repeated KCL10525740 S3 L002 input as technical repeats.

## Scientific status

This is a diagnostic/development audit. Historical outcomes have already been
viewed, so these rows are ineligible as confirmatory v3 holdout evidence and
cannot close `G5`, release, publication, truth-panel, comparator, calibration,
or package-superiority gates.

All rows and failures remain in the result set. No threshold, allocation method,
feature set, barcode set, or inclusion rule may be changed after outcome
inspection.

## Repository and preservation rules

- Work in `/exports/archive/hg-funcgenom-research/mdmanurung/ViralScan`.
- Preserve the existing `.living/**` changes and `docs/handoffs/` artifacts.
- Treat `/exports/archive/hg-funcgenom-research/evonk/viralscan/results` and
  Emma's installed v2.2.0 environment as read-only.
- Store large artifacts, environments, BUS intermediates, H5AD files, FASTQs,
  scheduler logs, and absolute-path manifests under the ignored
  `benchmark_runs/` tree.
- Store only relative identifiers, hashes, small tables, the tracker, protocol,
  and diagnostic report under `analysis/legacy_v2_v3/`.
- Do not write institutional absolute paths into tracked or public artifacts.
- Do not use `viralscan rerun-multimap` for the legacy trees because it copies
  complete source result directories. Use a read-only manifest-driven harness.

## Frozen cohort

The source inventory must resolve exactly:

- 44 technical result directories;
- 42 logical inputs;
- 3 EBV control runs: SRR12682296, SRR12682297, SRR12682298;
- 2 logical HIV controls: SRR6825024, SRR6825025;
- 38 skin technical outputs representing 37 unique skin input pairs;
- 44 archived chemistry whitelists: six 10x v2 and 38 10x v3.

Technical repeat groups:

- the active and `old/` SRR6825025 outputs;
- `WS_SKN__KCL10525740/S3/L001/WS` and
  `WS_SKN_KCL10525740/S3/L002/WS`, which point to the same L002 FASTQ names but
  have different retained BUS hashes.

Technical repeats are evaluated for reproducibility but never counted as
independent biological samples.

## Frozen comparison layers

For every archived BUS tree, report:

1. `legacy-v2-reconstructed`: the v2.2.0 detection matrix
   `counts_corrected + counts_original` and its strict `>1` gene threshold.
2. `v3-unique`: integer gene-unique CB-UMI molecules.
3. `v3-equal`: v3 molecule-aware equal allocation. This is the closest
   diagnostic analogue to legacy equal splitting and isolates the counting
   contract change.
4. `v3-host-conservative`: the current product-default allocation for a
   combined host-plus-virus reference.

Do not add EM or unique-weighted outcome arms to this audit. The implementation
may retain their automatically computed layers as diagnostics, but they are not
reported comparison endpoints.

## Implementation

### 1. Tracking and protocol

Create `analysis/legacy_v2_v3/` containing:

- `TRACKER.md`;
- `protocol.yaml`;
- `cohort_manifest.tsv` with relative identifiers only;
- final small result tables and `REPORT.md`.

Tracker states are `[ ]` pending, `[~]` active, `[x]` complete with evidence,
and `[!]` blocked or failed. A task is complete only after its command, input
hashes, output path, and validation evidence are recorded.

### 2. Inventory command

Add a private benchmark CLI in `scripts/compare_legacy_v2_v3.py` with
subcommands:

- `inventory`;
- `run-row`;
- `summarize`;
- `validate`.

`inventory` scans the read-only source root and emits a raw absolute-path
manifest under `benchmark_runs/` plus a sanitized tracked cohort manifest.
Required fields include run and logical IDs, technical-repeat group, chemistry,
sample class, expected target when externally documented, relative paths,
artifact sizes and SHA-256 values, reference hash IDs, and tool versions.

Inventory must fail closed if the expected row counts or required artifacts do
not match the frozen cohort.

### 3. Legacy reconstruction

Reconstruct v2.2.0 calls from each archived `adata_multimap.h5ad` using the
installed 2.2.0 semantics. Compare every reported gene and virus value with
`summary.txt` using absolute tolerance `1e-6`.

Missing or mismatched values are retained in `legacy_reproduction.tsv`; they
must not be corrected by editing the archived files.

### 4. v3 BUS harness

Generalize `scripts/benchmark_v3_multimap.py` so it:

- accepts arbitrary sample/run IDs;
- requires the per-tree archived whitelist;
- accepts the original t2g explicitly;
- reads the source tree without modifying or copying it;
- corrects and sorts raw BUS using the archived chemistry whitelist;
- preserves CB-UMI molecule identity and ignores BUS read multiplicity as
  molecule mass;
- writes v3 unique, equal, and host-conservative metrics in one pass;
- writes an H5AD, count audit, per-virus and per-cell tables, resource metrics,
  hashes, status, command, and failure records to a separate output tree;
- writes temporary resolved BUS text under `$TMPDIR`;
- retains or hashes the resolved binary BUS and removes only generated scratch
  text after successful output validation.

All successful rows must satisfy:

`input_molecules = unique_molecules + ambiguous_molecules + unresolved_molecules`

`selected_matrix_mass = unique_molecules + allocated_ambiguous_mass`

Matrices must be finite, non-negative, and aligned to the declared barcode and
feature universes.

### 5. Tests

Use vertical red-green TDD slices for:

- inventory and frozen-cohort validation;
- technical-repeat grouping;
- chemistry whitelist discovery;
- legacy summary reconstruction;
- the known original-count `1` plus corrected-count `1` legacy call of `2`;
- v3 molecule conservation;
- barcode/feature mismatch rejection;
- finite/non-negative matrix validation;
- deterministic row/chunk ordering;
- failure retention and sanitized tracked outputs.

Run focused tests first, then existing multimap/rerun/validation tests, Ruff,
and the full repository suite.

### 6. Real-data pilots

Pilot 1:
`WS_SKN_KCL10525738/S1/L001/WS`.

The pilot must reconstruct the known legacy EBNA-2 `1 + 1 = 2` call, run the v3
layers twice, prove deterministic matrices and audits, and leave the source
tree byte-identical.

Pilot 2:
`SRR12682296`.

Run with its archived 10x v2 whitelist. Compare against the retained v3 EBV
baseline while documenting that the earlier standalone benchmark omitted
whitelist correction. Measure wall time, peak RSS, scratch use, and output size.

No full array is submitted until both pilots pass.

### 7. Full BUS execution

Run one SLURM task per technical tree with maximum concurrency four.

Default resource tiers:

- raw BUS below 4 GiB: 8 cores, 32 GiB RAM, 6 hours;
- raw BUS at least 4 GiB: 8 cores, 128 GiB RAM, 24 hours.

Require at least three times the raw BUS size in `$TMPDIR`. Each task emits a
validated result or an explicit failure containing the stage, command, exit
code, stderr path, attempt ID, and unchanged scientific parameters.

### 8. BUS aggregation

Create sanitized:

- `legacy_reproduction.tsv`;
- `run_metrics.tsv`;
- `virus_metrics.tsv`;
- `cell_concordance.tsv`;
- `technical_repeat_audit.tsv`;
- `failures.tsv`.

Report paired count changes, molecule partitions, expected-target and off-target
signals, call-set Jaccard similarity, called-cell barcode overlap, lane
concordance, runtime, and memory. Do not perform inferential BCC-versus-normal
testing.

### 9. End-to-end controls

Query ENA for SRR12682296-98 FASTQ URLs, MD5 values, and byte sizes; archive the
metadata response, download resumably, verify MD5, then calculate SHA-256.

Hash and use the retained SRR6825024-25 FASTQs. Both packages must consume the
same frozen FASTQ files and original combined index/t2g.

Fresh v2.2.0 runs use Emma's read-only environment and new output directories.
Compare processed/pseudoaligned reads, matrix dimensions, genes, viruses,
barcodes, and totals with the archived outputs.

Fresh v3 runs use a source-bundle-pinned full-workflow environment. Run
host-conservative as product primary and export unique/equal diagnostic layers.
Every output must pass `viralscan validate-run`.

Where exact reads are available, run the v3 evidence workflow for the expected
target and the largest non-target candidate. Evidence remains diagnostic.

### 10. Interpretation and freeze

The report must answer separately:

- Did fresh v2.2.0 reproduce its archive?
- What changed from legacy v2 to molecule-aware equal on identical BUS evidence?
- What changed under host-conservative allocation?
- What changed between the two full package stacks on identical FASTQs?
- Were expected control targets recovered?
- Which skin candidates disappeared, persisted, or repeated across lanes?

Mandatory claim rules:

- abundant EBV recovery is qualitative positive-control evidence;
- failure to recover HIV remains a sensitivity failure;
- persistent unexpected calls remain specificity concerns;
- loss of a skin call supports an artifact explanation but does not prove viral
  absence;
- no BCC association, diagnostic-performance, or general package-superiority
  claim is allowed.

Freeze source, environment, input, reference, raw-result, and aggregate hashes.
Perform a second-pass audit of cohort membership, denominators, layers, repeat
handling, failures, and path sanitization before marking the tracker complete.

## Validation commands

```bash
git status --short --branch
git rev-parse HEAD
PYTHONPATH=src python -m pytest tests/test_legacy_v2_v3.py -q
PYTHONPATH=src python -m pytest tests/test_multimap.py tests/test_multimapping.py tests/test_rerun_multimap.py tests/test_validation.py -q
PYTHONPATH=src python -m ruff check scripts/compare_legacy_v2_v3.py scripts/benchmark_v3_multimap.py tests/test_legacy_v2_v3.py
PYTHONPATH=src python scripts/compare_legacy_v2_v3.py inventory --source-root "$EMMA_RESULTS_ROOT" --raw-manifest "$RUN_ROOT/source_manifest.tsv" --tracked-manifest analysis/legacy_v2_v3/cohort_manifest.tsv
PYTHONPATH=src python scripts/compare_legacy_v2_v3.py validate --run-root "$RUN_ROOT"
viralscan validate-run "$CONTROL_V3_RUN"
python3 scripts/check_data_governance.py
python3 scripts/validate_v3_protocol.py
```

Run the full repository test command recorded in `PLAN.md` after the focused
suite passes.

## Stop conditions

Stop before outcome execution if:

- the source inventory is not exactly 44 technical rows and 42 logical inputs;
- any required BUS, EC, transcript, whitelist, t2g, H5AD, or summary artifact is
  missing without an explicit frozen failure row;
- the v2.2.0 or v3 source/environment/reference cannot be fingerprinted;
- the v3 full-workflow environment fails its required doctor profile;
- either pilot fails a count, alignment, determinism, scratch, or source-integrity
  gate;
- the relevant source bundle changes after freeze;
- access to missing FASTQs requires unapproved network or credentials;
- a destructive operation outside generated scratch/output paths is proposed;
- implementation would modify Emma's tree or unrelated dirty-worktree changes.

Stop before aggregation or interpretation if successful v3 rows violate count
invariants. Failed rows remain in the denominator and report.

## Required final artifacts

- updated `analysis/legacy_v2_v3/TRACKER.md`;
- frozen protocol and sanitized cohort manifest;
- raw ignored run manifest, commands, logs, hashes, status, and failures;
- validated BUS comparison tables;
- five-control end-to-end reproduction tables;
- diagnostic report with explicit claim boundaries;
- concise restart handoff for any incomplete rows.
