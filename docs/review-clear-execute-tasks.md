# ViralScan v2.2.0 versus v3 legacy diagnostic tasks

## LVC-00 — Freeze execution packet

- [x] Read `docs/review-clear-execute-plan.md`, this task list, `PLAN.md`, and `docs/handoffs/handoff-CFtSVg.md`.
- [x] Run `git status --short --branch` and record unrelated changes that must be preserved.
- [x] Create `analysis/legacy_v2_v3/TRACKER.md` with task IDs LVC-00 through LVC-14 and status definitions.
- [x] Create `analysis/legacy_v2_v3/protocol.yaml` with the diagnostic-only claim boundary, frozen methods, cohort rules, endpoints, and failure policy.
- [x] Record in the tracker that historical outputs are excluded from confirmatory v3 validation and release claims.

## LVC-01 — Inventory source artifacts

- [x] Implement the `inventory` subcommand in `scripts/compare_legacy_v2_v3.py`.
- [x] Add the first failing inventory test, implement the minimum behavior, and repeat per TDD.
- [x] Inventory all legacy result trees into an ignored raw manifest.
- [x] Emit a sanitized relative-path `analysis/legacy_v2_v3/cohort_manifest.tsv`.
- [x] Verify exactly 44 technical rows, 42 logical inputs, six 10x v2 whitelists, and 38 10x v3 whitelists.
- [x] Mark the two SRR6825025 outputs as a technical-repeat group.
- [x] Mark the two KCL10525740 S3 L002-input outputs as a technical-repeat group.
- [x] Hash every required BUS, EC, transcript, whitelist, t2g, H5AD, summary, and reference artifact.
- [x] Stop and record a blocker if frozen cohort counts or required artifacts do not match.

## LVC-02 — Freeze source and environments

- [x] Record Emma's ViralScan 2.2.0 package metadata and installed code hashes.
- [x] Record v2.2.0 kb-python, kallisto, and bustools versions.
- [x] Hash the original index, t2g, and transcriptome.
- [x] Create a source bundle of the current fork and record repository HEAD, relevant diff hash, and bundle SHA-256.
- [x] Build a full-workflow v3 environment under the ignored run root.
- [x] Install the frozen source bundle into that environment.
- [x] Archive explicit package/tool versions and the built wheel hash.
- [x] Run and archive `viralscan doctor --profile full --json`.
- [x] Freeze the environment and source identifiers before outcome-producing runs.

## LVC-03 — Implement comparison behavior

- [x] Implement the `run-row` subcommand using public benchmark functions.
- [x] Implement v2.2.0 H5AD/summary reconstruction with tolerance `1e-6`.
- [x] Generalize `scripts/benchmark_v3_multimap.py` for arbitrary run IDs.
- [x] Require and apply the archived chemistry whitelist.
- [x] Require the original t2g explicitly.
- [x] Keep the source result tree read-only and write all products to a separate output root.
- [x] Emit v3 unique, molecule-aware equal, and host-conservative metrics in one BUS pass.
- [x] Emit count audits, H5AD, per-virus, per-cell, resource, command, status, hash, and failure artifacts.
- [x] Implement `summarize` without silently excluding failures.
- [x] Implement `validate` for manifests, output schemas, count invariants, matrices, paths, and hashes.
- [ ] Ensure tracked outputs contain no institutional absolute paths.

## LVC-04 — Verify the harness

- [x] Add tests for the exact 44/42 frozen cohort contract.
- [x] Add tests for technical-repeat grouping.
- [x] Add tests for 10x v2/v3 whitelist discovery.
- [x] Add a test reconstructing the legacy `1 + 1 = 2` call.
- [x] Add tests for legacy gene/virus summary agreement.
- [x] Add tests for v3 molecule and selected-matrix conservation.
- [x] Add tests rejecting barcode or feature misalignment.
- [x] Add tests rejecting negative or non-finite matrices.
- [x] Add tests for row-order and chunk-order invariance.
- [x] Add tests retaining failed and unsupported rows.
- [x] Run `PYTHONPATH=src python -m pytest tests/test_legacy_v2_v3.py -q`.
- [x] Run the existing multimap, rerun, and validation focused tests.
- [x] Run Ruff on the changed scripts and tests.
- [x] Run the full repository suite recorded in `PLAN.md`.

## LVC-05 — Minimal skin pilot

- [x] Inventory `WS_SKN_KCL10525738/S1/L001/WS` and verify all hashes.
- [x] Reconstruct the legacy EBNA-2 original `1`, corrected `1`, combined `2` result.
- [x] Run v3 BUS processing with the archived 10x v3 whitelist.
- [x] Run the same v3 row a second time with unchanged inputs.
- [x] Compare output matrices, audits, hashes, and metrics for determinism.
- [x] Verify molecule conservation, finite/non-negative matrices, and barcode/feature alignment.
- [x] Verify the source result tree hashes remain unchanged.
- [x] Record pilot commands, resource use, and exit-gate evidence in the tracker.

## LVC-06 — EBV scale pilot

- [x] Inventory SRR12682296 and verify its archived 10x v2 whitelist.
- [x] Run the full BUS harness on SRR12682296.
- [x] Compare the result with the retained v3 EBV baseline and document the earlier whitelist omission.
- [x] Run row/chunk invariance checks.
- [x] Record wall time, peak RSS, scratch use, and output size.
- [x] Confirm or revise the frozen SLURM resource tiers before full submission.
- [x] Stop if the scale pilot violates any invariant or scratch limit.

## LVC-07 — Execute all BUS rows

- [x] Generate a 44-row SLURM manifest from the frozen source manifest.
- [x] Separate tasks into below-4-GiB and at-least-4-GiB resource tiers.
- [x] Verify `$TMPDIR` capacity is at least three times each raw BUS size.
- [x] Submit bounded arrays with maximum concurrency four.
- [x] Record job IDs, task IDs, run IDs, command hashes, log paths, source hash, and environment hash.
- [x] Monitor all tasks to terminal state.
- [x] Preserve each failure with stage, exit code, stderr, attempt ID, and unchanged parameters.
- [x] Do not retry with modified scientific settings.

## LVC-08 — Validate BUS outputs

- [x] Run the comparison `validate` command over the complete run root.
- [x] Verify legacy reconstructed values against all reported summaries within `1e-6`.
- [x] Verify molecule conservation on every successful v3 row.
- [x] Verify selected-matrix conservation on every successful v3 row.
- [x] Verify matrices are finite, non-negative, and barcode/feature aligned.
- [x] Compare both technical-repeat groups without pooling them as independent samples.
- [x] Investigate every failed invariant and record the cause.
- [x] Retain unresolved failures in the final denominator and tracker.

## LVC-09 — Aggregate BUS results

- [x] Generate `analysis/legacy_v2_v3/legacy_reproduction.tsv`.
- [x] Generate `analysis/legacy_v2_v3/run_metrics.tsv`.
- [x] Generate `analysis/legacy_v2_v3/virus_metrics.tsv`.
- [x] Generate `analysis/legacy_v2_v3/cell_concordance.tsv`.
- [x] Generate `analysis/legacy_v2_v3/technical_repeat_audit.tsv`.
- [x] Generate `analysis/legacy_v2_v3/failures.tsv`.
- [x] Check paired count changes, molecule partitions, call-set Jaccard values, barcode overlaps, lane concordance, runtime, and memory.
- [x] Confirm no inferential BCC-versus-normal test was introduced.
- [x] Run path-sanitization and data-governance checks.

## LVC-10 — Acquire five control inputs

- [x] Query and archive ENA file metadata for SRR12682296, SRR12682297, and SRR12682298.
- [x] Download paired EBV FASTQs resumably after network approval.
- [x] Verify ENA MD5 values and calculate SHA-256 values.
- [x] Hash the retained SRR6825024 and SRR6825025 FASTQ pairs.
- [x] Validate read pairing, file integrity, chemistry, and immutable input paths.
- [x] Freeze one control-input manifest consumed by both packages.

## LVC-11 — Fresh v2.2.0 control reproduction

- [ ] Run all five controls with Emma's read-only v2.2.0 environment into new output directories.
- [ ] Use the frozen original index, t2g, chemistry, and identical FASTQ files.
- [ ] Record commands, tool versions, runtimes, memory, and output hashes.
- [ ] Compare processed and pseudoaligned reads with archived runs.
- [ ] Compare matrix dimensions, reported genes/viruses, totals, and barcodes.
- [ ] Classify every mismatch as input, dependency, nondeterminism, reference, package, or unresolved drift.
- [ ] Retain reproduction failures instead of modifying the archive.

## LVC-12 — Fresh v3 control execution

- [ ] Run all five controls with the frozen v3 source bundle and full-workflow environment.
- [ ] Use the same frozen FASTQs and original combined reference.
- [ ] Run host-conservative as product primary and export unique/equal diagnostic layers.
- [ ] Run `viralscan validate-run` for every control.
- [ ] Record commands, source/environment/reference hashes, runtime, memory, storage, and outputs.
- [ ] Run exact-read evidence for the expected target and largest non-target candidate where reads exist.
- [ ] Retain every failed control row and evidence step.

## LVC-13 — Interpret results

- [ ] Separate fresh-v2 reproduction findings from v2-to-v3 counting changes.
- [ ] Separate identical-BUS counting effects from full product-stack effects.
- [ ] Report each EBV and HIV control individually.
- [ ] Report skin results descriptively by biological sample and technical run.
- [ ] State whether each skin candidate disappeared, persisted, or repeated across lanes.
- [ ] State that EBV recovery is qualitative positive-control evidence.
- [ ] State that failure to recover HIV is a sensitivity failure.
- [ ] Keep persistent unexpected calls as specificity concerns.
- [ ] Avoid claims of viral absence, BCC association, calibrated diagnostic performance, or package superiority.
- [x] Draft `analysis/legacy_v2_v3/REPORT.md` from validated tables only.

## LVC-14 — Audit, freeze, and hand off

- [ ] Re-audit cohort membership, hashes, versions, methods, denominators, layers, repeat handling, and failed rows.
- [ ] Confirm there were no outcome-driven exclusions or threshold changes.
- [ ] Confirm tracked artifacts contain no private data or institutional absolute paths.
- [ ] Freeze source, environment, reference, input, raw-result, aggregate, and report hashes.
- [ ] Run `python3 scripts/check_data_governance.py`.
- [ ] Run `python3 scripts/validate_v3_protocol.py`.
- [ ] Mark tracker tasks complete only where exit evidence exists.
- [ ] Write a concise restart handoff for incomplete or failed rows.
- [ ] State explicitly that this diagnostic does not close v3 truth-panel, comparator, release, or publication gates.
