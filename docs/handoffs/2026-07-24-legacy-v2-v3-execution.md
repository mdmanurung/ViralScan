# Fresh Session Handoff: execute the ViralScan legacy diagnostic

## Objective

Implement and execute the frozen comparison between Emma Vonk's archived
ViralScan 2.2.0 outputs and the current v3 fork, covering all retained BUS
artifacts and five end-to-end public controls.

## Repository

The repository root containing this handoff.

Revised Plan Path: `docs/review-clear-execute-plan.md`

Task List Path: `docs/review-clear-execute-tasks.md`

## Frozen Plan

1. Freeze a diagnostic-only tracker and protocol.
2. Inventory exactly 44 technical result trees representing 42 logical inputs.
3. Fingerprint Emma's read-only v2.2.0 environment, the v3 source/environment,
   and the original combined reference.
4. Implement a manifest-driven inventory/run/summarize/validate CLI and
   generalize the existing BUS benchmark to require archived whitelists.
5. Verify behavior through vertical TDD slices and focused/full tests.
6. Pass the known low-count skin pilot and the SRR12682296 scale pilot.
7. Run and validate all 44 BUS rows through bounded SLURM arrays.
8. Aggregate legacy, v3 unique, v3 equal, and v3 host-conservative results.
9. Freeze matching FASTQs and rerun the three EBV and two HIV controls through
   fresh v2.2.0 and v3 workflows.
10. Audit, interpret conservatively, freeze hashes, and leave a restart record.

## Constraints and Non-Goals

- Preserve existing `.living/**` and `docs/handoffs/` changes.
- Never modify Emma's results or installed environment.
- Keep large artifacts and absolute paths under ignored `benchmark_runs/`.
- Do not use legacy outputs as confirmatory v3 evidence.
- Do not change thresholds, methods, features, barcodes, or inclusion after
  outcomes are inspected.
- Do not make viral-absence, BCC-association, calibrated diagnostic, release,
  publication, or package-superiority claims.
- Do not add EM or unique-weighted outcome arms.
- Do not run heavy work outside SLURM.

## Validation Commands

- `git status --short --branch`
- `git rev-parse HEAD`
- `PYTHONPATH=src python -m pytest tests/test_legacy_v2_v3.py -q`
- `PYTHONPATH=src python -m pytest tests/test_multimap.py tests/test_multimapping.py tests/test_rerun_multimap.py tests/test_validation.py -q`
- `PYTHONPATH=src python -m ruff check scripts/compare_legacy_v2_v3.py scripts/benchmark_v3_multimap.py tests/test_legacy_v2_v3.py`
- `PYTHONPATH=src python scripts/compare_legacy_v2_v3.py validate --run-root "$RUN_ROOT"`
- `viralscan validate-run "$CONTROL_V3_RUN"`
- `python3 scripts/check_data_governance.py`
- `python3 scripts/validate_v3_protocol.py`
- the full repository test command recorded in `PLAN.md`

## Stop Conditions

- Stop if repository drift affects the frozen plan or relevant source.
- Stop if the source inventory is not exactly 44 technical rows and 42 logical
  inputs.
- Stop if required artifacts, hashes, whitelists, reference, environment, or
  permissions are missing.
- Stop if either real-data pilot fails its invariants or source-integrity gate.
- Stop before network access without approval.
- Stop before destructive operations outside generated scratch/output paths.
- Stop before modifying Emma's tree or unrelated dirty-worktree changes.
- Stop before aggregation if any successful row violates v3 count invariants.

## Required Artifacts

- `analysis/legacy_v2_v3/TRACKER.md`
- frozen protocol and sanitized cohort manifest
- ignored raw run manifest, commands, statuses, logs, hashes, and failures
- validated BUS comparison tables
- five-control end-to-end reproduction tables
- `analysis/legacy_v2_v3/REPORT.md`
- a concise restart handoff for incomplete rows

## Fresh Executor Instruction

Read the revised plan and task list completely before editing. Execute tasks in
dependency order, use non-overlapping implementation subagents where useful,
mark a checkbox only after validation, and stop at the listed gates.
