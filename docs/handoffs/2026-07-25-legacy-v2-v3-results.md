# ViralScan legacy v2 versus v3 results handoff

Date: 2026-07-25

Branch: `codex/viralscan-v3`

Baseline commit: `a8f32876e6ec4a5d8cf29b21aaf325ccaa50e64f`

## Completed

- All 44 frozen identical-BUS comparison rows completed successfully.
- Full-root validation passed with 44 successes and no errors.
- Sanitized aggregate, cell-concordance, repeat-audit, candidate-persistence,
  validation, and report artifacts are under `analysis/legacy_v2_v3/`.
- The focused final analysis tests, Ruff, data governance, protocol
  validation, and diff check pass.
- No archived source or result was modified.

Primary interpretation: archived v2 values reproduce internally; v3 uses a
different molecule-aware count contract; EBV is recovered qualitatively in all
three positive controls; HIV recovery is not demonstrated because the frozen
v3 virus universe lacks an HIV-labelled feature; unexpected low-mass skin
entries remain specificity concerns.

## Active work

The five-control input gate is complete. All five pairs passed full-stream
structure, paired-identifier, chemistry, byte-count, MD5, and SHA-256 checks.
The sanitized manifest is `analysis/legacy_v2_v3/control_inputs.tsv`; the
executable copy remains under the ignored run tree.

Fresh-run attempt 1 arrays `25331024` and `25331026` failed before ViralScan
execution because Slurm copied the submitted shell script to its spool and the
script derived the packet root from that copy. The failure logs are retained
and no scientific output/status file was created.

Attempt 2 uses a tested explicit packet root and a byte-frozen, self-checking
runner packet. Array `25331035` contains the six EBV v2/v3 tasks with maximum
concurrency four. Dependency-gated array `25331037` contains the four HIV v2/v3
tasks. Both stacks consume identical per-sample FASTQ paths, sizes, hashes,
original index/t2g, 10x-v2 whitelist, chemistry, and eight-core policy.

Resume with:

```bash
squeue -j 25331035,25331037 -o '%.18i %.9T %.10M %.6D %R'
sacct -j 25331035,25331037 --format=JobID,State,ExitCode,Elapsed,MaxRSS
find benchmark_runs/legacy_v2_v3/fresh_control_packet_attempt2 -maxdepth 2 -type f -print
```

Do not mark a task complete from submission or output existence. Require
terminal exit evidence, nonempty expected artifacts, and `viralscan
validate-run` success for every v3 row. Retain failed rows without reuse or
overwrite.

## Next validation

1. Audit all ten terminal scheduler/status records.
2. Compare processed reads, pseudoalignment, barcodes, matrices, expected
   targets, and exact-read evidence only for validated outputs.
3. Freeze output hashes and update `REPORT.md` without changing thresholds.
4. Keep fresh-v2 reproduction findings separate from v2-to-v3 count-contract
   differences.

The current diagnostic does not close truth-panel, comparator, calibration,
release, publication, or package-superiority gates.
