# Matched ViralScan and STARsolo Benchmark Plan

> Historical pre-v3 execution packet. It must not supply v3 scientific claims.

## Objective

Execute a reproducible 12-row benchmark comparing combined human+virus references against separate host-first/virus-second alignment for ViralScan and STARsolo across HHV-6B, EBV, and HSV-1 datasets, using one auditable human source release and one all-virus panel: Serratus plus expanded anellovirus.

## Current State

The previous execution pass created a benchmark harness and a blocked run packet under `benchmark_runs/reference_strategy_2026-06-27/`. Treat those generated files as evidence and reusable scaffolding, not as submission-ready truth. The current `reference_audit.tsv` shows missing STAR genome directories, a missing STAR-compatible all-virus genome FASTA, and missing kallisto indices. No benchmark SLURM jobs should be submitted until a full reference audit succeeds.

The worktree is dirty with unrelated pre-existing files. Preserve unrelated changes. Do not overwrite historical EBV-only STARsolo/ViralScan artifacts such as `starsolo_p22_6*`, existing `logs/`, or existing manuscript edits unless the validation gates below pass.

## Frozen Plan

1. Re-establish the execution packet and preflight state.
   - Read `reference_manifest.json`, `reference_audit.tsv`, `benchmark_runs/reference_strategy_2026-06-27/reference_audit.tsv`, and the current `src/viralscan/reference_strategy.py` harness before editing.
   - Confirm the blocked reference audit rows are still true or update them from current filesystem evidence.
   - Keep the previous generated run directory as historical evidence unless deliberately creating a new dated run directory.

2. Fix the benchmark harness before regenerating commands.
   - Update manifest path auditing so every manifest path field is audited by schema, including fields named `gtf`, `t2g`, `kallisto_index`, `genome_dir`, `genome_fasta`, `genome_gtf`, host FASTA/GTF, and accession/provenance tables.
   - Require audit success before `prepare_reference_strategy_benchmark.py` writes or refreshes `commands.jsonl` and SLURM scripts.
   - Require non-empty SHA256 for all file artifacts when running the final audit; `--no-sha256` may be used only for quick blocker triage.
   - Add machine-checkable manifest/audit fields for build commands, source release/provenance, feature counts, expected-virus presence for HHV-6B/EBV/HSV-1, anellovirus expected/fetched/missing counts, sizes, mtimes, and hashes.
   - Replace shell-`eval` command execution with argv-safe execution, or materialize numeric thread counts in the wrapper before command serialization.
   - Add dataset-specific STARsolo geometry for `10xv2`, `10xv3`, and `DROPSEQ` using the same CB/UMI geometry logic used by ViralScan host filtering.

3. Build or locate the canonical references from the manifest.
   - Use GRCh38 2024-A human genome FASTA and GTF as the single human source for STARsolo.
   - Use the matching GRCh38 2024-A transcriptome/cDNA source for kallisto/ViralScan, derived from or documented against the same Cell Ranger 2024-A release.
   - Build or locate a STAR-compatible all-virus genome FASTA/GTF for Serratus plus expanded anellovirus. Do not use transcriptome FASTA as input to `STAR --runMode genomeGenerate`.
   - Build STAR genome directories for human-only, all-virus-only, and combined human+all-virus references.
   - Build kallisto references for human-only, all-virus-only, and combined human+all-virus transcriptomes.
   - Reuse existing artifacts only if hashes, source paths, and panel counts prove they match the manifest.
   - Stop if a canonical all-virus genome FASTA/GTF cannot be produced or audited.

4. Strengthen reference and command validation.
   - Fail the audit if any required path is missing, if any SHA256 is empty in final mode, if panel counts are absent, or if HHV-6B/EBV/HSV-1 cannot be found by a full-file or indexed check.
   - Fail if any manifest path or generated command references single-virus or historical paths, including `fasta_split`, `split_gtf`, `GRCh38_EBV`, `starsolo_p22_6`, `transcriptomev3_noHuman`, `Epstein_Barr_virus_NC_007605`, `Human_herpesvirus_1_NC_001806`, or `Human_herpesvirus_6_NC_001664`.
   - Add tests that cover all manifest path fields, missing path rejection, hash/provenance requirements, all three target-virus presence checks, anellovirus counts, forbidden references, STAR second-pass command generation, and denominator invariants.

5. Regenerate a fresh dated benchmark run directory after the audit passes.
   - Use a new run directory name if the existing run directory contains stale blocked artifacts.
   - Copy the final `reference_manifest.json` into that run directory.
   - Write `reference_audit.tsv`, `commands.jsonl`, `run_reference_strategy_array.sh`, and an initial `run_status.tsv`.
   - Generate exactly 12 primary rows: for each of HHV-6B `SRR20710641`, EBV `SRR12682296`, and HSV-1 `SRR8315713`, run STARsolo combined, STARsolo two-step, ViralScan combined, and ViralScan two-step.
   - Do not mark command/SLURM generation complete until the command manifest references only existing audited references.

6. Implement full two-step execution semantics.
   - ViralScan two-step must run kallisto host filtering and then viral-only ViralScan/kb on the host-unmapped reads.
   - STARsolo two-step must run a host STARsolo/STAR pass that writes paired host-unmapped FASTQs, verify mate order and read counts, then run a second all-virus STARsolo pass on those recovered FASTQs.
   - Record both host-pass and virus-pass commands for two-step rows in the command manifest.
   - Validate R1/R2 order, paired unmapped FASTQ counts, CB/UMI geometry, whitelist behavior, `soloFeatures`, `soloCellFilter`, and multimapper settings before row status is `complete`.

7. Submit and monitor via SLURM only after all preflight gates pass.
   - Use absolute log paths under the benchmark run directory.
   - Submit with `sbatch --parsable` and record job IDs, array task IDs, row IDs, command hashes, log paths, environment/tool versions, start/end times, status, and failure reasons in machine-readable artifacts.
   - Do not run heavy benchmark rows in background shell sessions.
   - If any row fails, keep the failure in `run_status.tsv`, stop before manuscript updates, and do not mix partial results into claims.

8. Summarize results deterministically from run outputs.
   - Implement or update result parsing so `scripts/summarize_reference_strategy.py` can derive `results/reference_strategy_benchmark.tsv` from the run directory, not merely validate a pre-existing TSV.
   - Parse ViralScan outputs from `results/viral_summary.tsv`, per-cell viral tables, and relevant count layers.
   - Parse STARsolo outputs from `Solo.out/<feature>/raw` and/or `filtered`, `Summary.csv`, feature/barcode matrices, and the two-step viral pass outputs.
   - Include row identity, reference hash IDs, command/job IDs, `count_layer`, target-virus metrics, related off-target-virus metrics, fixed barcode-universe denominator, method-called-cell denominator, shared-anchor denominator, and explicit combined-vs-two-step delta columns or companion delta table.
   - Enforce fixed per-dataset barcode universes with missing barcodes filled as zero; two-step viral-only called-cell denominators must come from the host/combined pass or external anchor, never from viral-only called cells.

9. Validate before manuscript updates.
   - Run the focused tests and any added tests.
   - Run final reference audit with SHA256 enabled.
   - Validate the generated benchmark TSV and denominator invariants.
   - Run the forbidden-reference scan against the new benchmark run directory and final TSV.
   - Update `BENCHMARK_COMPARISON.md` and `docs/manuscript_draft.md` only if all 12 rows are complete, all hashes match the manifest, and validations pass.
   - Mark old EBV-only or older-panel results as historical/superseded; do not combine them with new comparative claims.

## Constraints and Non-Goals

- No alevin-fry comparison.
- No negative-control dataset.
- No single-virus references.
- No heavy background jobs; use SLURM.
- Do not claim broad STARsolo superiority or inferiority.
- Frame results as within-method reference-strategy deltas and cautious cross-method comparisons.
- Preserve unrelated dirty worktree changes.
- Stop before destructive operations, restricted network access, or credential use without explicit permission.

## Validation Commands

```bash
PYTHONPATH=src python -m pytest tests/test_reference_strategy_benchmark.py tests/test_reference_audit.py -q
PYTHONPATH=src python scripts/audit_reference_strategy.py --run-dir <benchmark_run_dir> --expected-panel serratus_plus_expanded_anellovirus --fail-on-single-virus
PYTHONPATH=src python scripts/prepare_reference_strategy_benchmark.py --run-dir <benchmark_run_dir> --manifest <benchmark_run_dir>/reference_manifest.json
PYTHONPATH=src python scripts/summarize_reference_strategy.py --run-dir <benchmark_run_dir> --out results/reference_strategy_benchmark.tsv --validate results/reference_strategy_benchmark.tsv
rg -n "fasta_split|split_gtf|GRCh38_EBV|starsolo_p22_6|transcriptomev3_noHuman|Epstein_Barr_virus_NC_007605|Human_herpesvirus_1_NC_001806|Human_herpesvirus_6_NC_001664" <benchmark_run_dir> results/reference_strategy_benchmark.tsv
```

## Stop Conditions

- Stop if required reference sources cannot be located or built from auditable provenance.
- Stop if final reference audit has any missing path, empty SHA256, missing build command/provenance, missing panel counts, or missing HHV-6B/EBV/HSV-1 presence.
- Stop if the human source cannot be matched across STARsolo and ViralScan references.
- Stop if the Serratus plus expanded-anellovirus panel cannot be built or audited.
- Stop if any benchmark command references a single-virus or historical reference path.
- Stop if STAR-compatible all-virus genome FASTA/GTF is missing.
- Stop if STARsolo two-step does not include both host-unmapped recovery and a second all-virus STARsolo pass.
- Stop before submitting SLURM jobs if command generation references non-existent references.
- Stop before manuscript updates if any of the 12 primary rows is failed, blocked, or incomplete.
