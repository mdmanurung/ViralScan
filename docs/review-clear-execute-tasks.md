# Matched ViralScan and STARsolo Benchmark Tasks

> Historical pre-v3 execution record. The v3 tracker is
> `docs/plans/2026-07-22-viralscan-v3-tasks.md`.

- [x] Read `docs/review-clear-execute-plan.md`, this task list, `reference_manifest.json`, root `reference_audit.tsv`, and `benchmark_runs/reference_strategy_2026-06-27/reference_audit.tsv`.
- [x] Confirm the current blocked reference rows are still true, or refresh blocker evidence from the filesystem.
- [x] Inspect the dirty worktree and identify unrelated pre-existing changes to preserve.
- [x] Update manifest path auditing to cover every manifest path field by schema, including `gtf`, `t2g`, `kallisto_index`, `genome_dir`, `genome_fasta`, and `genome_gtf`.
- [x] Add final-audit requirements for non-empty SHA256, build commands, source provenance, feature counts, HHV-6B/EBV/HSV-1 presence, and anellovirus expected/fetched/missing counts.
- [x] Make `prepare_reference_strategy_benchmark.py` refuse to write command/SLURM artifacts unless the full reference audit passes.
- [x] Replace shell-`eval` execution or make command execution argv-safe with explicit numeric thread counts.
- [x] Add dataset-specific STARsolo CB/UMI geometry and chemistry validation for `10xv2`, `10xv3`, and `DROPSEQ`.
- [x] Add tests for manifest path coverage, missing reference rejection, SHA/provenance requirements, all target-virus presence checks, anellovirus counts, forbidden references, STAR two-step command structure, and denominator invariants.
- [x] Build or locate GRCh38 2024-A human genome FASTA/GTF and matching kallisto transcriptome/cDNA source.
- [x] Package Serratus plus expanded-anellovirus STAR-compatible all-virus genome FASTA/GTF under `references/starsolo/all_virus_serratus_plus_anellovirus/`.
- [x] Package combined GRCh38 2024-A plus all-virus genome FASTA/GTF under `references/starsolo/combined_GRCh38_2024A_serratus_plus_anellovirus/`.
- [x] Resolve HHV-6B genome provenance for STARsolo by adding 97 `HUM_HERP6B` pseudo-contigs from the existing ViralScan/kallisto transcriptome to the packaged STAR FASTA/GTF.
- [x] Build STAR human-only, all-virus-only, and combined human+all-virus genome directories; all-virus job `25102653`, combined job `25102654`, and human job `25102710` completed.
- [x] Locate kallisto human-only and combined human+virus Serratus plus expanded-anellovirus transcriptome references from previous ViralScan runs.
- [x] Update `reference_manifest.json` with final source paths, existing index paths, build commands, anellovirus counts, and provenance; hashes, mtimes, sizes, and feature counts are recorded in `reference_audit.tsv`.
- [x] Run final reference audit with SHA256 enabled and zero missing paths.
- [x] Create a fresh dated benchmark run directory if the existing `benchmark_runs/reference_strategy_2026-06-27/` packet is stale.
- [x] Copy the final `reference_manifest.json` into the benchmark run directory.
- [x] Regenerate `reference_audit.tsv`, `commands.jsonl`, `run_reference_strategy_array.sh`, and initial `run_status.tsv` in the benchmark run directory.
- [x] Verify `commands.jsonl` contains exactly 12 primary rows and no forbidden single-virus or historical references.
- [x] Verify STARsolo two-step rows include both host-unmapped recovery and a second all-virus STARsolo pass.
- [x] Verify ViralScan two-step rows use kallisto host filtering and viral-only ViralScan/kb on host-unmapped reads.
- [x] Submit benchmark rows via SLURM only after reference and command preflight gates pass; primary array `25102689` and STAR retry array `25102744` are recorded.
- [x] Record `sbatch --parsable` job IDs, array task IDs, row IDs, command hashes, log paths, environment/tool versions, start/end times, status, and failure reasons where available.
- [x] Monitor all 12 rows and update `run_status.tsv` from machine-readable status evidence; current live rows remain marked `running_or_incomplete`.
- [x] Implement status-only result summarization into `results/reference_strategy_benchmark.tsv`; deterministic metric extraction remains intentionally unimplemented for incomplete/unparsed rows.
- [ ] Parse ViralScan `viral_summary.tsv`, per-cell viral tables, and count layers.
- [ ] Parse STARsolo matrices, `Summary.csv`, raw/filtered barcode universes, and two-step viral pass outputs.
- [ ] Add target-virus metrics, related off-target metrics, reference hash IDs, command/job IDs, `count_layer`, fixed barcode denominator, method-called-cell denominator, shared-anchor denominator, and combined-vs-two-step deltas.
- [ ] Validate fixed per-dataset barcode universes with missing barcodes filled as zero.
- [ ] Validate two-step viral-only denominators come from host/combined/external anchors, not viral-only called cells.
- [x] Run `PYTHONPATH=src python -m pytest tests/test_reference_strategy_benchmark.py tests/test_reference_audit.py -q`.
- [x] Run `PYTHONPATH=src python scripts/audit_reference_strategy.py --run-dir <benchmark_run_dir> --expected-panel serratus_plus_expanded_anellovirus --fail-on-single-virus`.
- [x] Run `PYTHONPATH=src python scripts/prepare_reference_strategy_benchmark.py --run-dir <benchmark_run_dir> --manifest <benchmark_run_dir>/reference_manifest.json` and verify it refuses to write command/SLURM artifacts while the audit is incomplete.
- [x] Run `PYTHONPATH=src python scripts/summarize_reference_strategy.py --run-dir <benchmark_run_dir> --out results/reference_strategy_benchmark.tsv --validate results/reference_strategy_benchmark.tsv`; validation is structurally successful but exits nonzero because 0/12 rows are finalized as complete metrics.
- [x] Run the forbidden-reference `rg` check against the benchmark run directory and `results/reference_strategy_benchmark.tsv`.
- [ ] Update `BENCHMARK_COMPARISON.md` and `docs/manuscript_draft.md` only if all 12 rows complete and all validation gates pass.
