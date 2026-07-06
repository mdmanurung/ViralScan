# Reference Strategy Run Packet: 2026-06-28 fresh12

This packet preserves lightweight, clean-clone-visible provenance for the 12-row
reference-strategy benchmark attempt (combined vs two_step × ViralScan vs STARsolo
× three target-virus SRA runs). The bulky run tree lives under the gitignored
`benchmark_runs/reference_strategy_2026-06-28_fresh12/`; only the small provenance
files are tracked here.

Tracked files:

- `commands.jsonl` — command manifest used by the SLURM array (one JSON object per row).
- `fastq_audit.tsv` — raw-byte FASTQ audit with sizes, hashes, record counts, and source URLs.
- `reference_audit.tsv` — reference artifact checksums and feature counts.
- `failure_summary.tsv` — concise status of the 8 rows that did not produce final
  ViralScan/STARsolo outputs, derived row-by-row from
  `results/reference_strategy_benchmark.tsv` (`row_id`, `status`, `slurm_job_id`, `reason`).

Row status (2026-07-02): 4 complete, 3 blocked, 2 incomplete, 2 running_or_incomplete,
1 failed. See `../../REFERENCE_STRATEGY_BENCHMARK.md` for the aggregate table and the
two valid comparisons the 4 complete rows support.

Bulky files intentionally remain ignored under `benchmark_runs/`: raw logs,
intermediate matrices, STAR indices, kallisto outputs, and FASTQs.

**Publication status:** this run packet is provenance only. It is *not* a complete
benchmark result and must not be cited as final performance evidence. The blocked
rows failed because `kallisto`/`kb` were absent on the compute node, and one EBV
STARsolo row failed on a FASTQ quality-string length error — both operational, not
methodological. Reintroduce this benchmark into manuscript claims only after all 12
rows complete and `scripts/summarize_reference_strategy.py` exits successfully.
