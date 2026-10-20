# CPU-budget qualification status — 2026-10-20

The CPU-budget candidate bounds native worker requests and Snakemake reservations.
Implementation approval does not resolve the raw-byte compatibility gate below.
No runtime code changed during this evidence reconciliation.

## Completed local checks

- Full default suite: 5841 passed, 15 skipped, 163 deselected.
- Native/integration suite: 118 passed, zero skipped, 5901 deselected.
- Review fixes committed in `3985b7a`: CLI-to-native worker plumbing, rule-expression
  reservations tested above the scheduler clamp, and runtime numerical-library pins.
- Original DAG matrix: 48 paired cases / 96 dry runs. Overnight matrix: 160 cases,
  covering 1/6/8/32 cores, whitelist branches, optional stages, reruns and overrides.
- The overnight native matrix matched scientific bytes in all completed arms.
  Operational metadata, logged worker arguments and harness core inputs differ.
  A/A variability in an operational file does not authorize every A/B value change.

## Unresolved raw-index discrepancy

The earlier 42-arm native matrix retained a one-core `_finish_reference_build`
index with 55 changed bytes, beginning at zero-based offset 76. Both files are
1750 bytes; prepared FASTA, cDNA, t2g, inspection and reproducibility files match.

- Raw-identical baseline A/A index SHA-256:
  `abe8a1be12682f5baec6d3edfadc1c284bea03b74ad18aaae6219bb195108595`
- Retained one-core candidate index SHA-256:
  `77468ce78a6029d016c43e6f0e6808a9cdd4b425d6031375d89b0a8da8f75d2d`

The overnight one-core index matches the baseline. Forty further direct native
builds using the retained cDNA (20 each at one/eight workers, alternating) also
match the baseline. These repeated matches do not explain the earlier difference
or prove that the worker count caused it. The failed comparison remains evidence.

Four real downstream count runs consumed the retained baseline/differing indexes:
baseline A/A at eight workers and candidate A/A at one worker. Complete output
file sets and scientific bytes match; only `kb_info.json` and `run_info.json`
differ. This fixture contains two processed, pseudoaligned reads and one called
barcode with one counted molecule. It does not establish production equivalence
or waive raw-index identity. No bytes, goldens or comparators were normalized.

## Delivery boundary

A draft PR may run CI and expose the change for review. Merge remains blocked by
the unresolved discrepancy and the user's existing requirement for a separate
changed-byte compatibility decision. Production G02 stage attribution, whole-run
performance and sustained CPU utilization remain unmeasured. No SLURM submission,
download, scientific setting change or retained-data deletion was performed here.

Detailed monitoring receipts live in the shared research checkout at
`analysis/computational_bottlenecks_2026-10-17/followups/20261019-performance-monitor/`.
They identify both historical native matrices, exact raw hashes, downstream
commands and the 40-repeat native probe. The original failed artifacts are retained.

## PR CI fixture correction

The first PR CI run found one failing CLI-plumbing test on workers without Rscript.
Its default emptyDrops preflight exited before the mocked reference builder.
The failure was reproduced locally with Rscript discovery disabled. The fixture
now explicitly selects knee calling and checks the successful stop-hook exit,
retaining the assertion that the reference builder receives two cores. The affected
CLI/worker suites pass 184 tests (16 integration cases deselected) with Rscript
hidden after the correction; Ruff and formatting pass. Production
code, dependency checks and scientific defaults are unchanged. Final CI results
are recorded on the PR for the corrected revision.
