# ViralScan legacy v2 versus v3 diagnostic

This directory contains the small, tracked control plane for a diagnostic
comparison of archived ViralScan 2.2.0 results with the current v3
molecule-counting implementation.

The comparison is intentionally historical and developmental. The archived
outcomes have already been inspected, so they are not confirmatory holdout
evidence and cannot close v3 truth-panel, comparator, calibration, release,
publication, or package-superiority gates.

Tracked artifacts are limited to relative identifiers, hashes, the frozen
protocol, the tracker, sanitized manifests, small result tables, and the final
diagnostic report. Raw manifests containing absolute paths, large BUS/H5AD
artifacts, FASTQs, environments, scheduler logs, and run products belong under
the ignored `benchmark_runs/` tree.

Start with:

- [`protocol.yaml`](protocol.yaml) for the frozen scientific and execution
  contract;
- [`TRACKER.md`](TRACKER.md) for evidence-backed task status;
- [`REPORT.md`](REPORT.md) for the validated identical-BUS conclusions;
- `cohort_manifest.tsv` for the sanitized frozen denominator;
- `validation.json` for the final 44-row validation result; and
- the TSV files in this directory for legacy reproduction, run/virus metrics,
  cell concordance, technical-repeat audit, and candidate persistence.
