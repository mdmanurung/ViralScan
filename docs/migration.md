# Migrating from earlier versions

ViralScan v3 does **not** migrate legacy outputs. There is no converter, and
none is planned. The v3 counting contract, identity table and evidence tiers
are not reconstructable from pre-v3 result files.

- **Rebuild only.** Rebuild the reference with `viralscan build-ref` and rerun
  v3 from the original FASTQ files.
- **Do not compare or carry over old values.** Result values from earlier
  versions (counts, percentages, virus lists) must not be compared with v3
  output or copied into v3 tables, reports or manuscripts. Differences reflect
  changed definitions, not just changed numbers.
- **Old run directories** are not accepted by `viralscan validate-run` as v3
  runs; keep them as archives and write new results to a fresh `--output`.

See the [output reference](output_reference.md) for the v3 schema and the
[FAQ](faq.md) for which comparisons across strategies are valid.

## Command-line changes

Quantification is now the `run` subcommand. A bare `viralscan [OPTIONS]` still runs it, with a
deprecation warning. The spellings below were removed without aliases; each one now fails with
the replacement named. Configuration keys, `config.yaml` and `run_fingerprint` are unchanged.

| Subcommand | Old | New |
|---|---|---|
| (quantification) | `viralscan -s1 R1 -s2 R2 ...` | `viralscan run --r1 R1 --r2 R2 ...` |
| `run`, `check-whitelist`, `check-chemistry` | `-s1` / `--sample1`, `-s2` / `--sample2` | `--r1`, `--r2` |
| `run` | `-ref`, `-reference`, `-gtf`, `-fasta`, `-f1`, `-umap`, `-acc`, `-mm` | `--reference`, `--gtf`, `--fasta`, `--f1`, `--umap`, `--ncbi-accession`, `--multimapping` |
| `run` | `-v` | `--visual` / `--no-visual` |
| `run` | `--reference True`, `--visual False`, `--multimapping False` (value form) | `--reference`, `--no-visual`, `--no-multimapping` |
| `run` | `--data-cache-dir` | `--cache-dir` |
| `run` | `--cell-types CSV` | `--cell-type-labels CSV` |
| `run` | `--hostresponse-control-mito` | `--hostresponse-mito-control` |
| `hostresponse` | `--hostresponse-min-umi`, `--detection-threshold` | `--min-umi` |
| `hostresponse` | `-o RUN_DIR` (an input run directory) | `--run-dir RUN_DIR`; `--output DIR` now optionally names where hostresponse writes (default `<run-dir>/posthoc/hostresponse`) |
| `hostresponse` | `--no-use-hvg` (still valid) | also `--use-hvg` |
| `validate-run` | positional `RUN_DIR` | `--run-dir RUN_DIR` |
| `validate-run` | `--json-output PATH` | `--json PATH` |
| `doctor` | `--json` (switch; the report was always printed) | `--json PATH` |

### Post-hoc outputs move to `<run-dir>/posthoc/` (layout version 2)

Post-hoc verbs no longer write into, or require an output outside of, the run's own files.
`--output` still overrides every default.

| Subcommand | Before | Now |
|---|---|---|
| `evidence` | `--output` required | optional; default `<run-dir>/posthoc/evidence/<virus>` |
| `rerun-multimap` | `-o` required, outside the run | optional; default `<run-dir>/posthoc/multimap_<method>`; still must be empty or absent |
| `rerun-programs` | wrote `results/gene_program_*.tsv`, `log/gene_programs.done` and re-stamped `run_complete.json` in the run | writes `<run-dir>/posthoc/programs/results/gene_program_*.tsv`; the run is untouched |
| `hostresponse` | `<run-dir>/hostresponse` and re-stamped `run_complete.json` | `<run-dir>/posthoc/hostresponse`; marker untouched |
| `evidence --review-candidates` | `results/evidence_review.tsv` in the run | `evidence_review.tsv` in the evidence output directory |

`run_manifest.json` now records `layout_version: 2`. `validate-run` accepts versions 1 and 2 (and
manifests without the field) and ignores `posthoc/`; a newer version is an error. Version 1 runs keep
working with every verb. Anything reading `<run>/hostresponse`, `results/gene_program_*.tsv` or
`results/evidence_review.tsv` after a post-hoc command must read the `posthoc/` path (or pass `--output`).
After `rerun-programs`, a later `--resume` runs the gene-programme rule itself, since the run's
`log/gene_programs.done` is no longer written.

Long options no longer match by prefix (`--multi` is an error, not `--multimapping`).
`build-ref` gains `-c` for `--cores`, and `doctor` and `validate-run` gain `--verbose`/`--quiet`.
