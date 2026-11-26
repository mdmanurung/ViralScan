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

## Migrating from 3.0.0rc1

`3.0.0rc1` was never published. If you ran a development build of 3.0, these changes apply
(details in the [changelog](changelog.md)).

- **Finish or restart runs begun by another version.** `--resume` refuses a run made by a different
  ViralScan version, and says which version made it. Start the run again in a new `--output`, or
  finish it with the version that began it. `validate-run` still reads the finished output.
- **Pass `--strand` for 5' and unconfirmed chemistries.** `run` without `--strand` exits 2 unless
  tchem confirms a 3' chemistry. Use `--strand forward|reverse|unstranded|auto`, or `--kb-params auto`
  to let the built-in detector, tchem and the strand pilot choose and stop on disagreement.
- **The infected-cell rule is stricter.** `--infected-min-molecules` defaults to 2.0, so `pct_infected`
  and the infected-cell counts fall for cells with a single molecule. To reproduce an earlier run, pass
  `--infected-min-molecules 1`. The new default has not been validated against the truth panel.
- **`--host-filter starsolo` uses the barcode on-list** for named 10x chemistries. Rerun a host-filtered
  run to get the new counts.
- **`--cell-calling none` and `knee` cannot be combined with `--host-index`.** Use `emptydrops` or
  `external`.
- **Post-hoc verbs take the run root.** `evidence`, `rerun-programs`, `hostresponse` and `host-loss`
  used to take one sample directory as `--run-dir`. Pass the run root (the `-o` directory of `run`)
  as `--run-dir`, and `--sample NAME` when it has several samples. The old form still works, prints a
  `viralscan: warning:` line that names the replacement, and is removed in 3.1.0. A root without
  `run_complete.json` is refused. Outputs and their default locations (`<sample>/posthoc/<verb>/`)
  are unchanged.
- **Four `run` options that changed nothing are gone (PLAN `FLAGS-01`).** They have no alias. Each one
  stops with exit 2 and a message that says why.

  | Removed | Why | What to do |
  |---|---|---|
  | `--f1`, `-f1` | Nothing read it. `--reference` writes its own `cdna.fa` under `<output>/index`. | Drop the flag. |
  | `--multimap-primary-call` | It accepted one value, `selected-method`, the v3 count contract. | Drop the flag. |
  | `--se-threshold` | It never changed a call or a count. It only set the cut-off line of the exploratory high-burden plot. | Drop the flag. The cut-off is fixed at 10 molecules, the old default. |
  | `--positive-control-scope virus_key` | The value needs an approved transfer calibration (`VAL-RA-CAL`) that this release lacks. `RunConfig` already refused it. | Use `exact_sequence` with `--positive-control-virus-key`, or `panel_mechanics`. |

  `run_fingerprint` is unchanged for a run that used only kept flags at the same values: the two
  fixed options still enter the hash at their old values. A run that used `--f1` hashed that FASTA, so
  its manifest is no longer the same run; `--resume` already refuses manifests from another version.
- **Renamed options keep the old spelling as a deprecated alias (PLAN `FLAGS-01`).** The old name
  prints a `viralscan: warning:` line to stderr and is removed in 3.1.0. Defaults, values and outputs
  are unchanged. Output columns and config keys that carry `umi` (`marker_min_umi`, `programme_min_umi`,
  `hostresponse_min_umi`) keep their names.

  | Verb | Old | New | Why |
  |---|---|---|---|
  | `run`, `rerun programs` | `--programme-min-umi` | `--programme-min-molecules` | It sets a molecule count, as `--programme-productive-min-molecules` does. |
  | `run` | `--hostresponse-min-umi` | `--hostresponse-min-molecules` | It is a threshold on the viral molecule estimate. |
  | `rerun hostresponse` | `--min-umi` | `--min-molecules` | Same threshold, same name as the `run` option. |
  | `discover run` | `--host-subtract` | `--host-filter` | `run` calls the same STAR step `--host-filter`; `--host-filter-star-params` already did. |
  | `discover run` | `--called-cells` | `--called-cells-file` | `run` names its called-cell list `--called-cells-file`. |
  | `translated` | `--cells` | `--called-cells-file` | Same list, one name. |

- **`build ref --virus-accessions` takes commas as well as spaces** (`A,B` equals `A B`), like the other list options.
- **`run --help` lists 16 options.** `viralscan run --help-all` lists all of them. Every option parses
  either way.
- **Commands are grouped (PLAN `CMDGROUP-01`).** The post-hoc, check and build verbs moved under
  five groups. The old names still run the same handler and print a `viralscan: warning:` line that
  names the new one. They are removed in 3.1.0. Options, outputs and exit codes are unchanged.

  | Old | New |
  |-----|-----|
  | `check-chemistry` | `check chemistry` |
  | `check-whitelist` | `check whitelist` |
  | `validate-run` | `check run` |
  | `doctor` | `check env` |
  | `discover` | `discover run` |
  | `discover-collapse` | `discover collapse` |
  | `discover-states` | `discover states` |
  | `rerun-multimap` | `rerun multimap` |
  | `rerun-programs` | `rerun programs` |
  | `hostresponse` | `rerun hostresponse` |
  | `build-ref` | `build ref` |
  | `star-index` | `build star-index` |
  | `translated-index` | `build translated-index` |
  | `host-loss` | `dev host-loss` (hidden group, not a supported interface) |

  `run`, `evidence`, `evidence-plot`, `translated`, `status` and `data fetch` stay as they are.
  The post-hoc record (`posthoc/records/`) keeps the old verb name in its `verb` field, whichever
  spelling you typed, so existing records and `viralscan status` read the same.
- **Output files.** `adata_multimap.h5ad` always carries `uns["multimap_diagnostics"]`; the layout
  (version 3) and the schema (3.0.0) are unchanged.

## Command-line changes

Quantification is now the `run` subcommand. A bare `viralscan [OPTIONS]` still runs it, with a
deprecation warning. The spellings below were removed without aliases; each one now fails with
the replacement named. Configuration keys, `config.yaml` and `run_fingerprint` are unchanged.

| Subcommand | Old | New |
|---|---|---|
| (quantification) | `viralscan -s1 R1 -s2 R2 ...` | `viralscan run --r1 R1 --r2 R2 ...` |
| `run`, `check-whitelist`, `check-chemistry` | `-s1` / `--sample1`, `-s2` / `--sample2` | `--r1`, `--r2` |
| `run` | `-ref`, `-reference`, `-gtf`, `-fasta`, `-umap`, `-acc`, `-mm` | `--reference`, `--gtf`, `--fasta`, `--umap`, `--ncbi-accession`, `--multimapping` |
| `run` | `-f1`, `--f1` | removed in 3.0.0b1 (it had no effect); see "Migrating from 3.0.0rc1" |
| `run` | `-v` | `--visual` / `--no-visual` |
| `run` | `--reference True`, `--visual False`, `--multimapping False` (value form) | `--reference`, `--no-visual`, `--no-multimapping` |
| `run` | `--data-cache-dir` | `--cache-dir` |
| `run` | `--cell-types CSV` | `--cell-type-labels CSV` |
| `run` | `--hostresponse-control-mito` | `--hostresponse-mito-control` |
| `hostresponse` | `--hostresponse-min-umi`, `--detection-threshold` | `--min-molecules` |
| `hostresponse` | `-o RUN_DIR` (an input run directory) | `--run-dir RUN_DIR`; `--output DIR` now optionally names where hostresponse writes (default `<run-dir>/posthoc/hostresponse`) |
| `evidence`, `rerun-programs`, `hostresponse`, `host-loss` | `--run-dir SAMPLE_DIR` | `--run-dir RUN_ROOT [--sample NAME]`; the sample-directory form is deprecated (warning, removed in 3.1.0) |
| `status` | `--run-dir` a sample directory | `--run-dir RUN_ROOT`; the sample-directory form is deprecated (warning, removed in 3.1.0) |
| `hostresponse` | `--no-use-hvg` (still valid) | also `--use-hvg` |
| `validate-run` | positional `RUN_DIR` | `--run-dir RUN_DIR` |
| `validate-run` | `--json-output PATH` | `--json PATH` |
| `doctor` | `--json` (switch; the report was always printed) | `--json PATH` |

### Layout v3 changes (breaking, before 3.0.0 stable)

`run_manifest.json` records `layout_version: 3`. Version 1 and 2 runs still validate, and every
reader and post-hoc verb falls back to their paths. `rerun-multimap` and `run --resume` on an older run
copy `log/analysis.txt` to the new path first, so nothing is recomputed.

| What | Before (v1/v2) | Now (v3) |
|---|---|---|
| Index viral gene list | `<sample>/log/analysis.txt` | `<sample>/results/analysis_genes.txt` |
| Detected genes with counts | `<sample>/log/found_genes.txt` | `<sample>/results/found_genes.txt` |
| Super-expressor tables | `<sample>/plots/SuperExpressor_<virus>.{tsv,json}` | `<sample>/results/superexpressor/` (the PNG stays in `plots/`) |
| In-run host-response output | `<sample>/hostresponse/` | `<sample>/results/hostresponse/` |
| `discover` / `run --discover` default | `<run>/posthoc/discovery/<sample>/` | `<run>/<sample>/posthoc/discovery/` |
| `translated` default | `<run>/posthoc/translated/<sample>/` | `<run>/<sample>/posthoc/translated/` |

Post-hoc verbs now have one level, `<sample>/posthoc/<verb>/`. `rerun-multimap` is the exception: it
copies the whole run, so it stays at `<run>/posthoc/multimap_<method>/`. `log/` keeps logs and `.done`
sentinels only. New files: `<sample>/README.txt` (what is where, and which file answers which
question), a run-level `README.txt`, and `<sample>/results/artifact_inventory.json` (path, size, sha256
per file). `kb-python/` and the other work dirs stay at the sample root; the README marks them as
intermediates. Scripts that read the old paths must read the new ones, or try the new path first.

### Alias layers dropped from `adata_multimap.h5ad`

Four layers duplicated another matrix and are no longer written. Existing h5ads keep them. Scripts that read a new run must use the new name; the package itself reads the new name and falls back to the old one.

| Old layer | New name |
|---|---|
| `counts_original` | `counts_unique` |
| `counts_corrected` | `counts_ambiguous_allocated` |
| `counts_combined` | `adata.X` |
| `counts_host_viral_selected_sibling_weighted` | `counts_host_viral_selected_equal` |

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
