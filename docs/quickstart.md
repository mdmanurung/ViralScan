# Run ViralScan

This is the recommended path for one routine host-aware run: build one combined
host + virus reference, quantify, validate the output, and open the results.
Every command below assumes the full installation
([Installation](installation.md)) and paired-end FASTQ files in the current
directory.

## Before you run

Check the toolchain before starting a long job:

```bash
viralscan --help
viralscan doctor --profile full
kb info
snakemake --version
```

`viralscan doctor --profile full` prints a JSON dependency report and exits
nonzero when a required full-workflow tool is missing. A green `pip` profile does
not imply that `kb`, `snakemake` or STAR are installed.

| Requirement | What to check |
|---|---|
| Reads | Paired-end FASTQ, optionally gzipped; both `--r1` (R1) and `--r2` (R2) are required. R1 carries the cell barcode and UMI. |
| Chemistry | 10x 3′ and 5′ are supported, and Drop-seq is supported. ViralScan detects the chemistry from the first 100k R1 reads; an explicit `-x` that the reads contradict stops the run. A GEM-X 5′ library needs its on-list supplied with `-w`. |
| Reference | A kallisto `index.idx` and `t2g.txt`. Step 1 builds a combined host + virus reference; an existing compatible index can be reused instead. |
| Cell calling | The default `--cell-calling auto` runs EmptyDrops, which needs `Rscript` with `DropletUtils`. Supply an external list with `--called-cells-file` when you already have one. |
| Resources | `-c/--cores` defaults to 6. Building a combined reference downloads the Ensembl host cDNA and NCBI viral records, then runs `kb ref`; allow disk for `combined.fa`, `combined.gtf` and the index. |

ViralScan refuses a non-empty output directory. Reuse is explicit: `--resume`
requires an exact fingerprint match, and `--overwrite` replaces the resolved
output target after confirmation.

## Step 1 — Build a combined host + virus reference

This is the recommended reference for routine host-aware analysis. Host and
viral features share one index, so host–virus ambiguity stays visible instead of
being deleted before viral quantification.

```bash
viralscan build-ref \
  --host human \
  --virus-accessions NC_045512.2 NC_002021.3 \
  --output ref_human_virus/ \
  --ncbi-email you@example.org
```

Replace the accessions with the viruses you expect. Run
`viralscan build-ref --list-species` for the supported host names. The build
writes:

- `ref_human_virus/index.idx` — pass to `-i`
- `ref_human_virus/t2g.txt` — pass to `-t`
- `ref_human_virus/index.idx.build_manifest.json` — the index's gene-set receipt, read by every run
- `ref_human_virus/reference_manifest.json` — sources, digests and build receipt
- `ref_human_virus/combined.fa` and `combined.gtf` — the concatenated inputs

Expanded anellovirus sequences are never included implicitly. Add `--anellovirus`
only for an explicit screening reference. Builds fail closed on missing
accessions, duplicate identifiers or sequences, and unavailable tools.

**Already have a compatible index?** Skip this step and pass its `index.idx` and
`t2g.txt` to Step 2. An index from `viralscan build-ref` carries
`<index>.build_manifest.json` (beside `index.idx`) and needs nothing else. For an
index without that manifest, run `viralscan data fetch` once so the viral gene
list comes from the bundled panel; custom `--gtf` runs and `--ncbi-accession` runs do not
need it.

## Step 2 — Quantify

```bash
viralscan run \
  -t ref_human_virus/t2g.txt \
  -i ref_human_virus/index.idx \
  -o output/ \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz
```

The sample directory is inferred from the R1 filename before the first
underscore, so this run writes to `output/sample/`.

Defaults worth knowing:

- `--multimap-method conservative` (recorded as `host-conservative`) keeps host–virus ambiguous molecule mass
  out of primary viral estimates. This is the conservative choice for a combined
  reference. `--multimap-method equal` is an explicit equal-allocation
  comparison, not a better default; `em` is iterated cell-local allocation.
- `--cell-calling auto` uses an external list when `--called-cells-file` is
  supplied and otherwise runs EmptyDrops on the host transcriptome only (never on viral reads). The called-cell set is the denominator
  behind every `*_called` field.
- Multiple samples in one run: pass comma-separated, positionally matched paths
  to `--r1` and `--r2`, for example
  `--r1 A_R1.fastq.gz,B_R1.fastq.gz --r2 A_R2.fastq.gz,B_R2.fastq.gz`.

### The called-cell requirement

Cell calling decides which barcodes the per-cell rates are over. EmptyDrops needs
a realistic number of barcodes; it is the default because an all-barcode
denominator is diluted by empty droplets.

- If you already have called cells, pass `--called-cells-file cells.tsv` with
  `--cell-calling auto` (or `external`). The full external list sets the
  denominator, including cells with no viral molecules.
- `--cell-calling knee` is a sensitivity-only check, not a source of reported
  numbers.
- `--cell-calling none` is for deliberate toy data only. It reports over all
  barcodes, so every count includes empty droplets. It is **not** a routine
  workaround for failed cell calling; fix the chemistry, whitelist or external
  list instead.

## Step 3 — Validate the run

```bash
viralscan validate-run --run-dir output/
```

Validation prints a JSON report. A successful run has `"ok": true` and no
`error`-level entries:

```json
{
  "ok": true,
  "run_dir": "/path/to/output",
  "issues": []
}
```

Every issue has a `level` (`error` or `warning`), a `code`, a message and a
path. The command exits nonzero when `ok` is false. It checks technical
consistency:

- the run manifest is present and matches the v3 schema;
- the run fingerprint and recorded input hashes still match the files on disk;
- `run_complete.json` is present when the manifest requires it, belongs to this
  run, and its recorded artifact hashes match;
- the count H5AD is readable and satisfies the v3 count contract, including the
  molecule-conservation invariants.

A warning is reported when an older manifest declares no completion marker. Use
`--no-verify-inputs` only when the inputs are intentionally unavailable (for
example after moving a run); it skips input re-hashing.

Validation does **not** establish biological validity. It cannot tell you whether
the reference covers your virus, whether a call is infection, or whether a
negative is informative; see [Results](results.md).

## Step 4 — Open the results

For `sample_R1.fastq.gz`, the main artifacts are:

```text
output/sample/report.html
output/sample/results/viral_summary.tsv
output/sample/results/per_cell_viral.tsv
output/sample/results/called_cells.tsv
output/run_complete.json
output/run_manifest.json
```

[Results](results.md) explains which file answers which question, what the
molecule estimates and denominators mean, how to read one row of
`viral_summary.tsv`, and the five checks to run before a candidate becomes a
claim.

## Other reference modes

The workflow above is the recommended one. The other routes (a pre-built index,
`--reference --fasta ... --gtf ...`, `--ncbi-accession`, and the extended
Anelloviridae panel) and the host-control choices beyond the combined reference
(STARsolo pre-subtraction, genome D-list masking) are in the
[Reference panel](reference_panel.md); read
{ref}`Host control and reference strategies <host-control>` before choosing one.
Without `--genome-dlist`, `build-ref` warns that host homology was not measured
(which is not the same as zero), and it warns about reference records that are
near-copies of each other on either strand; `--strict-reference` fails the build
on such a pair instead (see [Panel integrity](panel_integrity.md)).
