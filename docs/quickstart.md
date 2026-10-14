# Quickstart

This guide shows the shortest path from paired FASTQ files to ViralScan
results. All examples assume you have activated the `viralscan` conda
environment and that your files are in the current directory.

---

## Before you run

Check that the command-line tools are available:

```bash
viralscan --help
kb --version
snakemake --version
```

Choose one reference mode:

| Use this mode | When you have |
|---------------|---------------|
| Pre-built index | `index.idx` and `t2g.txt` from `kb ref` |
| FASTA + GTF | Viral or host+virus FASTA/GTF files and want ViralScan to run `kb ref` |
| NCBI accession | RefSeq/GenBank accessions and want ViralScan to download the viral reference |

ViralScan writes one output folder per sample under `--output`. The sample
folder name is inferred from the R1 filename before the first underscore. For
example, `sample_R1.fastq.gz` writes results under `output/sample/`.

## Mode 1 — Pre-built index

Use this when you already have a kallisto index (`index.idx`) and
transcript-to-gene map (`t2g.txt`) created by `kb ref`.

```bash
viralscan \
  -t t2g.txt \
  -i index.idx \
  -o output/ \
  -s1 sample_R1.fastq.gz \
  -s2 sample_R2.fastq.gz
```

Key results are written to:

```text
output/sample/results/viral_summary.tsv
output/sample/results/per_cell_viral.tsv
output/sample/report.html
output/sample/kb-python/counts_unfiltered/adata_multimap.h5ad
```

Confirm the run is complete and its artifacts match their schemas:

```bash
viralscan validate-run output/
```

A good run prints a JSON report with `"ok": true` and an empty `"issues"` list.
The run-level `output/run_manifest.json` and `output/run_complete.json` sit beside
the per-sample folders.

The default `--cell-calling auto` runs emptyDrops on the host transcriptome only, which needs a realistic number
of barcodes. On a toy input with only a couple of barcodes it stops with a
`CellCallingError`; add `--cell-calling none` to report over all barcodes
deliberately (every ViralScan count then includes empty droplets).

---

## Mode 2 — Build index from FASTA + GTF

Provide FASTA and GTF files and ViralScan will run `kb ref` for you.

```bash
viralscan \
  --reference \
  -fasta viral_genome.fasta \
  -gtf   viral_annotation.gtf \
  -o     output/ \
  -s1    sample_R1.fastq.gz \
  -s2    sample_R2.fastq.gz
```

Multiple FASTA or GTF files can be supplied as comma-separated paths
(no spaces):

```bash
-fasta host.fasta,virus.fasta -gtf host.gtf,virus.gtf
```

The index is built under `output/index/` (`index.idx`, `t2g.txt`,
`index.idx.build_manifest.json`, `reference_manifest.json`); pass those two files
to Mode 1 to reuse the reference without rebuilding it.

---

## Mode 3 — Download from NCBI by accession

ViralScan fetches the FASTA and GTF automatically from NCBI, builds the index,
and then runs quantification. NCBI requires a contact email whenever a sequence has
to be downloaded; pass `--ncbi-email` or set `NCBI_EMAIL`. A run whose accessions are
all in the local cache (`~/.cache/viralscan/ncbi/`) needs no email.

```bash
viralscan \
  -acc   NC_002021.3 \
  -o     output/ \
  -s1    sample_R1.fastq.gz \
  -s2    sample_R2.fastq.gz \
  --ncbi-email you@example.org
```

Multiple accessions are comma-separated:

```bash
-acc NC_002021.3,NC_001512.1
```

Downloaded references are cached in `~/.cache/viralscan/ncbi/` and reused
on subsequent runs.

---

## Building a combined host + virus reference

Use `viralscan build-ref` to create a combined kallisto index from an
Ensembl host transcriptome and NCBI viral genomes. This is the recommended
host-aware workflow for routine runs:

```bash
viralscan build-ref \
  --host human \
  --virus-accessions NC_045512.2 NC_002021.3 \
  --output ref_human_sars/ \
  --ncbi-email you@example.org
```

Then use the resulting files with Mode 1:

```bash
viralscan \
  -t ref_human_sars/t2g.txt \
  -i ref_human_sars/index.idx \
  -o output/ \
  -s1 sample_R1.fastq.gz \
  -s2 sample_R2.fastq.gz
```

The default `--multimap-method conservative` (recorded as `host-conservative`) keeps host-virus ambiguous
equivalence-class mass out of primary viral estimates — the conservative choice for
combined host+virus references, where host-virus cross-homology can inflate
candidate viral evidence. Use `--multimap-method equal` for an explicit
equal-allocation comparison.

Add `--multimap-molecule-assignments` to also write
`results/molecule_assignments.tsv.gz`: one row per corrected cell-barcode/UMI molecule with
its gene set, status and per-gene weights under the selected method. It is off by default,
changes no count, and the file is large on deep samples; see the
{ref}`output reference <molecule-assignments>`.

Run `viralscan build-ref --list-species` to see supported host names.
Expanded anellovirus sequences are never included implicitly; add
`--anellovirus` only for an explicit screening reference. Reference builds
fail closed on missing accessions, duplicate identifiers/sequences, or a
requested tool that is unavailable, and record provenance in
`reference_manifest.json`.

For libraries at risk of intronic/intergenic host homology, build with the full
host genome as a D-list:

```bash
viralscan build-ref --host human --virus-accessions NC_002021.3 \
  --genome-dlist GRCh38.primary_assembly.genome.fa -o ref_human_ebv/
```

This both masks shared host-genomic k-mers in the kallisto index and writes the
raw per-virus homology measurements. It is compute- and memory-intensive.

---

## Optional cell-type enrichment

If you have cell annotations from Seurat, scanpy, or another workflow, write a
CSV with two columns:

```text
barcode,cell_type
AAACCCAAGAGT-1,T cell
AAACCCAGTGCA-1,Monocyte
```

Then add it to any quantification command:

```bash
viralscan \
  -t t2g.txt \
  -i index.idx \
  -o output/ \
  -s1 sample_R1.fastq.gz \
  -s2 sample_R2.fastq.gz \
  --cell-types cell_types.csv
```

The enrichment table is written to
`output/sample/results/cell_type_enrichment.tsv`.

---

## Useful flags

| Flag | Description |
|------|-------------|
| `-c N` | Use N cores (default: 6) |
| `--umap` | Generate UMAP plot (increases runtime) |
| `--no-multimapping` | Skip multimapping correction |
| `--multimap-method METHOD` | Host/viral ambiguity handling; default is `conservative` (`equal` is an explicit equal-allocation comparison, `em` iterated cell-local allocation) |
| `--multimap-molecule-assignments` | Also write `results/molecule_assignments.tsv.gz`, one row per corrected molecule (off by default; changes no count; large on deep samples) |
| `--detection-threshold N` | Min estimated viral molecules in called cells for a virus to pass the detection gate (default: 1) |
| `--detection-min-cells N` | Distinct called cells that must be infected to pass the gate (default: 1); other viruses stay listed as candidates |
| `--infected-min-molecules X` | Summed allocated molecules for a cell to count as infected (default: 1.0); weaker cells are counted as sub-threshold |
| `--cell-types CSV` | Add per-virus cell-type enrichment to the report |
| `--host-filter starsolo --host-index PATH` | Optional genome-level host pre-subtraction before quantification |
| `--verbose` | Enable debug logging |

See the [CLI reference](cli_reference.md) for all options.
