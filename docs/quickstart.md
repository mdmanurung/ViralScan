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

---

## Mode 3 — Download from NCBI by accession

ViralScan fetches the FASTA and GTF automatically from NCBI, builds the index,
and then runs quantification. NCBI requires a contact email; pass
`--ncbi-email` or set `NCBI_EMAIL`.

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

The default `--multimap-method host-conservative` keeps host-virus ambiguous
equivalence-class mass out of primary viral estimates — the conservative choice for
combined host+virus references, where host-virus cross-homology can inflate
candidate viral evidence. Use `--multimap-method equal` for an explicit
equal-allocation comparison.

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

```csv
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
| `--multimap-method METHOD` | Host/viral ambiguity handling; default is `host-conservative` (`equal` is an explicit equal-allocation comparison) |
| `--detection-threshold N` | Min estimated viral molecule support to report a candidate virus (default: 1) |
| `--cell-types CSV` | Add per-virus cell-type enrichment to the report |
| `--host-filter starsolo --host-index PATH` | Optional genome-level host pre-subtraction before quantification |
| `--verbose` | Enable debug logging |

See the [CLI reference](cli_reference.md) for all options.
