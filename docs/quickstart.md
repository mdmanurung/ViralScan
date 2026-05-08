# Quickstart

This guide shows three common ways to run ViralScan. All examples assume
you have activated the `viralscan` conda environment and your files are in
the current directory.

---

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

Results are written to `output/results/`.

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

ViralScan fetches the FASTA and GTF automatically from NCBI RefSeq.

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
Ensembl host transcriptome and NCBI viral genomes:

```bash
viralscan build-ref \
  --host human \
  --virus-accessions NC_045512.2 \
  --output ref_human_sars/ \
  --ncbi-email you@example.org
```

Then use the resulting `ref_human_sars/` directory with Mode 1.

---

## Useful flags

| Flag | Description |
|------|-------------|
| `-c N` | Use N cores (default: 6) |
| `--umap` | Generate UMAP plot (increases runtime) |
| `--no-multimapping` | Skip multimapping correction |
| `--detection-threshold N` | Min viral UMI to call a virus detected (default: 1) |
| `--verbose` | Enable debug logging |

See the [CLI reference](cli_reference.md) for all options.
