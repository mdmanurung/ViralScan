# CLI Reference

Run `viralscan --help` or `viralscan build-ref --help` to see all options
for the current installed version. This page documents every flag as of
ViralScan **2.2.0**.

---

## `viralscan` — viral quantification (default mode)

```
viralscan [OPTIONS]
```

### Input / output

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--output PATH` | `-o` | *(required)* | Output directory |
| `--sample1 PATH` | `-s1` | *(required)* | Forward FASTQ (comma-separated for multiple) |
| `--sample2 PATH` | `-s2` | *(required)* | Reverse FASTQ (comma-separated for multiple) |

### Reference (choose one of three modes)

| Flag | Short | Description |
|------|-------|-------------|
| `--index PATH` | `-i` | Pre-built kallisto index |
| `--transcripts PATH` | `-t` | t2g.txt transcript-to-gene map |
| `--f1 PATH` | `-f1` | cDNA FASTA (lamanno/nucleus/kite mode) |
| `--reference` | `-ref` | Build index from `-fasta` + `-gtf` |
| `--fasta PATH` | `-fasta` | FASTA file(s), comma-separated |
| `--gtf PATH` | `-gtf` | GTF file(s), comma-separated |
| `--ncbi-accession ACC` | `-acc` | NCBI accession(s), comma-separated; fetch + build |
| `--ncbi-email EMAIL` | | Contact email for NCBI (or `$NCBI_EMAIL`) |

### Analysis parameters

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--technology STRING` | `-x` | `10xv3` | Single-cell technology (`kb --list` for all) |
| `--whitelist PATH` | `-w` | *(bundled)* | Barcode whitelist file |
| `--cores N` | `-c` | `6` | CPU cores |
| `--multimapping` / `--no-multimapping` | `-mm` | on | Multimapping correction |
| `--umap` | `-umap` | off | Generate UMAP plot |
| `--visual` / `--no-visual` | `-v` | on | Generate visualisations |

### Detection thresholds

| Flag | Default | Description |
|------|---------|-------------|
| `--detection-threshold N` | `1` | Min total viral UMI to call a virus detected |
| `--se-threshold N` | `10` | UMI count to call a cell a "super-expressor" |
| `--min-counts N` | `1000` | Min total UMI per cell (UMAP QC) |
| `--min-genes N` | `200` | Min detected genes per cell (UMAP QC) |
| `--cell-types PATH` | *(none)* | CSV (barcode,cell_type) for per-type enrichment |

### Verbosity

| Flag | Description |
|------|-------------|
| `--verbose` | Enable DEBUG-level logging |
| `--quiet` | Suppress INFO; show warnings and errors only |

---

## `viralscan build-ref` — reference builder

```
viralscan build-ref [OPTIONS]
```

Build a combined host + virus kallisto reference without running a full
analysis.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--host SPECIES` | | *(none)* | Host species (see `--list-species`) |
| `--virus-accessions ACC [ACC ...]` | | *(none)* | NCBI accession(s) |
| `--output PATH` | `-o` | `viralscan_ref` | Output directory |
| `--ncbi-email EMAIL` | | *(none)* | Contact email for NCBI |
| `--ncbi-api-key KEY` | | *(none)* | NCBI API key |
| `--cache-dir PATH` | | `~/.cache/viralscan/` | Download cache root |
| `--no-kb-ref` | | off | Stop after writing FASTA + GTF; skip `kb ref` |
| `--list-species` | | off | Print supported host species and exit |
| `--verbose` | | off | DEBUG-level logging |
| `--quiet` | | off | Warnings + errors only |

### Supported host species

Run `viralscan build-ref --list-species` for the current list.
Examples: `human`, `mouse`, `rat`, `zebrafish`, `chicken`, `macaque`, `pig`.
