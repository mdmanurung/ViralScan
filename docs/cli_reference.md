# CLI Reference

Run `viralscan --help`, `viralscan data fetch --help`, or
`viralscan build-ref --help` to see options for the installed version. This
page documents the public CLI as of ViralScan **2.2.0**.

---

## `viralscan` — viral quantification (default mode)

```
viralscan [OPTIONS]
```

Use this command with paired FASTQ files. ViralScan validates inputs, prepares
or reuses the reference, then dispatches the Snakemake workflow.

### Input / output

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--output PATH` | `-o` | *(required)* | Root output directory. Each FASTQ pair is written to a sample subdirectory. |
| `--sample1 PATHS` | `-s1` | *(required)* | R1 FASTQ path, or comma-separated R1 paths for multiple samples. |
| `--sample2 PATHS` | `-s2` | *(required)* | R2 FASTQ path, or comma-separated R2 paths matching `--sample1`. |

### Reference (choose one of three modes)

| Flag | Short | Description |
|------|-------|-------------|
| `--index PATH` | `-i` | Pre-built kallisto index from `kb ref`. Required in pre-built mode. |
| `--transcripts PATH` | `-t` | `t2g.txt` transcript-to-gene map from `kb ref`. Required in pre-built mode. |
| `--f1 PATH` | `-f1` | Optional cDNA FASTA passed through for kb workflows that need it. |
| `--reference` | `-ref` | Build index from `-fasta` + `-gtf` |
| `--fasta PATH` | `-fasta` | FASTA file(s), comma-separated |
| `--gtf PATH` | `-gtf` | GTF file(s), comma-separated |
| `--ncbi-accession ACC` | `-acc` | NCBI accession(s), comma-separated; fetch + build |
| `--ncbi-email EMAIL` | | Contact email for NCBI (or `$NCBI_EMAIL`) |

Reference modes are mutually exclusive:

- Pre-built: use `-i` and `-t`.
- FASTA/GTF: use `--reference -fasta ... -gtf ...`.
- NCBI: use `-acc ...` and provide an NCBI email.

### Analysis parameters

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--technology STRING` | `-x` | `10xv3` | Single-cell technology (`kb --list` for all) |
| `--whitelist PATH` | `-w` | *(bundled)* | Barcode whitelist file |
| `--cores N` | `-c` | `6` | CPU cores |
| `--multimapping` / `--no-multimapping` | `-mm` | on | Multimapping correction |
| `--multimap-method METHOD` | | `equal` | Multimapper allocation: `equal`, `host-conservative`, or `unique-weighted` |
| `--multimap-pseudocount FLOAT` | | `1.0` | Positive pseudocount for `unique-weighted` |
| `--multimap-primary-call MODE` | | `legacy` | Viral calling policy: `legacy`, `unique-only`, or `confidence` |
| `--umap` | `-umap` | off | Generate UMAP plot |
| `--visual` / `--no-visual` | `-v` | on | Generate visualisations |
| `--host-filter ALIGNER` | | *(none)* | Optional host subtraction before quantification. Choices: `starsolo`, `kallisto`. |
| `--host-index PATH` | | *(none)* | Required with `--host-filter`. STAR genome directory for `starsolo`; kallisto cDNA index for `kallisto`. |

### Detection thresholds

| Flag | Default | Description |
|------|---------|-------------|
| `--detection-threshold N` | `1` | Min total viral UMI to call a virus detected |
| `--se-threshold N` | `10` | UMI count to call a cell a "super-expressor" |
| `--cell-types PATH` | *(none)* | CSV with `barcode,cell_type` columns for per-virus cell-type enrichment |

### UMAP / QC parameters

These flags affect only the optional `--umap` workflow.

| Flag | Default | Description |
|------|---------|-------------|
| `--min-counts N` | `1000` | Min total UMI per cell (UMAP QC) |
| `--min-genes N` | `200` | Min detected genes per cell (UMAP QC) |
| `--hvg-min-mean X` | `0.0125` | Scanpy highly-variable-gene `min_mean` |
| `--hvg-max-mean X` | `3.0` | Scanpy highly-variable-gene `max_mean` |
| `--hvg-min-disp X` | `0.5` | Scanpy highly-variable-gene `min_disp` |
| `--umap-n-neighbors N` | `15` | Neighbor count for Scanpy graph construction |

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
analysis. The command downloads the host cDNA FASTA/GTF from Ensembl, viral
FASTA/GTF from NCBI, concatenates them, and runs `kb ref` unless
`--no-kb-ref` is supplied.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--host SPECIES` | | *(none)* | Host species (see `--list-species`) |
| `--virus-accessions ACC [ACC ...]` | | *(none)* | One or more NCBI accessions separated by spaces |
| `--output PATH` | `-o` | `viralscan_ref` | Output directory |
| `--ncbi-email EMAIL` | | *(none)* | Contact email for NCBI |
| `--ncbi-api-key KEY` | | *(none)* | NCBI API key |
| `--cache-dir PATH` | | `~/.cache/viralscan/` | Download cache root |
| `--no-kb-ref` | | off | Stop after writing FASTA + GTF; skip `kb ref` |
| `--list-species` | | off | Print supported host species and exit |
| `--verbose` | | off | DEBUG-level logging |
| `--quiet` | | off | Warnings + errors only |

Example:

```bash
viralscan build-ref \
  --host human \
  --virus-accessions NC_045512.2 NC_002021.3 \
  --output ref_human_virus/ \
  --ncbi-email you@example.org
```

Expected outputs when `kb ref` succeeds:

| File | Use |
|------|-----|
| `combined.fa` | Concatenated host + virus FASTA |
| `combined.gtf` | Concatenated host + virus GTF |
| `index.idx` | Pass to `viralscan -i` |
| `t2g.txt` | Pass to `viralscan -t` |
| `cdna.fa` | cDNA FASTA produced by `kb ref -f1` |

### Supported host species

Run `viralscan build-ref --list-species` for the current list.
Examples: `human`, `mouse`, `rat`, `zebrafish`, `chicken`, `macaque`, `pig`.

---

## `viralscan data fetch` — bundled viral panel

```
viralscan data fetch [OPTIONS]
```

Download the bundled viral GTF panel from Zenodo, verify checksums, and unpack
the GTF files into the local ViralScan cache.

| Flag | Default | Description |
|------|---------|-------------|
| `--cache-dir PATH` | `~/.cache/viralscan/` | Cache root. GTF files are written to `PATH/data/`. |
| `--url URL` | Zenodo archive URL | Override archive URL, mainly for tests or mirrors. |
| `--sha256 DIGEST` | *(none)* | Optional expected SHA-256 digest for the downloaded archive. |
| `--force` | off | Re-download and replace cached GTF files. |
| `--verbose` | off | DEBUG-level logging. |
| `--quiet` | off | Warnings + errors only. |

Typical use:

```bash
viralscan data fetch
```
