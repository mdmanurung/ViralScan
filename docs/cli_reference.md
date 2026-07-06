# CLI Reference

Run `viralscan --help` or `viralscan <subcommand> --help` to see options for
the installed version. This page documents the public CLI as of ViralScan
**2.5.0**. Available subcommands: `build-ref`, `data fetch`, `evidence`,
`rerun-multimap`, `hostresponse`, `check-whitelist`.

---

## `viralscan` — viral quantification (default mode)

```
viralscan [OPTIONS]
```

Use this command with paired FASTQ files. ViralScan validates inputs, prepares
or reuses the reference, then dispatches the Snakemake workflow.

For host-aware viral detection, the recommended workflow is to build a combined
host+virus reference with `viralscan build-ref`, then pass its `index.idx` and
`t2g.txt` to this command. The default multimapping method is
`host-conservative`, which keeps host-virus ambiguous EC mass out of primary
viral counts — the safer choice for combined host+virus references where
host-virus cross-homology can inflate viral calls. Use `--multimap-method equal`
for a fast, unbiased first pass, or `em` for iterated maximum-likelihood allocation.

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
| `--data-cache-dir PATH` | | Viral annotation cache root for the bundled panel (`PATH/data/`) |

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
| `--multimap-method METHOD` | | `host-conservative` | Multimapper allocation: `host-conservative`, `equal`, `unique-weighted`, or `em` |
| `--multimap-pseudocount FLOAT` | | `1.0` | Positive pseudocount for `unique-weighted` |
| `--multimap-primary-call MODE` | | `legacy` | Viral calling policy: `legacy`, `unique-only`, or `confidence` |
| `--umap` | `-umap` | off | Generate UMAP plot |
| `--visual` / `--no-visual` | `-v` | on | Generate visualisations |
| `--host-filter ALIGNER` | | *(none)* | Optional advanced host subtraction before quantification. Choices: `starsolo`, `kallisto`. |
| `--host-index PATH` | | *(none)* | Required with `--host-filter`. STAR genome directory for `starsolo`; kallisto cDNA index for `kallisto`. |

### Detection thresholds

| Flag | Default | Description |
|------|---------|-------------|
| `--detection-threshold N` | `1` | Min total viral UMI to call a virus detected |
| `--se-threshold N` | `10` | UMI count to call a cell a "super-expressor" |
| `--cell-types PATH` | *(none)* | CSV with `barcode,cell_type` columns for per-virus cell-type enrichment |

### Host-response analysis (optional)

Associate viral presence with host gene expression by providing a host h5ad.
When `--host-h5ad` is supplied, ViralScan trains per-virus L2 logistic
regression models (Luebbert et al. 2026) and runs randomized Lasso stability
selection; results are written to `<output>/hostresponse/`. Pathway enrichment
requires `pip install "viralscan[enrichment]"`.

| Flag | Default | Description |
|------|---------|-------------|
| `--host-h5ad PATH` | *(none)* | Host gene-expression h5ad (cells × genes). Triggers host-response analysis when provided. |
| `--hostresponse-n-seeds N` | `6` | Number of random seeds for the multi-seed L2 regression |
| `--hostresponse-n-stab-iter N` | `100` | Iterations for randomized Lasso stability selection |
| `--hostresponse-use-hvg` / `--no-hostresponse-use-hvg` | on | Use highly variable genes as features (disable to use all genes) |
| `--hostresponse-stab-min-prob PROB` | `0.6` | Min selection probability to call a gene stably associated |
| `--hostresponse-top-n-genes N` | `50` | Top N stable genes to pass to pathway enrichment |
| `--enrichment` | off | Run pathway enrichment via gget (requires `viralscan[enrichment]`) |
| `--enrichment-db DB` | `GO_Biological_Process_2023` | gget.enrichr database |

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
`--no-kb-ref` is supplied. This is the preferred host-aware setup for routine
runs.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--host SPECIES` | | *(none)* | Host species (see `--list-species`) |
| `--virus-accessions ACC [ACC ...]` | | *(none)* | One or more NCBI accessions separated by spaces |
| `--output PATH` | `-o` | `viralscan_ref` | Output directory |
| `--ncbi-email EMAIL` | | *(none)* | Contact email for NCBI |
| `--ncbi-api-key KEY` | | *(none)* | NCBI API key |
| `--cache-dir PATH` | | `~/.cache/viralscan/` | Download cache root |
| `--no-kb-ref` | | off | Stop after writing FASTA + GTF; skip `kb ref` |
| `--anellovirus` / `--no-anellovirus` | | on | Include the packaged Anelloviridae accession table in host+virus references |
| `--reference-panel anellovirus` | | *(none)* | Build a predefined Anelloviridae panel (bundled FASTA if cached, else NCBI download) |
| `--no-mask` | | off | Skip dustmasker low-complexity masking for anellovirus references |
| `--cluster` | | off | Cluster anellovirus sequences at 95% identity (cd-hit-est) after masking |
| `--list-species` | | off | Print supported host species and exit |
| `--verbose` | | off | DEBUG-level logging |
| `--quiet` | | off | Warnings + errors only |

By default, `build-ref` includes the packaged Anelloviridae accession table in
the combined host+virus reference. Pass `--no-anellovirus` when you want a
reference containing only the host and the explicit `--virus-accessions`.

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

---

## `viralscan evidence` — trace reads behind viral calls

```
viralscan evidence [OPTIONS]
```

Trace the reads whose (barcode, UMI) were assigned to viral genes in a
completed ViralScan run, extract them, optionally re-align to a viral genome
for IGV, and score evidence quality (genome coverage + BLAST identity). Use
this to confirm a viral hit is real rather than host cross-homology.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--run-dir PATH` | | *(required)* | A completed ViralScan run output directory |
| `--output PATH` | `-o` | *(required)* | Directory for evidence outputs |
| `--viral-fasta PATH` | | *(none)* | Viral genome FASTA to align extracted reads against (enables BAM/coverage/BLAST). Omit for read-extraction only. |
| `--virus STRING` | | *(none)* | Restrict to viral genes whose ID contains this substring (e.g. `EPSTEIN`, `HERP6B`) |
| `--blast` | | off | BLAST a sample of extracted reads against the viral reference (requires `blast+`) |
| `--cores N` | `-c` | `4` | Threads for minimap2/samtools/blast |
| `--verbose` | | off | Enable DEBUG-level logging |
| `--quiet` | | off | Suppress INFO messages |

Example:

```bash
viralscan evidence \
  --run-dir output/sample/ \
  --viral-fasta viruses.fa \
  --virus EBV \
  --blast \
  -o output/sample/evidence/
```

---

## `viralscan rerun-multimap` — swap multimapping method

```
viralscan rerun-multimap [OPTIONS]
```

Re-run the multimapping step with a different algorithm for an existing run,
skipping the expensive pseudoalignment (kb count). For `equal`,
`host-conservative`, and `unique-weighted` methods, all three count layers are
pre-stored in every multimap h5ad, so the swap is instant (no bus-file
reprocessing). For `em`, the bus file is reprocessed (slower, but still skips
kb count). Detection and UMAP are always re-run after the swap.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--output PATH` | `-o` | *(required)* | The base output directory of the original viralscan run |
| `--multimap-method METHOD` | | *(required)* | Multimapping resolution method: `host-conservative`, `equal`, `unique-weighted`, or `em` |
| `--cores N` | `-c` | `6` | Number of cores for snakemake workers |
| `--multimap-em-max-iter N` | | *(preserve)* | (em only) Maximum EM iterations; omit to keep the existing config value |
| `--multimap-em-tol TOL` | | *(preserve)* | (em only) EM convergence tolerance; omit to keep the existing config value |
| `--verbose` | | off | Enable DEBUG-level logging |
| `--quiet` | | off | Suppress INFO messages |

Examples:

```bash
viralscan rerun-multimap -o out/ --multimap-method host-conservative
viralscan rerun-multimap -o out/ --multimap-method em --cores 8
```

**Tip:** the default (`host-conservative`) is the safe choice; swap to `equal`
for a fast unbiased first pass or `em` for iterated allocation with
`rerun-multimap` — no need to repeat pseudoalignment (non-EM swaps are instant).

---

## `viralscan hostresponse` — host-response analysis

```
viralscan hostresponse [OPTIONS]
```

Associate viral presence with host gene expression via logistic regression
(Luebbert et al. 2026 approach) on a completed viralscan run. For each
detected virus, ViralScan trains an L2 logistic regression model predicting
virus-positive vs. virus-negative cells from host gene expression, then runs
randomized Lasso stability selection to identify robustly associated host genes.
Optional gget pathway enrichment (`--enrichment`) identifies enriched biological
processes in the stable gene set.

The viralscan output directory must already contain `config.yaml` and a
completed multimap step (`log/multimap.done`). Results are written to
`<output>/hostresponse/`.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--output PATH` | `-o` | *(required)* | Existing viralscan sample output directory (contains `config.yaml`) |
| `--host-h5ad PATH` | | *(required)* | Host gene-expression h5ad (cells × genes, matched barcodes) |
| `--n-seeds N` | | from config or `6` | Random seeds for multi-seed L2 regression |
| `--n-stab-iter N` | | from config or `100` | Stability-selection iterations |
| `--no-use-hvg` | | off | Use all genes instead of highly variable genes as features |
| `--stab-min-prob P` | | from config or `0.6` | Min selection probability to call a gene stably associated |
| `--top-n-genes N` | | from config or `50` | Top N stable genes to pass to pathway enrichment |
| `--detection-threshold N` | | from config or `1` | Min UMI count to call a cell virus-positive (raw label) |
| `--label {raw,cpm,fraction}` | | `raw` | Positive-call label (see **Depth-confound controls** below) |
| `--depth-match` | | off | Restrict each virus to a depth-matched cohort so any signal is depth-independent by design |
| `--mito-control` / `--no-mito-control` | | on | Add per-cell %mito as a covariate to the per-gene E-values |
| `--gene-symbols` | | off | Annotate output CSVs with HGNC symbols from Ensembl IDs via mygene.info (network) |
| `--differential` | | off | Write a genome-wide depth/%mito-adjusted differential table (`<virus>_differential.csv`) |
| `--enrichment` | | off | Run pathway enrichment via gget (requires `pip install "viralscan[enrichment]"`) |
| `--enrichment-db DB` | | `GO_Biological_Process_2023` | gget.enrichr database |
| `--verbose` | | off | Enable DEBUG-level logging |
| `--quiet` | | off | Suppress INFO messages |

#### Depth-confound controls (important)

The default positive-call label (`--label raw`, `counts >= detection-threshold`)
**tracks sequencing depth**: deeper cells carry more viral *and* more host counts,
so a naive host-gene AUC is partly a library-size artifact rather than biology.
On the EBV showcase run, sequencing depth *alone* predicts virus status at
AUC 0.97 — higher than the 0.87 host-gene model — so the raw headline is
depth-confounded. To make this visible and correctable, every run reports a
depth-confound baseline, and two opt-in controls give a depth-independent estimate:

- **Always reported** (no flag needed): `hostresponse_metrics.csv` includes
  `depth_alone_auc_mean` (the AUC from sequencing depth alone, under the identical
  split) next to the model AUC, and `<virus>_depth_diagnostics.csv` gives a
  depth-adjusted odds ratio + [E-value](https://doi.org/10.7326/M16-2607) per stable
  gene. If `depth_alone_auc` ≈ the model AUC, treat the headline as depth-confounded.
- **`--label cpm`** (or `fraction`): a depth-normalized, prevalence-matched label —
  positives are the cells with the highest viral UMI *per host UMI*, keeping the same
  number of positives so the AUC stays comparable. Decouples the label from depth.
- **`--depth-match`**: builds a coarsened-exact depth-matched case/control cohort so
  depth alone can no longer discriminate; any surviving host-gene signal is
  depth-independent by construction.

On the EBV showcase data these controls move the model AUC from a confounded
**0.87** (depth-alone 0.97) to an honest **~0.67** (depth-alone ~0.5).

Examples:

```bash
# Minimal — run with defaults from config (raw label; still reports the depth baseline)
viralscan hostresponse -o output/sample/ --host-h5ad host_genes.h5ad

# Depth-independent estimate: CPM label + %mito control + gene symbols
viralscan hostresponse \
  -o output/sample/ \
  --host-h5ad host_genes.h5ad \
  --label cpm --gene-symbols

# Strongest de-confounding: depth-matched cohort + genome-wide differential table
viralscan hostresponse \
  -o output/sample/ \
  --host-h5ad host_genes.h5ad \
  --depth-match --differential

# With stability tuning and pathway enrichment
viralscan hostresponse \
  -o output/sample/ \
  --host-h5ad host_genes.h5ad \
  --n-stab-iter 200 \
  --enrichment
```

> **Pathway enrichment**: requires the optional `gget` dependency.
> Install it with `pip install "viralscan[enrichment]"` before using `--enrichment`.

---

## `viralscan check-whitelist` — chemistry/whitelist preflight

```
viralscan check-whitelist -s1 R1.fastq.gz -w WHITELIST -x TECHNOLOGY
```

Sample the first reads of R1, extract the cell barcode using the technology's
barcode geometry, and report the fraction that match the whitelist. A low match
rate means `--technology` / `--whitelist` do **not** match the library chemistry —
in which case `bustools` silently discards most reads and `kb count` produces an
all-empty-droplet matrix with no error. (In one real case a GEM-X 5′ library
mislabeled as `10xv3` lost 96.5% of its reads before this was caught.)

Run this before a long quantification job when you are unsure of the chemistry.
The main `viralscan` run also performs this check automatically as a **warning**
whenever an explicit `--whitelist` is supplied.

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--sample1 R1` | `-s1` | *(required)* | R1 FASTQ (the barcode read; optionally `.gz`) |
| `--whitelist PATH` | `-w` | *(required)* | Barcode whitelist, one barcode per line (optionally `.gz`) |
| `--technology X` | `-x` | `10xv3` | Single-cell technology (sets the barcode length/offset) |
| `--min-match-rate F` | | `0.5` | Minimum acceptable barcode match rate |
| `--n-sample N` | | `100000` | Number of R1 reads to sample |

Exits `0` when the match rate is acceptable and `1` on a likely mismatch, so it
can gate a pipeline. Example:

```bash
viralscan check-whitelist \
  -s1 sample_R1.fastq.gz \
  -w 3M-february-2018.txt \
  -x 10xv3
```
