# Frequently Asked Questions

---

## Installation

### Why do I get an error about `connection_pool` when installing with pip?

`snakemake` has a transitive dependency (`connection_pool`) that fails to build
with recent versions of setuptools.  Install `snakemake` and `kb-python` via
conda first, then `pip install ViralScan`:

```bash
conda install -c conda-forge -c bioconda snakemake kb-python
pip install ViralScan
```

Or use the provided `environment.yml`:

```bash
conda env create -f environment.yml
```

---

## Running ViralScan

### My run exits with "no reads pseudoaligned" — what does that mean?

`kb count` reported that none of the sequencing reads aligned to the reference.
Common causes:

- The reference was built for a different organism than the sample.
- The technology flag (`-x`) is wrong (e.g. `10xv2` vs `10xv3`).
- The FASTQ files are corrupted or empty.
- The sample files are swapped (R1 ↔ R2).

Useful checks:

```bash
kb --list | less
zcat sample_R1.fastq.gz | head
zcat sample_R2.fastq.gz | head
```

Confirm that `--technology` matches the library chemistry and that R1 is the
barcode/UMI read for your 10x data.

A fast, direct check is the built-in whitelist preflight, which reports the
fraction of R1 barcodes that match the whitelist for a given chemistry:

```bash
viralscan check-whitelist -s1 sample_R1.fastq.gz -w whitelist.txt -x 10xv3
```

A low match rate confirms a chemistry/whitelist mismatch — the single most common
cause of a silent all-empty matrix (e.g. a GEM-X 5′ library mislabeled `10xv3`).
Try other `-x` values until the match rate is high. The main `viralscan` run also
emits this as a warning automatically when an explicit `--whitelist` is given.

### Can I process multiple samples in one run?

Yes — provide comma-separated paths to `-s1` and `-s2`:

```bash
viralscan -t t2g.txt -i index.idx -o output/ \
  -s1 A_R1.fastq.gz,B_R1.fastq.gz \
  -s2 A_R2.fastq.gz,B_R2.fastq.gz
```

This creates separate run directories under `output/`, for example
`output/A/` and `output/B/`.

### The UMAP step takes a very long time. Can I skip it?

Yes — omit `--umap` (it is off by default).  The detection and reporting
steps run without it.

---

## Results

### What does `viral_neighbor_enrichment` measure?

It is a permutation test that asks: are cells with candidate viral molecule support over-represented
among each other's nearest neighbours in the UMAP embedding? A low p-value
indicates spatial clustering of the candidate-support labels beyond what would be expected
after shuffling the viral labels.

### How is `pct_infected` calculated? Why are there two infection percentages?

```
pct_infected        = candidate-support cells / ALL barcodes        × 100
pct_infected_called = candidate-support cells / CALLED cells       × 100
```

The field names are retained for schema compatibility; they do not establish
biological infection. The all-barcode value is diluted by empty droplets. In v3,
`cell_calling=auto` uses an external
CellRanger/STARsolo list when supplied and otherwise runs EmptyDrops. Knee
calling is available only when explicitly requested; there is no silent
fallback. A nonzero viral molecule is candidate evidence, not by itself proof
of infection.

`--detection-threshold` (default 1) is a sample-level threshold for reporting a
candidate virus. It does not change which cells have nonzero molecule support.

### My `hostresponse` model AUC is high — is the host-response signal real?

Check `depth_alone_auc_mean` in `hostresponse_metrics.csv` first. The default
`counts >= threshold` candidate-support label **tracks sequencing depth** (deeper cells
carry more viral *and* more host counts), so a high model AUC can be a library-size
artifact rather than biology. If `depth_alone_auc` is close to the model AUC, treat
the headline as depth-confounded.

To reduce and diagnose measured depth imbalance:

- `--label cpm` — a prevalence-matched label based on the viral molecule estimate
  divided by the total molecule estimate;
- `--depth-match` — a depth-matched case/control cohort;
- per-gene **E-values** in `<virus>_depth_diagnostics.csv` (≥ 2 = robust to moderate
  confounding), reported automatically. `%mito` is controlled by default.

None of these controls proves that all depth, technical, or biological confounding
has been removed. Inspect cohort balance and the depth-only baseline.

### What units is `viral_molecules_per_10k_est` in?

It is a normalized selected-method molecule estimate, not a raw UMI count:

```
viral_molecules_per_10k_est =
    viral_molecules_total_est / molecules_total_est × 10 000
```

---

## Reference panel

### How do I add a virus that is not in the bundled panel?

Use `--ncbi-accession` to fetch any RefSeq nucleotide record:

```bash
viralscan -acc NC_045512.2 -o output/ -s1 R1.fastq.gz -s2 R2.fastq.gz
```

Or use `--reference -fasta my_virus.fasta -gtf my_virus.gtf` with your own
reference files.

### The virus name in the output shows the gene_id prefix (e.g. "HUM_SARS").
### How do I get the full name?

The `VIRUS_NAME_MAP` in `src/viralscan/constants.py` maps prefixes to full
names.  If your virus prefix is not listed, open a GitHub issue or submit a
pull request to add it.

---

## Reducing false positives

### Why might ViralScan report a virus that isn't really there?

ViralScan uses kallisto pseudo-alignment, which is k-mer based.  Any 31-mer
shared between a host transcript and a viral genome will cause host reads to
land on the viral feature.  The main culprits are:

- **Endogenous viral elements (EVEs/HERVs)** — ~8 % of the human genome is
  ancient integrated retroviral sequence, which shares k-mers with many
  exogenous viruses (especially retroviruses).
- **Low-complexity / poly-A regions** — repetitive sequences produce shared
  k-mers that confuse pseudo-alignment.
- **Viral homologs of host genes** — some viruses encode genes with strong
  host homology (e.g. viral IL-10).

The recommended mitigation is a combined host+virus reference built with
`viralscan build-ref` (competitive mapping). By default, multimapping uses the
`host-conservative` method: if a multi-gene equivalence class is compatible
with both host and viral genes, its molecule mass is allocated only among the
compatible host genes. This prevents that mixed ambiguity from creating viral
molecule support while preserving diagnostic evidence in
`results/multimap_evidence.tsv`; formal false-positive performance remains a
truth-panel question.

The `equal` method remains available with `--multimap-method equal`. V3 always
uses the complete selected-method matrix for summaries and reports unique,
virus–virus ambiguous, and host–virus ambiguous molecule evidence separately.

---

### How do I use host pre-subtraction to reduce false positives?

Host pre-subtraction is an optional advanced filter. V3 maps reads to the full host
genome **before** viral quantification and discards fragments that
align. The remaining reads are then passed to `kb count`.

For routine host-aware analysis, prefer the combined host+virus reference plus
the default `--multimap-method host-conservative`. Use pre-subtraction when you
want an extra conservative filter or need to remove host-aligned reads before
viral quantification.

V3 supports one exact-fragment host filter via `--host-filter`:

| Aligner | Flag value | What it does |
|---|---|---|
| STARsolo | `starsolo` | Full genome alignment; unmapped reads collected from STAR's `--outReadsUnmapped Fastx` output |

**STARsolo host subtraction**

Requires a STAR genome directory.  If you do not already have one, build it once:

```bash
# Build STAR genome index (human GRCh38 example; needs ~30 GB RAM, ~30 min)
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_hg38/ \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v44.primary_assembly.annotation.gtf \
     --runThreadN 16
```

Then pass it to ViralScan:

```bash
viralscan \
  -t t2g.txt -i index.idx -o output/ \
  -s1 R1.fastq.gz -s2 R2.fastq.gz \
  --host-filter starsolo \
  --host-index /path/to/star_hg38/
```

**What happens internally**

1. A new Snakemake rule (`host_filter`) runs before `kb_count`.
2. Filtered FASTQ files are written to the sample run directory:
   `{output}/{sample}/host_filtered/R1.fastq.gz` and `R2.fastq.gz`.
3. `kb_count` automatically uses those files instead of the originals — no
   further changes to your command are needed.
4. The original FASTQ files are never modified.
5. Pair synchronization is validated before and after filtering. Retained read
   IDs are written to `fragment_lineage.tsv.gz`, with aggregate counts in
   `host_filter_audit.tsv`.

When `--host-filter` is not supplied, no pre-subtraction is performed.

---

### Which option should I choose?

- **STARsolo** catches genome-level host reads, including intronic and
  intergenic reads, but requires substantial RAM and a pre-built genome index.
- Kallisto host subtraction is not exposed in v3 because its BUS output lacks
  exact source read IDs. Removing every fragment sharing a mapped CB–UMI can
  delete unrelated viral evidence.

For most 10x experiments, start with a combined host+virus reference and the
default host-conservative multimapping. Use STARsolo subtraction only when its
irreversible information loss is acceptable.

---

## Development

### How do I run the test suite?

```bash
PYTHONPATH=src python -m pytest tests/ -v
```

Network-dependent tests are gated by the `network` marker:

```bash
PYTHONPATH=src python -m pytest tests/ -v -m network
```

### How do I contribute?

1. Fork the repository on GitHub.
2. Create a feature branch.
3. Make your changes (follow the conventions in `CLAUDE.md`).
4. Run `ruff check . && ruff format . && pytest tests/`.
5. Open a pull request against `main`.
