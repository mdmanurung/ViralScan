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

It is a permutation test that asks: are viral-infected cells over-represented
among each other's nearest neighbours in the UMAP embedding? A low p-value
indicates spatial clustering of infected cells beyond what would be expected
after shuffling the viral labels.

### How is `pct_infected` calculated?

```
pct_infected = (cells with any UMI assigned to this detected virus) / (total cells in the count matrix) × 100
```

`--detection-threshold` (default 1) is a sample-level threshold for calling a
virus detected. It does not change the per-cell infected-cell count.

### What units is `umi_per_10k` in?

It is the total viral UMI normalised to 10 000 total UMI (CPM-equivalent):

```
umi_per_10k = total_viral_umi / total_all_umi × 10 000
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
with both host and viral genes, its ambiguous mass is not assigned to the viral
gene for primary counts. This reduces false positives while preserving
diagnostic evidence in `results/multimap_evidence.tsv`.

The legacy `equal` method is still available with `--multimap-method equal` for
backward-compatible comparisons. For even stricter calling, use
`--multimap-primary-call unique-only`.

---

### How do I use host pre-subtraction to reduce false positives?

Host pre-subtraction is an optional advanced filter. It maps reads to the host
genome or transcriptome **before** viral quantification and discards reads that
align. The remaining reads are then passed to `kb count`.

For routine host-aware analysis, prefer the combined host+virus reference plus
the default `--multimap-method host-conservative`. Use pre-subtraction when you
want an extra conservative filter or need to remove host-aligned reads before
viral quantification.

Two aligners are supported via `--host-filter`:

| Aligner | Flag value | What it does |
|---|---|---|
| STARsolo | `starsolo` | Full genome alignment; unmapped reads collected from STAR's `--outReadsUnmapped Fastx` output |
| kallisto | `kallisto` | Pseudo-alignment against a host cDNA index; unmapped read pairs identified via BUS file subtraction |

**Option A — STARsolo (most comprehensive host subtraction)**

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

**Option B — kallisto (faster; stays within the kb ecosystem)**

Requires a kallisto index built from a host cDNA FASTA. Download the host cDNA
FASTA from Ensembl or another trusted source, then build the index once:

```bash
kallisto index -i host.idx Homo_sapiens.GRCh38.cdna.all.fa.gz
```

Then pass it to ViralScan:

```bash
viralscan \
  -t t2g.txt -i index.idx -o output/ \
  -s1 R1.fastq.gz -s2 R2.fastq.gz \
  --host-filter kallisto \
  --host-index host.idx
```

**What happens internally**

1. A new Snakemake rule (`host_filter`) runs before `kb_count`.
2. Filtered FASTQ files are written to the sample run directory:
   `{output}/{sample}/host_filtered/R1.fastq.gz` and `R2.fastq.gz`.
3. `kb_count` automatically uses those files instead of the originals — no
   further changes to your command are needed.
4. The original FASTQ files are never modified.

When `--host-filter` is not supplied, no pre-subtraction is performed.

---

### Which option should I choose?

- **STARsolo** catches more host reads (genome-level, including intronic and
  intergenic reads) but requires ~30 GB of RAM and a pre-built genome index.
- **kallisto** is faster and uses less memory, but only subtracts reads whose
  cDNA sequence pseudo-aligns to an annotated host transcript; reads from
  unannotated loci or introns are not removed.

For most 10x experiments, start with a combined host+virus reference and the
default host-conservative multimapping. If you add pre-subtraction, choose
kallisto for speed and STARsolo when genome-level host depletion matters.

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
