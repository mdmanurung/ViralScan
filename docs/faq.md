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

### Can I process multiple samples in one run?

Yes — provide comma-separated paths to `-s1` and `-s2`:

```bash
viralscan -t t2g.txt -i index.idx -o output/ \
  -s1 A_R1.fastq.gz,B_R1.fastq.gz \
  -s2 A_R2.fastq.gz,B_R2.fastq.gz
```

### The UMAP step takes a very long time. Can I skip it?

Yes — omit `--umap` (it is off by default).  The detection and reporting
steps run without it.

---

## Results

### What does `viral_neighbor_enrichment` measure?

It is a Fisher's exact test that asks: are viral-infected cells
over-represented among each other's nearest neighbours (in gene-expression
space)?  A low p-value indicates spatial clustering of infected cells beyond
what would be expected by chance.

### How is `pct_infected` calculated?

```
pct_infected = (cells with ≥ detection_threshold viral UMI) / (total cells passing QC) × 100
```

Adjust `--detection-threshold` (default 1) to change the minimum UMI count
required to call a cell infected.

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
