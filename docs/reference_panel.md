# Reference Panel

ViralScan uses **195 viral GTF annotation files** distributed through Zenodo
under DOI `10.5281/zenodo.20112332`. Fetch the panel once after installation:

```bash
viralscan data fetch
```

The command verifies the Zenodo checksum and unpacks the GTF files to
`~/.cache/viralscan/data/`. These annotations were generated from RefSeq
GenBank entries using the `Viral_GTF_maker.py` script
(whole-genome-as-gene model: one gene per chromosome, covering the full
replicon length).

---

## Anelloviridae — extended reference

The Zenodo panel's Anelloviridae coverage is intentionally minimal (20 RefSeq
genomes). ViralScan ships a **~2,000-accession Anelloviridae reference** derived
from the complete set of human anellovirus genomes available in NCBI GenBank,
reconciled with the curated accession/taxonomy table from the
[clareaulab/anellovirus\_reference](https://github.com/clareaulab/anellovirus_reference)
project (Lareau *et al.*, 2023). That repository curates representative genomes
across the modern Anelloviridae genera (*Alpha-*, *Beta-*, *Gamma-*, *Mem-*,
*Samektorquevirus*, …) and served as the accession + taxonomy source; the
underlying GenBank sequences are re-fetched from NCBI to ensure license
clarity.

### Building the anellovirus reference

```bash
viralscan build-ref --anellovirus \
    --output anellovirus_ref/ \
    --ncbi-email you@example.org
```

This downloads all ~2,000 accessions from NCBI (cached under
`~/.cache/viralscan/ncbi/`), hard-masks low-complexity regions with
`dustmasker` (`-window 64 -level 30`), regenerates a whole-genome GTF
(`gene_id "{accession}_geneN"`), and optionally runs `kb ref` to build
a kallisto index.

| Flag | Effect |
|------|--------|
| `--no-mask` | Skip `dustmasker` hard-masking |
| `--cluster` | Run `cd-hit-est` at 95% identity after masking (off by default — clareaulab set is already clustered) |
| `--no-kb-ref` | Skip `kb ref`; produce only FASTA + GTF |
| `--ncbi-api-key KEY` | NCBI API key for higher download throughput |
| `--cache-dir PATH` | Root for NCBI download cache (default `~/.cache/viralscan/`) |

### Prerequisites

| Tool | Version | Purpose |
|------|---------|---------|
| `dustmasker` | any (NCBI BLAST+) | Hard-masking of low-complexity regions |
| `cd-hit-est` | ≥ 4.8.1 | Redundancy clustering (optional) |
| `kb` (kb-python) | ≥ 0.27 | kallisto index construction |

Both `dustmasker` and `cd-hit-est` are gracefully skipped with a warning
if they are not on `PATH`. Install via conda:

```bash
conda install -c bioconda blast cd-hit
```

### Virus naming

Detected genes (e.g. `NC_002076.2_gene1`) are resolved to genus-level labels
(*Alphatorquevirus*, *Betatorquevirus*, *Gammatorquevirus*, …) via the packaged
`anellovirus_accessions.tsv` + `merged_name_map()`. The legacy `TTV`-prefix
gene_ids from the Zenodo panel still resolve correctly.

### Citation / provenance

The accession/taxonomy table is derived from
[clareaulab/anellovirus\_reference](https://github.com/clareaulab/anellovirus_reference).
Please cite the clareaulab repository and the underlying Lareau *et al.* study
when publishing results that use this reference.

---

## When you need this panel

Run `viralscan data fetch` if you want to use ViralScan's bundled viral
annotations or helper code that expects the cached panel. You do not need the
cache when every run supplies custom GTF files with `--reference -gtf ...` or
when a run uses `--ncbi-accession` to fetch a reference on demand.

Check the cache location:

```bash
viralscan data fetch
ls ~/.cache/viralscan/data/
```

Use a different cache root when running on a shared filesystem:

```bash
viralscan data fetch --cache-dir /shared/viralscan-cache
```

The GTF files will be under `/shared/viralscan-cache/data/`. Use the same
cache during quantification with either:

```bash
viralscan --data-cache-dir /shared/viralscan-cache ...
```

or by setting `VIRALSCAN_CACHE=/shared/viralscan-cache` for both fetch and run
commands.

---

## Included viruses

The panel covers the following virus families and genera (non-exhaustive):

- **Adenoviridae** — Human adenovirus
- **Anelloviridae** — *Alphatorquevirus*, *Betatorquevirus*, *Gammatorquevirus*,
  *Memtorquevirus*, *Samektorquevirus* (Torque teno virus spp.); ~2,000 accessions
  via `viralscan build-ref --anellovirus`; 20 RefSeq genomes in the Zenodo panel
- **Arenaviridae** — Lassa virus, Lymphocytic choriomeningitis virus
- **Bunyaviridae / Hantaviridae** — Hantaan virus, Seoul virus, Bunyamwera virus, La Crosse virus
- **Caliciviridae** — Norwalk virus, Sapporo virus
- **Coronaviridae** — SARS-CoV, MERS-CoV, Human coronavirus 229E / HKU1 / NL63 / OC43
- **Filoviridae** — Ebolavirus, Lake Victoria marburgvirus
- **Flaviviridae** — Dengue virus, Hepatitis C virus, West Nile virus, Yellow fever virus, Zika virus
- **Herpesviridae** — EBV, CMV, HSV-1/2, HHV-6/7/8, VZV
- **Orthomyxoviridae** — Influenza A, B, C
- **Papillomaviridae** — Human papillomavirus 1, 2, 16, 18
- **Paramyxoviridae** — Measles virus, Mumps virus, Nipah virus, Hendra virus
- **Parvoviridae** — Human parvovirus B19
- **Picornaviridae** — Poliovirus, Rhinovirus, Coxsackievirus, Echovirus, Enterovirus 68/70
- **Polyomaviridae** — BK polyomavirus, JC polyomavirus, Merkel cell polyomavirus
- **Poxviridae** — Vaccinia virus, Cowpox virus, Monkeypox virus, Variola virus
- **Reoviridae** — Rotavirus A/B/C, Banna virus
- **Rhabdoviridae** — Rabies virus, Australian bat lyssavirus
- **Togaviridae** — Chikungunya virus, Rubella virus, Alphavirus spp.
- **Other** — Hepatitis A/B/D/E, Adeno-associated virus (AAV), and more

---

## Adding custom references

You can supplement the cached panel with your own annotation files using the
`-gtf` flag (comma-separated list) together with `--reference`:

```bash
viralscan --reference \
  -fasta custom_virus.fasta \
  -gtf   custom_virus.gtf \
  -o output/ \
  -s1 R1.fastq.gz -s2 R2.fastq.gz
```

Alternatively, use `--ncbi-accession` to fetch and build a reference for any
RefSeq nucleotide accession on-the-fly.

For host-aware competitive mapping, prefer `viralscan build-ref` and include
both the host transcriptome and viral accessions in one index. The default
`--multimap-method host-conservative` then keeps host-virus ambiguous signal
out of primary viral counts:

```bash
viralscan build-ref \
  --host human \
  --virus-accessions NC_045512.2 NC_002021.3 \
  --output ref_human_virus/ \
  --ncbi-email you@example.org
```

---

## GTF format

Each panel GTF file follows the RefSeq GTF convention:

```
NC_001477.1  RefSeq  gene  1  10735  .  +  .  gene_id "DENV_DV1_gp1"; ...
```

The `gene_id` prefix (e.g. `DENV`) is used by `detection.py` to look up the
human-readable virus name in `VIRUS_NAME_MAP` (`constants.py`).

Anellovirus GTFs produced by `viralscan build-ref --anellovirus` use accession-
based gene_ids (e.g. `NC_002076.2_gene1`). These are resolved to genus labels
via the bundled `anellovirus_accessions.tsv` + `merged_name_map()`.

For custom GTFs, use stable and descriptive `gene_id` values. If the prefix is
not present in `VIRUS_NAME_MAP`, ViralScan may display the prefix rather than a
curated full virus name.

---

## Cache location

By default, `viralscan data fetch` writes to `~/.cache/viralscan/data/`.
Set `VIRALSCAN_CACHE` or use `viralscan data fetch --cache-dir PATH` to
populate another cache root. Quantification reads that cache when you pass
`--data-cache-dir PATH` or set the same `VIRALSCAN_CACHE` environment variable.
