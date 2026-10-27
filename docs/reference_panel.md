# Reference Panel

ViralScan uses **195 viral GTF annotation files** distributed through Zenodo
under DOI `10.5281/zenodo.20112332`. It is required only for an index that has no
`<index>.build_manifest.json`; an index built by `viralscan build-ref` carries its
own viral gene list. When it is needed, fetch the panel once after installation:

```bash
viralscan data fetch
```

The command verifies the Zenodo checksum and unpacks the GTF files to
`~/.cache/viralscan/data/`. The annotations were generated from RefSeq
GenBank entries using the `Viral_GTF_maker.py` script
(whole-genome-as-gene model: one gene per chromosome, covering the full
replicon length).

---

## When you need this panel

Run `viralscan data fetch` to use ViralScan's bundled viral annotations or helper
code that expects the cached panel. Run it also for an index that has no
`<index>.build_manifest.json`, because the analysis step then reads the viral gene
list from the panel. The cache is not needed in three cases: the index has a build
manifest (every index from `viralscan build-ref`), every run supplies custom GTF
files with `--reference --gtf ...`, or a run uses `--ncbi-accession` to fetch a
reference on demand.

Check the cache location:

```bash
viralscan data fetch
ls ~/.cache/viralscan/data/
```

Use a different cache root when running on a shared filesystem:

```bash
viralscan data fetch --cache-dir /shared/viralscan-cache
```

The GTF files are written to `/shared/viralscan-cache/data/`. Quantification
reads the same cache with either:

```bash
viralscan run --cache-dir /shared/viralscan-cache ...
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

(host-control)=
## Host control and reference strategies

Host signal can reach viral features through shared k-mers. There are three
supported ways to handle it, and they trade different things. For routine
host-aware work, use the combined host + virus index with the default
`--multimap-method host-conservative`.

### Combined host + virus index (recommended)

Build one index from the host transcriptome and the viral accessions:

```bash
viralscan build-ref \
  --host human \
  --virus-accessions NC_045512.2 NC_002021.3 \
  --output ref_human_virus/ \
  --ncbi-email you@example.org
```

Host and virus compete in one k-mer space, so a fragment compatible with both
stays visible as a host–virus ambiguous molecule. Nothing is deleted: the
ambiguity appears in `multimap_evidence.tsv` and in the count layers
(`counts_host_viral_ambiguous`, `counts_host_viral_selected`), and the default
`host-conservative` method keeps mixed host–virus mass out of primary viral
estimates. The cost is that some viral molecules remain ambiguous and are
reported as such rather than as unique.

### STARsolo pre-subtraction (two-step)

`--host-filter starsolo --host-index PATH` maps reads to the full host genome
before viral quantification and passes only unmapped fragments to `kb count`.
This catches genome-level host reads, including intronic and intergenic reads,
but it irreversibly removes host-homologous and ambiguous fragments before the
viral pass. The viral-only reference never sees them and cannot report that
ambiguity, and the called-cell set changes, which is why
`pct_infected_comparable` exists. Use it when the extra conservative filter is
worth that information loss; it needs a pre-built STAR genome index and
substantial RAM and disk. Kallisto host subtraction is not exposed in v3 because
its BUS output lacks exact source read IDs; removing every fragment sharing a
mapped CB–UMI can delete unrelated viral evidence.

### Genome D-list masking

`--genome-dlist FASTA` passes the full host genome to `kb ref` as a D-list and
writes raw per-virus host-homology measurements. A kallisto D-list masks shared
k-mers by **exact match**, so it cannot see host sequence that has diverged even
slightly — and divergent host sequence is exactly what produces host-homologous
viral signal in the first place.

The warning below comes from one measured comparison, not a universal ranking.
On the same covid PBMC sample (`LUM-SJ-x213-g`, 1,203,332,091 reads), three ways
of removing host signal gave:

| | no filter | `--genome-dlist` | `--host-filter starsolo` |
|---|---:|---:|---:|
| reads reaching viral quant | 1,203,332,091 | 1,203,332,091 | 120,405,696 |
| `p_pseudoaligned` | 6.4 % | 4.5 % | 5.8 % |
| **`p_unique`** | **2.1 %** | **0.6 %** | **4.4 %** |
| quantified molecules | 21,613,840 | 8,404,326 | 5,308,302 |
| called cells | 143,243 | 28,922 | 28,921 |

D-list masking left most of the anellovirus artefact in place — the published
covid *Alphatorquevirus* signal was an artefact of poly-G no-signal reads
(≈90 %) and host-homologous reads (≈10 %); low-level divergent anellovirus is
not excluded — while cutting `p_unique` **3.5×** (2.1 % → 0.6 %). It destroyed
uniquely placed molecules to achieve very little.

The table has three arms only. It does not include a combined host + virus
*genome* reference arm, so it does not rank that option against the other two,
and no ordering beyond these measured three is claimed. Prefer
`--host-filter starsolo` over `--genome-dlist` when the goal is host subtraction;
use `--genome-dlist` only when the host genome FASTA is already on hand and its
host-homology measurements are the actual goal. With it,
`host_homology_annotations.tsv` retains maximum identity, query coverage, aligned
bases and best host target for every viral sequence, and the genome path and
SHA-256 are frozen in the manifest.

### Choosing

| Goal | Option |
|---|---|
| Routine host-aware quantification | Combined host + virus index, default `host-conservative` allocation |
| Extra conservative removal of genome-level host reads, accepting irreversible loss | `--host-filter starsolo` with a pre-built STAR index |
| Raw host-homology measurements for the reference | `--genome-dlist`, read as observations, not as a host-control performance claim |

---

## Papillomaviridae: what the HPV panel supports

The HPV entries carry named ORFs, so ViralScan can report HPV transcript
presence and an oncogene-versus-capsid contrast: E6/E7 (early-region
oncoproteins) against L1/L2 (late-region, productive infection). It does **not**
support per-genotype attribution from L1 or other late-region reads: L1 is the
most conserved coding region in the genus (pan-HPV PCR primers target it), so its
reads cross-map between genotypes and a genotype label on an L1 count is not
independent evidence of which type is present. Low-level calls for an off-type
HPV are likewise not evidence of that type; they reflect cross-mapping or the
library noise floor. E6/E7 are more type-divergent, but a cross-mapped count is
still not a transcript count. The panel is transcriptomic, so a silent
integrated genome is invisible rather than negative.

---

## Known limitations

- **HSV-1 latency (LAT).** HSV-1 LAT is modelled as the two-exon 6,391 nt transcript
  (660 nt + 5,731 nt). The stable major LAT intron is not polyadenylated, so sensitivity to latent
  HSV-1 in poly(A)-selected, 3' or 5' 10x libraries is not established either way.
  `panel_completeness` for HSV-1 stays `partial`.
- **`kmer_twin_of` exclusions.** A candidate excluded as `kmer_twin_of:<panel id>` (fewer than 5% of its
  31-mers lie outside the panel) is excluded on a genomic redundancy heuristic. That is not proof that
  the captured transcripts are equivalent: a small divergent region could dominate captured RNA.
  Catalogue-only additions do not change the sensitivity of the shipped index.
- **Blank `risk_class`.** A blank `risk_class` in the identity table means unknown, not "no risk"
  (see `docs/output_reference.md`).

---

## Adding custom references

Supplement the cached panel with your own annotation files. Pass the `--gtf` flag
(comma-separated list) together with `--reference`:

```bash
viralscan run --reference \
  --fasta custom_virus.fasta \
  --gtf   custom_virus.gtf \
  -o output/ \
  --r1 R1.fastq.gz --r2 R2.fastq.gz
```

Alternatively, use `--ncbi-accession` to fetch and build a reference for
any RefSeq nucleotide accession on the fly; multiple accessions are
comma-separated:

```bash
viralscan run --ncbi-accession NC_002021.3 -o output/ \
  --r1 R1.fastq.gz --r2 R2.fastq.gz --ncbi-email you@example.org
```

NCBI requires a contact email when a sequence has to be downloaded; pass
`--ncbi-email` or set `NCBI_EMAIL`. Accessions already in
`~/.cache/viralscan/ncbi/` are reused without a download.

For host-aware competitive mapping, use `viralscan build-ref` and include
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

Custom GTFs need stable and descriptive `gene_id` values. If the prefix is
not present in `VIRUS_NAME_MAP`, ViralScan may display the prefix rather than a
curated full virus name.

---

(anellovirus-reference)=
## Optional: extended Anelloviridae reference

The general panel above is enough to start. The rest of this page is an optional
reference for anellovirus screening; you do not need it for a routine run.

The Zenodo panel's Anelloviridae coverage is intentionally minimal (20 RefSeq
genomes). ViralScan ships a **~2,000-accession Anelloviridae reference** derived
from the complete set of human anellovirus genomes available in NCBI GenBank,
reconciled with the accession list of the Clareau lab's
[clareaulab/human\_anellovirus\_pangenome](https://github.com/clareaulab/human_anellovirus_pangenome)
repository (formerly `anellovirus_reference`; commit `3ed77e1`, 2026-09-18). That
repository clusters NCBI complete human Anelloviridae genomes (taxon 687329, 3,545
genomes) into 2,023 representatives across the modern genera (*Alpha-*, *Beta-*,
*Gamma-*, *Mem-*, *Samektorquevirus*, …); 2,020 of them are in the table, and the
GenBank sequences and their CDS annotation are re-fetched from NCBI, because the
repository's FASTA is pre-masked, its GTF has one whole-genome feature per genome, and
it publishes no licence.

### Building the anellovirus reference

```bash
viralscan build-ref --reference-panel anellovirus \
    --output anellovirus_ref/ \
    --ncbi-email you@example.org
```

This uses the bundled anellovirus FASTA when cached, or downloads the packaged
accession list from NCBI (cached under `~/.cache/viralscan/ncbi/`). It hard-masks
low-complexity regions with
`dustmasker` (windows 64 and 30, level 30), then masks remaining low-complexity
k-mer windows. The GTF uses the packaged CDS/exon catalogue where available.
Unannotated genomes get a whole-genome placeholder (`{accession}_geneN`).
The command runs `kb ref` unless `--no-kb-ref` is supplied.

Both native builders share the masking gate, GTF emission, and index receipt
boundary. `reference_manifest.json` records `source_adapters` (`ensembl_cdna`,
`ncbi_genbank`, or `local_fasta`) and whether a successful index command applied
the genome D-list (`build_receipt.genome_dlist_applied`). These software receipts
do not establish scientific performance or licence clearance.

| Flag | Effect |
|------|--------|
| `--no-mask` | Skip `dustmasker` hard-masking of the viral sequence (every build path); the k-mer gate then allows a few low-complexity k-mers instead of none |
| `--cluster` | Run `cd-hit-est` at 95% identity after masking (off by default — the accession list is already clustered at 95 % ANI) |
| `--no-kb-ref` | Skip `kb ref`; produce only FASTA + GTF |
| `--ncbi-api-key KEY` | NCBI API key for higher download throughput |
| `--cache-dir PATH` | Root for NCBI download cache (default `~/.cache/viralscan/`) |

### Prerequisites

| Tool | Version | Purpose |
|------|---------|---------|
| `dustmasker` | any (NCBI BLAST+) | Hard-masking of low-complexity regions |
| `cd-hit-est` | ≥ 4.8.1 | Redundancy clustering (optional) |
| `kb` (kb-python) | ≥ 0.27 | kallisto index construction |

Requested masking, clustering, or indexing stops if its tool is missing.
Use `--no-mask`, omit `--cluster`, or use `--no-kb-ref` to skip that operation.
Install masking and clustering tools via conda:

```bash
conda install -c bioconda blast cd-hit
```

### Virus naming

Detected genes (e.g. `NC_002076.2_gene1`) are resolved to genus-level labels
(*Alphatorquevirus*, *Betatorquevirus*, *Gammatorquevirus*, …) via the packaged
`anellovirus_accessions.tsv` + `merged_name_map()`. The legacy `TTV`-prefix
gene_ids from the Zenodo panel still resolve correctly.

### Citation / provenance

The accession list is derived from
[clareaulab/human\_anellovirus\_pangenome](https://github.com/clareaulab/human_anellovirus_pangenome)
at commit `3ed77e1`. The repository publishes no licence, paper or DOI, so cite the URL and commit,
and note that ViralScan redistributes none of its sequences or annotations (accessions only).

---

## Maintenance, integrity and comparison pages

These pages are kept for citation and are not in the sidebar:

- [Reference updates and manifest compatibility](reference_maintenance.md): a
  frozen reference is never changed in place; identity is the checksum; manifest
  `schema_version` rules.
- [Panel integrity checks](panel_integrity.md): GTF corpus checks, catalogue
  structure, near-duplicate and host-homology screens at build time.
- [Workflow comparison](workflow_comparison.md): the five workflow families the
  validation protocol compares, with their defined strengths and limits.
- [ViralScan 3.0 molecule-counting contract](v3_counting_contract.md): how
  corrected molecules become gene counts, and the conservation invariants.

---

## Cache location

By default, `viralscan data fetch` writes to `~/.cache/viralscan/data/`.
Set `VIRALSCAN_CACHE` or use `viralscan data fetch --cache-dir PATH` to
populate another cache root. Quantification reads that cache when you pass
`--cache-dir PATH` or set the same `VIRALSCAN_CACHE` environment variable.
