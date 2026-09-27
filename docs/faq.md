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

There is now a **third** denominator, `pct_infected_comparable`, for comparing
runs that used different host-filtering strategies. `pct_infected_called` divides
by the called cells *of this run*, and cell calling happens after host
subtraction, so the denominator moves. Measured on one covid PBMC sample:

| strategy | called cells | Alphatorquevirus UMI | `pct_infected_called` |
|---|---:|---:|---:|
| no host filter | 143,243 | 1,167,103 | 56.64 % |
| `--host-filter starsolo` | 28,921 | 57,715 | 62.89 % |

Viral molecules fell 20.2x and the reported prevalence went *up*, because the
denominator collapsed 5.0x faster. Use `pct_infected_comparable` across
strategies; it uses an absolute host-UMI floor that host filtering cannot move.

### ViralScan found nothing. Is there really nothing there?

**Usually you cannot tell from the output alone — which is why every run now
writes `results/sensitivity.tsv` and states the limit in `summary.txt`.** Three
separate terms decide whether a virus that *is* present gets reported, and only
the first is measurable from inside a run:

1. **Depth.** Molecules arrive as a thinning Poisson process, so the abundance
   resolved with 95 % probability is about **3 molecules**. A routine 10x run
   quantifies 5–20 M molecules, giving an LOD95 of roughly 0.001–0.006 viral UMI
   per 10k host UMI. Depth is almost never the binding constraint: the three
   covid configurations above all landed in the `informative` band.

2. **k-mer capture.** Pseudoalignment needs an *exact* 31-mer match. A 90 bp
   fragment at 15 % divergence is captured with probability 0.32; at 20 %, 0.06;
   at 30 %, 0.001. **This is the term that decides your negative**, and it does
   not appear anywhere in a count matrix.

3. **Allocation survival.** The default `host-conservative` multimap method
   credits host-virus-ambiguous molecules *zero* to the virus.

So: read `informative_negative` in `results/sensitivity.tsv`. It is `false`
unless depth is sufficient **and** a k-mer capture term was *measured*. Since
capture cannot be measured without a control, that column is `false` on almost
every run — deliberately. To make a negative certifiable:

```bash
# plant a spike-in at a known abundance, then require it
viralscan ... --positive-control-gene SPIKEIN_gp1 \
              --positive-control-molecules 1000 \
              --require-positive-control
```

The recovered fraction is the capture term, and `results/positive_control.json`
reports it along with the sequence divergence it implies. Alternatively, align
the reads to the target directly with `viralscan evidence`.

A concrete case: the bundled 20-genome Torque teno virus panel shares **1.36 %**
of the 31-mer space of the 2,042 real human anellovirus genomes, and 85.8 % of
those genomes share *zero* 31-mers with it. A TTV negative from that panel is
not a statement about the sample. Build with
`viralscan build-ref --anellovirus` and the panel captures the whole
anellovirus sequence space (1.21x index inflation, 6 MB).

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

### Can we tell whether EBV is latent or lytic? Can we do this for other viruses?

**Yes, as a second layer, for nine viruses — but not by comparing gene totals, and
the result is only meaningful for four of them.**

Run it with `--gene-programs`. It writes `results/gene_program_summary.tsv` and
`results/gene_program_cells.tsv`, one row per cell, over the viruses the first
layer already detected.

**Why a naive version is wrong.** Summing latent genes and summing lytic genes
on the bundled EBV LCL run (`SRR12682296`, a cell line latently infected by
construction) gives 236,342 latent against 247,633 lytic — an aggregate ratio
of **1.15**. `EBNA-1`, which is expressed from every latent episome and must be
present in every infected cell, is 920 UMI, about 155x below `BHLF1`. That is
not biology: EBV's latent transcripts come from a region packed with nested and
antisense lytic ORFs, so reads cross-map both ways. Per-gene aggregate totals
are simply uninformative here.

**What ViralScan does instead.** It calls per cell, and requires breadth across
distinct **non-overlapping overlap groups** rather than a count of genes — in EBV
the whole latent EBNA locus (`EBNA-1`, `EBNA-2`, `EBNA-LP`) is one group, and
`BTRF1` shares a group with `BcLF1`, so neither pair is independent evidence.
Evidence comes from the *uniquely-placing* molecule layer, where `BZLF1` — the
canonical lytic marker — has **zero** molecules, which is the right answer for a
latent cell line, and `BARF1.2` has 13,668 where the multimap-allocated layer has
none.

> **Correction (2026-09-27).** `BARF1` is latent only in epithelial cancers (NPC, EBV-gastric; PMID 32708965) and `BaRF1.1` is the lytic ribonucleotide reductase, so neither is a latency marker in a B-cell line. Both have been removed from the catalogue. The per-marker figures and the 2,240 / 1,277 cell counts in this section were measured with them included. Re-measured without them on the same run: **895** cells latent on the uniquely-placing layer versus **856** on the allocated layer, so most of the apparent latent-sensitivity gain came from `BARF1.2` (PLAN `PROG-07`).

**What that buys, precisely.** Per-cell calling is directionally consistent on
both layers — the two never disagree in the dangerous direction (0 cells go
latent-on-unique to productive-on-allocated). The gain is **sensitivity**:
**2,240** cells called latent versus **1,277** on the allocated layer, because
1,263 fall to `indeterminate` there once cross-mapping has drained their latent
signal. `gene_program_summary.tsv` reports both so you can see this rather than
take it on trust.

**Which viruses, and how far to trust it.**

| | viruses |
|---|---|
| `panel_completeness=complete` | EBV, HHV-6A, HHV-7 |
| `panel_completeness=partial` | CMV, HSV-1, HSV-2, HHV-6B, VZV, KSHV |

CMV is partial for a different reason: single-cell HCMV latency shows no restricted latency programme but a late-lytic one at much lower levels (PMID 29535194), so marker presence cannot separate the states. For the other partial viruses the latency anchor set is too thin to support an absence
claim — HSV-1's only latency transcript is `LAT` — so
`latency_observable_in_rna` is `false` and the `latent` state is **unreachable by
construction**. Those rows can only ever read `productive` or `indeterminate`.
That is the honest answer, but it does mean an HSV-1 sample cannot be shown to
be latent with this tool.

**Two things it will not tell you.** It is a transcriptomic assay throughout: a
silent HIV provirus, a transcriptionally silent integrated HPV genome and the
HBV cccDNA pool all produce no reads, so they are **invisible, not latent**. And
`indeterminate` is not a negative — it means the virus was detected but no
programme met its threshold, which is a statement about breadth, not virology.

Already have a run? `viralscan rerun-programs --run-dir <dir>` works in place,
because layer 2 only reads the counts and layer 1's summary.

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

Viral gene IDs are resolved in two tiers, both in
`src/viralscan/virus_grouping.py`:

1. `VIRUS_NAME_MAP` — underscore/digit-boundary prefix match. This is strict on
   purpose: it is what stops `EPSTEIN_HHV4_BORF1` being read as Orf virus and
   `BUNYAMW_...` as Bunyavirus La Crosse.
2. `VIRUS_GENE_ID_ALIASES` — plain prefix match, consulted **only** when tier 1
   matches nothing. This exists for panel schemes that write the virus token and
   the gene token with no separator (`Ydvgp129`, `TTVgp1`, `HHV1gp00p39`,
   `HHV5wtgp045`), which tier 1 rejects by design. Being consulted second, it
   is strictly additive and can never rename a gene that already resolves.

If your prefix is in neither, the raw gene ID is reported — visible, but split
one row per gene, which also breaks `accession_breadth`, sibling cross-mapping
and `eve_risk` for that virus. Add it to `VIRUS_NAME_MAP` if the token is
underscore-delimited, or to `VIRUS_GENE_ID_ALIASES` if it is concatenated. A
token with no confident virus assignment should be left out on purpose: a wrong
name is worse than the raw-ID fallback.

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
