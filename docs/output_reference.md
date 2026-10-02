# Output Reference

ViralScan writes one sample directory under the path passed to `--output / -o`.
The sample directory is inferred from the R1 FASTQ filename before the first
underscore. For `sample_R1.fastq.gz`, the run directory is `output/sample/`.

---

## Directory layout

```
output/
└── sample/
    ├── config.yaml
    ├── summary.txt
    ├── report.html
    ├── log/
    │   ├── analysis.txt
    │   ├── detection.done
    │   ├── kb.done
    │   ├── multimap.done
    │   └── umap.done
    ├── kb-python/
    │   ├── counts_unfiltered/
    │   │   ├── adata.h5ad
    │   │   └── adata_multimap.h5ad
    │   ├── output.bus
    │   ├── output.resolved.sorted.bus
    │   ├── output.resolved.sorted.bus.txt
    │   ├── run_info.json
    │   └── ...
    ├── plots/
    │   ├── <virus>_histogram.png
    │   ├── SuperExpressor_<virus>.png
    │   ├── umap_binary.html
    │   └── umap_continuous.html
    └── results/
        ├── virus_identity.tsv
        ├── viral_summary.tsv
        ├── sensitivity.tsv
        ├── positive_control.json
        ├── gene_program_summary.tsv
        ├── gene_program_cells.tsv
        ├── per_cell_viral.tsv
        ├── multimap_evidence.tsv
        ├── cell_type_enrichment.tsv
        └── reference_provenance.json
```

`reference_provenance.json` records the viral reference used (index/t2g/GTF,
technology, multimap settings, and the viral accessions in the reference and
detected) so results are traceable to their annotation.
`gene_program_*.tsv` are present only when `--gene-programs` is supplied; see
`results/gene_program_summary.tsv` below.
`cell_type_enrichment.tsv` is present only when `--cell-types` is supplied.
`multimap_evidence.tsv` is present only when multimapping is enabled.
UMAP files are present only when `--umap` is supplied. `host_filtered/` is
present only when `--host-filter` is supplied.

`virus_identity.tsv` is written by every run; see
`results/virus_identity.tsv` below.

The raw `output.bus` is not a v3 counting input. ViralScan barcode-corrects it
when a whitelist is supplied, sorts it with bustools, and retains the resolved
sorted BUS plus text representation used for molecule counting and read lineage.

The v3 STAR filter also writes `host_filter_audit.tsv` with input, retained,
and removed fragment totals, plus `fragment_lineage.tsv.gz` containing each
retained exact read ID and its `host_unmapped` decision. ViralScan validates
mate synchronization before and after filtering.

---

## `run_complete.json`

Written at the run root (beside `run_manifest.json`) only after every sample has finished, and again after `rerun-multimap`. It is the completion marker: a run without it was interrupted or is still in progress. Fields: `schema_version`, `run_fingerprint` (must equal the manifest's), `viralscan_version`, `samples`, `completed_at`, and `artifacts`, a map of `<sample>/<path>` to sha256 for the narrow set `results/viral_summary.tsv`, `results/virus_identity.tsv`, `results/multimap_evidence.tsv` and `kb-python/counts_unfiltered/adata_multimap.h5ad` (each only if present).

The marker is deleted whenever a run starts, resumes or overwrites, and re-computed by in-place commands (`rerun-programs`, `hostresponse`) when one existed. New runs set `completion_marker: true` in `run_manifest.json`; `viralscan validate-run` then treats a missing marker, a fingerprint mismatch or an artifact hash mismatch as an error. Runs whose manifest lacks the field only get a warning for a missing marker. Schema: `schemas/v3/run_complete.schema.json`.

---

## `virus_identity.tsv`

Tab-separated, one row per gene of the index: whether it is viral and which
virus it belongs to. The `analysis` step builds it from the index t2g, whose
fifth column is each transcript's genome accession, joined to the packaged
virus catalogue.

| Column | Description |
|--------|-------------|
| `gene_id` | Gene ID as in the index |
| `genome_accession` | t2g column 5: the genome (GTF seqname) of the gene; empty for a pre-v3 t2g |
| `status` | `catalogued` (accession in the catalogue), `uncatalogued` (in the viral GTF set, not catalogued; the run warns), `host`, or `legacy_prefix` (pre-v3 t2g, named by prefix maps) |
| `viral` | `true` / `false` |
| `virus_key` | One virus: `taxid:<n>`; `taxid:<n>\|<strain>` for one isolate of a segmented virus whose taxid several isolates share; `genus:<genus>` for anelloviruses; `accession:<acc>` for uncatalogued genomes; `name:<name>` for the legacy fallback |
| `virus_name` | Display name of the virus: the curated name, else the NCBI organism |
| `organism`, `taxid`, `species`, `genus`, `family` | NCBI organism and taxonomy of the genome |
| `segment`, `strain` | Segment and strain/isolate of the genome, when recorded |
| `sibling_group` | Near-identical viruses that allocation can move molecules between (e.g. `HHV-4` for EBV types 1 and 2) |
| `risk_class` | `eve` for families with endogenous or ubiquitous-commensal risk |
| `role` | `decoy` for lab-contaminant decoy genomes |

A run with no viral gene in its index stops at this step.

## `viral_summary.tsv`

Tab-separated, one row per detected virus.

| Column | Description |
|--------|-------------|
| `virus_name` | Human-readable virus name |
| `viral_molecules_total_est` | Unique viral molecules plus allocated ambiguous molecule mass |
| `infected_cells` | Legacy-named schema field: cells with nonzero selected-method candidate molecule support after the sample-level reporting threshold; not confirmed infection |
| `total_cells` | Total cells in the count matrix (**all** barcodes) |
| `pct_infected` | `infected_cells / total_cells × 100` (all-barcode denominator) |
| `viral_molecules_per_10k_est` | Viral molecule estimate divided by full-matrix molecule estimate × 10,000 |
| `n_called_cells` | Number of **called** cells (real, non-empty droplets) — see cell-calling below |
| `infected_called` | Candidate-support cells restricted to the called-cell set |
| `pct_infected_called` | `infected_called / n_called_cells × 100`; a called-cell candidate-support rate, not a biological infection rate. **Within-run only** — see below |
| `infected_comparable` | Candidate-support cells over the strategy-independent denominator |
| `n_comparable_cells` | Barcodes clearing an absolute host-UMI floor (200 molecules), intersected with the called set |
| `pct_infected_comparable` | `infected_comparable / n_comparable_cells × 100`. Use this to compare runs that used different host-filtering strategies |

**Two denominators.** The all-barcode `pct_infected` field is diluted by empty
droplets; `pct_infected_called` uses only the declared called-cell set and is the
appropriate denominator for a per-cell candidate-support rate. The field names
are retained for schema compatibility and do not establish infection. The v3
`auto` default uses an external called-cell list when one is supplied and
otherwise runs DropletUtils EmptyDrops. The approximate knee caller is explicit
only and is never a fallback. It is sensitivity-only and must not supply
reported numbers: on real libraries it calls empty droplets as cells. With `cell_calling=none`, the `*_called` columns
equal the all-barcode values.

**Three denominators, and why the middle one is not comparable across runs.**
`pct_infected_called` divides by the called cells *of this run*, and cell calling
runs after host subtraction, so the denominator moves with the strategy. Measured
on the same covid PBMC sample:

| strategy | called cells | Alphatorquevirus UMI | `pct_infected_called` |
|---|---:|---:|---:|
| no host filter | 143,243 | 1,167,103 | 56.64 % |
| `--host-filter starsolo` | 28,921 | 57,715 | 62.89 % |

Viral molecules fell **20.2×** and the reported prevalence rose, because the
denominator collapsed 5.0× faster than the numerator. Comparing those two runs
via `pct_infected_called` inverts the result. `pct_infected_comparable` uses an
absolute host-UMI floor that host filtering cannot move, so it is the field to
use across strategies; `pct_infected_called` remains the within-run primary.

**Count layer.** V3 summaries use `adata.X`, the complete selected-method
molecule matrix. `counts_unique` and `counts_ambiguous_allocated` are disjoint
and sum to `X`. Nonzero molecule support is candidate evidence; biological
interpretation requires calibrated evidence and may still require orthogonal
confirmation. The separate read-level workflow supplies diagnostics rather than
an automatic infection call.

---

## `sensitivity.tsv`

Tab-separated, one row per virus, written on **every** run. Answers the question
a zero otherwise cannot: *is there nothing there, or did we not look hard enough?*

| Column | Description |
|--------|-------------|
| `virus_name` | Virus the row describes |
| `observed_molecules` | Molecules attributed to the virus; `0` means a negative |
| `detection_threshold` | The sample-level UMI gate that decided the call |
| `capture` | Fraction of true viral molecules surviving exact k-mer matching |
| `capture_measured` | `true` only if the capture term came from a positive control, not a default |
| `lod95_per_10k` | Estimated viral molecules per 10k host molecules at which the virus would be reported with 95 % probability |
| `lod95_molecules` | Expected true molecules at that limit — always ≈ 3 × `detection_threshold` |
| `lod_interpretation` | `informative` / `adequate` / `shallow` / `insufficient-depth` |
| `depth_sufficient` | Molecular depth alone resolves `lod95_per_10k` |
| `informative_negative` | `depth_sufficient` **and** `capture_measured` |
| `expected_molecules_at_1_per_10k` | Expected observed molecules for a virus at 1 estimated molecule per 10k host molecules (cf. `viral_molecules_per_10k_est`) |
| `p_detect_at_1_per_10k` | Probability of clearing the threshold at that abundance |
| `p_zero_at_lod95` | ≈ 0.05 by construction; reported so the arithmetic is checkable |
| `notes` | Why the LOD is a bound rather than an estimate, when a row is zero, etc. |

**`informative_negative` is the column that matters, and it is false almost
always.** Depth is not the limiting term in practice: the three real covid
configurations produced LOD95 values of 0.0003–0.0056 estimated viral molecules per 10k host molecules, all in
the `informative` band, and the covid samples called SARS-CoV-2 = 0 at 21.6 M
quantified molecules. What cannot be measured from inside a run is the k-mer
**capture** term, which falls to 0.32 at 15 % sequence divergence and 0.06 at
20 %. Without a measured capture term a negative cannot be certified at any
depth — which is why `informative_negative` requires both conditions.

Depth is the sum of the count matrix, **not** raw reads: only quantified
molecules can be detected, and the two differ substantially — in this repo's own
benchmarks, pseudoaligned-read counts exceed quantified-molecule counts by
roughly 2x, so substituting reads for molecules would understate the LOD95 by
about the same factor. The Poisson floor was verified by molecule-level
downsampling of the bundled EBV LCL run: P(detect) stayed 1.0000 down to 1,270
downsampled reads and first reached 0 at 127.

## `positive_control.json`

Written on every run. Present so a negative can be audited.

| Field | Description |
|-------|-------------|
| `status` | `not-configured` / `measured` / `failed` / `over-recovered` / `gene-not-in-reference` |
| `certifies_negatives` | `true` only for `measured` |
| `gene`, `expected_molecules`, `observed_molecules`, `capture` | The recovery ratio |
| `implied_divergence` | Per-base divergence whose capture matches, by bisection; `null` when unidentifiable |

Supply a control with `--positive-control-gene` and `--positive-control-molecules`
(both required together). `failed` means the planted control was invisible, which
makes every negative in that run uninterpretable. `over-recovered` means more was
recovered than planted — the control is not spike-in-specific — and is treated as
no measurable loss, which is the optimistic direction.

Add `--require-positive-control` to fail the run outright when nothing is
detected and no capture term could be measured.

---

## `gene_program_summary.tsv` and `gene_program_cells.tsv`

Second layer, over the viruses the first layer detected. Present only when
`--gene-programs` is supplied.

### What it is for, and why a naive version is wrong

The obvious implementation is to sum the latent genes, sum the lytic genes, and
compare. On the bundled EBV LCL run (`SRR12682296`), which is latently
infected by construction, that gives:

| | UMI |
|---|---:|
| LATENT (10 genes) | 236,342 |
| LYTIC (12 genes) | 247,633 |

A latent:lytic aggregate ratio of **1.15** in a cell line defined by latency is
not biology. `EBNA-1` — expressed from every latent episome, so present in every
infected cell — is 920 UMI, ~155× below `BHLF1`. The cause is pervasive
overlapping-ORF cross-mapping: EBV's latent transcripts are transcribed from a
region densely packed with nested and antisense lytic ORFs, so reads cross-map
in both directions. **Per-gene aggregate totals are uninformative**, and no
amount of care applied to them recovers an answer.

> **Correction (2026-09-27).** `BARF1` is latent only in epithelial cancers (NPC, EBV-gastric; PMID 32708965) and `BaRF1.1` is the lytic ribonucleotide reductase, so neither is a latency marker in a B-cell line. Both have been removed from the catalogue. The per-marker figures and the 2,240 / 1,277 cell counts in this section were measured with them included. Re-measured without them on the same run: **895** cells latent on the uniquely-placing layer versus **856** on the allocated layer, so most of the apparent latent-sensitivity gain came from `BARF1.2` (PLAN `PROG-07`).

The per-marker breakdown on that same run shows the mechanism directly:

| marker | programme | uniquely-placing | multimap-allocated |
|---|---|---:|---:|
| `BARF1.2` | latent | **13,668** | 0 |
| `BNLF2a` | latent | **0** | 49,662 |
| `BNLF2b` | latent | **0** | 45,108 |
| `BZLF1` | productive | **0** | 9,308 |
| `BMRF1` | productive | **0** | 45,005 |
| `BcLF1` | productive | **6,609** | 59 |

The uniquely-placing layer puts **zero** molecules on `BZLF1` — the canonical
lytic marker — which is the correct answer for a latent cell line, and puts
13,668 on `BARF1.2` where the allocated layer put none.

Two defences, and it is worth being precise about what each buys:

1. **Breadth in overlap groups, never a gene count.** Overlap groups are
   computed in `extras/build_gene_programs.py` by exonic interval intersection.
   In EBV the whole latent EBNA locus (`EBNA-1`, `EBNA-2`, `EBNA-LP`) is one
   group, and `BTRF1` shares a group with `BcLF1` — so detecting one is not
   independent evidence for the others. This is what removes the need for an
   aggregate comparison at all.
2. **Uniquely-placing evidence.** Per-cell breadth calling is directionally
   consistent on *both* layers: on the real run the two never disagree in the
   dangerous direction (**0** cells go latent-on-unique to
   productive-on-allocated). What the unique layer buys is **sensitivity** —
   **2,240** cells called latent versus **1,277** on the allocated layer,
   because 1,263 cells fall to `indeterminate` there once cross-mapping has
   drained their latent signal onto lytic ORFs. The failure mode is lost
   sensitivity, not an inverted call.

A third defence is about honesty rather than arithmetic:

3. **No absence claim from a thin anchor set.** A virus whose latency side is a
   single transcript (HSV-1's `LAT`) cannot support "not detected, therefore
   latent", so `latency_observable_in_rna` is false and `latent` is unreachable.

### States

| State | Meaning |
|---|---|
| `productive` | ≥ `min_breadth` distinct non-overlapping productive overlap groups carry uniquely-placing molecules |
| `latent` | ≥1 latent overlap group, productive below `min_breadth`, **and** `latency_observable_in_rna` |
| `mixed` | both of the above in the same cell |
| `indeterminate` | the virus was detected but no programme met its threshold — not a negative |
| `not_applicable` | no programme model exists for this virus (summary rows only) |

### Coverage and honesty fields

| Column | Description |
|---|---|
| `panel_completeness` | `complete` (EBV, HHV-6A, HHV-7 — a real latency *and* reactivation split), `partial` (CMV, HSV-1/2, HHV-6B, VZV, KSHV), `not_applicable`. CMV is partial because single-cell latency mirrors a low-level late-lytic programme (PMID 29535194), so no marker's presence separates the states |
| `latency_observable_in_rna` | Whether a `latent` call is reachable. `false` ⇒ `n_cells_latent` and `n_cells_mixed` are 0 **by construction** |
| `evidence_layer` | Always `counts_unique_viral` |
| `min_breadth` | The `--programme-min-breadth` used |
| `productive_breadth_median` / `latent_breadth_median` | Breadth on the unique layer |
| `n_cells_latent_selected_layer` / `n_cells_productive_selected_layer` | What the same rule would have called on the multimap-allocated layer — the honest comparison |
| `selected_*_breadth_median` | Breadth on the allocated layer, for reference |
| `layer1_molecules` | Layer 1's molecule total, so the two layers need not be joined by hand |
| `caveat` | Why a row is weaker than it looks |

### Scope

Only viruses **detected by layer 1** appear. A virus layer 1 did not call gets
no programme row — that is a detection-limit problem (`sensitivity.tsv`), not
something layer 2 can repair. This is a transcriptomic assay throughout:
DNA-level latency (a silent HIV provirus, a transcriptionally silent integrated
HPV genome, the HBV cccDNA pool) produces no reads and is **invisible, not
latent**.

---

## `per_cell_viral.tsv`

Tab-separated, one row per cell × detected virus combination.

| Column | Description |
|--------|-------------|
| `barcode` | Cell barcode |
| `virus_name` | Virus name |
| `viral_molecules_total_est` | Viral molecule estimate for this cell from the selected-method matrix |
| `molecules_total_est` | Full selected-method molecule estimate for this cell |
| `viral_fraction` | Viral molecule estimate divided by full molecule estimate |
| `is_called_cell` | Boolean flag: whether the barcode is in the primary called-cell denominator (see cell-calling above) |

---

## `report.html`

A self-contained HTML file with:

- Run metadata (date, parameters, sample names)
- QC summary table
- Per-virus detection table with normalised metrics
- Multimapping evidence table with unique and ambiguous viral support
- Embedded histogram plots (base64 PNG)
- Interpretation guidance

UMAP clustering annotations are written to the interactive UMAP HTML files
when `--umap` is supplied.

Open in any modern browser — no internet connection required.

---

## `cell_type_enrichment.tsv`

Tab-separated, one row per detected virus and labeled cell type. Written only
when `--cell-types cell_types.csv` is supplied.

| Column | Description |
|--------|-------------|
| `virus` | Virus name |
| `cell_type` | Cell-type label from the CSV |
| `n_infected` | Legacy-named field: labeled cells with candidate molecule support in this cell type, using the selected-method matrix |
| `n_total` | Total labeled cells of this type |
| `pct` | `n_infected / n_total × 100` |
| `OR` | One-sided Fisher exact odds ratio |
| `pvalue` | Raw Fisher exact p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |

Input CSV requirements:

```csv
barcode,cell_type
AAACCCAAGAGT-1,T cell
AAACCCAGTGCA-1,Monocyte
```

Barcodes must match `adata.obs_names`. If no barcodes overlap, ViralScan skips
the enrichment table and logs a warning.

---

## `multimap_evidence.tsv`

Tab-separated, one row per viral gene in the reference. This table is additive:
it does not replace `viral_summary.tsv` or change its default schema.

| Column | Description |
|--------|-------------|
| `virus_name` | Human-readable virus name or accession fallback |
| `gene_id` | Viral gene/accession ID |
| `viral_molecules_unique` | Integer molecules resolving only to this viral gene |
| `viral_molecules_ambiguous_allocated` | Fractional ambiguous molecule mass assigned here |
| `host_virus_ambiguous_molecules` | Molecules compatible with host and viral genes |
| `viral_molecules_total_est` | Unique molecules plus allocated ambiguous mass |
| `viral_molecules_upper_bound` | Unique signal plus all viral-compatible ambiguous mass |
| `n_unique_viral_cells` | Cells with unique viral signal |
| `n_ambiguous_viral_cells` | Cells with viral-compatible ambiguous signal |
| `multimap_method` | `equal`, `host-conservative`, `unique-weighted`, `em-global`, or `em-cell` |
| `evidence_tier` | `candidate_unique`, `candidate_virus_ambiguous`, `candidate_host_virus_ambiguous`, or `not_detected` |

The default `multimap_method` is `host-conservative` (keeps host-virus ambiguous
mass off viral genes); use `equal` for an equal-allocation comparison. These are
molecule-evidence tiers only. `probable` and `strong` require calibrated
read-level evidence and are never assigned from molecule counts alone.

---

<!-- viralscan-claim:v3-evidence-workflow status=provenance_incomplete -->
## `viralscan evidence` output

Evidence is generated in the explicit `--output` directory for one exact
accession, registered alias, or canonical call. With `--viral-fasta` it also
requires a full `--host-fasta`; alignment and BLAST are competitive rather
than virus-only.

The implementation and test sources cover this workflow, but no immutable
successful execution receipt is registered; its execution claim therefore
remains provenance-incomplete.

| Output | Description |
|--------|-------------|
| `read_lineage.tsv.gz` | Exact read ID, CB, UB, EC, compatible genes, ambiguity class, assigned weight, method, tier, and exclusion reason |
| `evidence_manifest.json` | Schema version, exact target, allocation method, run fingerprint, reference hashes, and output hashes |
| `competitive_reads.raw.bam` | Sorted/indexed raw competitive host-plus-target alignments |
| `competitive_reads.<mode>_dedup.bam` | Separate UMI, coordinate-marked, or non-deduplicated evidence BAM |
| `competitive_reads.deduplicated.tagged.bam` | Optional indexed BAM with CB/UB tags for per-cell IGV grouping |
| `coverage.raw.tsv`, `coverage.deduplicated.tsv` | Raw and deduplicated coverage summaries |
| `coverage.raw_vs_deduplicated.png` | Direct depth-track comparison |
| `alignment_qc.tsv` | Per-reference breadth at 1x/3x/10x, depth, intervals, hotspots, strands, mapping quality, identity, host competition, cells, molecules, and duplicate fraction |
| `per_cell_alignment_qc.tsv` | Per-cell host/virus competitive reads, molecules, strands, mapping quality, and identity |
| `blast_identity.tsv` | Best host and viral hit with identity, query coverage, E-value, bit score, score difference, and sequence-complexity flag |
| `blast_sampling.json` | Deterministic sampling strategy, seed, counts, and fraction |
| `interpretation_flags.tsv` | Transparent host-homology, low-complexity, ambiguity, and hotspot diagnostics; all are diagnostic only |
| `viralscan_evidence.igv.xml` | IGV session containing raw, deduplicated, and optional CB/UB-tagged BAMs |

Contamination, expected 3-prime bias, and subgenomic-RNA-like patterns remain
`not_assessed` unless a suitable negative-control or target-specific model is
available. No diagnostic flag is an automatic biological conclusion.

---

## `hostresponse/` — host-response outputs

Written by `viralscan hostresponse` (or a main run with `--host-h5ad`) under
`<output>/hostresponse/`. See the CLI reference for the flags; the key point is
that the raw candidate-support label can be **depth-confounded**, so interpret
the model AUC together with the depth baseline, cohort balance, and sensitivity
metrics.

**`hostresponse_metrics.csv`** — one row per virus:

| Column | Description |
|--------|-------------|
| `virus`, `n_positive` | Virus name; number of candidate-positive cells under the chosen label |
| `label`, `depth_matched`, `mito_controlled` | Which de-confounding design was used (`raw`/`cpm`/`fraction`; depth-matched cohort; %mito covariate) |
| `auc_mean` / `sensitivity_*` / `specificity_*` / `balanced_acc_*` / `mcc_*` | Held-out model metrics (mean/SD over seeds) |
| `depth_alone_auc_mean` | **AUC from sequencing depth ALONE** under the identical split. If this ≈ `auc_mean`, the model AUC is a depth artifact |
| `n_stable_genes`, `n_genes_evalue_ge2` | Stable genes, and how many survive depth adjustment with an E-value ≥ 2 (robust) |
| `n_differential_fdr05` | (with `--differential`) genes significant at FDR < 0.05 in the genome-wide test |

**`<virus>_gene_weights.csv`** — per-fold-HVG L2 coefficients (`weight_mean/sd`,
`n_folds_selected`). **`<virus>_stability.csv`** — per-gene randomized-Lasso
selection probability (`stab_prob`, `stable`). Both gain a `symbol` column with
`--gene-symbols`.

**`<virus>_depth_diagnostics.csv`** — per stable gene, the depth- (and, by
default, %mito-) adjusted odds ratio (`adj_OR`) and its `E_value` (the confounder
strength needed to explain the association away under the stated sensitivity
model). It does not prove absence of residual or unmeasured confounding.

**`<virus>_differential.csv`** (with `--differential`) — a genome-wide,
depth/%mito-adjusted differential table over **all** features: `partial_r`,
`p_value`, `fdr` (Benjamini-Hochberg), `direction` (`up`/`down`).

---

## AnnData files (`.h5ad`)

The AnnData objects can be loaded with [scanpy](https://scanpy.readthedocs.io/):

```python
import scanpy as sc

adata = sc.read_h5ad("output/sample/kb-python/counts_unfiltered/adata_multimap.h5ad")
print(adata)
# Layers: counts_unique, counts_ambiguous_allocated, plus diagnostic layers
```

Key layers:

| Layer | Description |
|-------|-------------|
| `counts_unique` | Bustools-resolved one-gene molecule counts |
| `counts_ambiguous_allocated` | Selected-method ambiguous molecule allocation |
| `counts_multimap_equal` | Equal-split ambiguous molecule allocation |
| `counts_multimap_host_conservative` | Allocation of mixed host-virus molecule mass only among compatible host genes |
| `counts_multimap_unique_weighted` | Heuristic allocation weighted by unique-gene evidence plus pseudocount |
| `counts_unique_viral` | Unambiguous viral molecule evidence retained for auditing |
| `counts_host_viral_ambiguous` | Viral-compatible host-virus ambiguous signal |
| `counts_host_viral_selected` | Portion of the selected allocation assigned to viral genes from host-virus ambiguous ECs |

`adata.X` = `counts_unique + counts_ambiguous_allocated` (complete selected-method matrix).

---

## Multiple samples

When `--sample1` and `--sample2` contain comma-separated FASTQ lists,
ViralScan processes each pair separately:

```bash
viralscan \
  -t t2g.txt -i index.idx -o output/ \
  -s1 A_R1.fastq.gz,B_R1.fastq.gz \
  -s2 A_R2.fastq.gz,B_R2.fastq.gz
```

Expected directories:

```text
output/A/
output/B/
```
