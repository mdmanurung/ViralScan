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
        ├── viral_summary.tsv
        ├── per_cell_viral.tsv
        ├── multimap_evidence.tsv
        ├── cell_type_enrichment.tsv
        └── reference_provenance.json
```

`reference_provenance.json` records the viral reference used (index/t2g/GTF,
technology, multimap settings, and the viral accessions in the reference and
detected) so results are traceable to their annotation.
`cell_type_enrichment.tsv` is present only when `--cell-types` is supplied.
`multimap_evidence.tsv` is present only when multimapping is enabled.
UMAP files are present only when `--umap` is supplied. `host_filtered/` is
present only when `--host-filter` is supplied.

The raw `output.bus` is not a v3 counting input. ViralScan barcode-corrects it
when a whitelist is supplied, sorts it with bustools, and retains the resolved
sorted BUS plus text representation used for molecule counting and read lineage.

The v3 STAR filter also writes `host_filter_audit.tsv` with input, retained,
and removed fragment totals, plus `fragment_lineage.tsv.gz` containing each
retained exact read ID and its `host_unmapped` decision. ViralScan validates
mate synchronization before and after filtering.

---

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
| `pct_infected_called` | `infected_called / n_called_cells × 100`; a called-cell candidate-support rate, not a biological infection rate |

**Two denominators.** The all-barcode `pct_infected` field is diluted by empty
droplets; `pct_infected_called` uses only the declared called-cell set and is the
appropriate denominator for a per-cell candidate-support rate. The field names
are retained for schema compatibility and do not establish infection. The v3
`auto` default uses an external called-cell list when one is supplied and
otherwise runs DropletUtils EmptyDrops. The approximate knee caller is explicit
only and is never a fallback. With `cell_calling=none`, the `*_called` columns
equal the all-barcode values.

**Count layer.** V3 summaries use `adata.X`, the complete selected-method
molecule matrix. `counts_unique` and `counts_ambiguous_allocated` are disjoint
and sum to `X`. Nonzero molecule support is candidate evidence; biological
interpretation requires calibrated evidence and may still require orthogonal
confirmation. The separate read-level workflow supplies diagnostics rather than
an automatic infection call.

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

## `viralscan evidence` output

Evidence is generated in the explicit `--output` directory for one exact
accession, registered alias, or canonical call. With `--viral-fasta` it also
requires a full `--host-fasta`; alignment and BLAST are competitive rather
than virus-only.

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
