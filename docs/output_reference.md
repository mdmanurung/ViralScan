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
        └── cell_type_enrichment.tsv
```

`cell_type_enrichment.tsv` is present only when `--cell-types` is supplied.
`multimap_evidence.tsv` is present only when multimapping is enabled.
UMAP files are present only when `--umap` is supplied. `host_filtered/` is
present only when `--host-filter` is supplied.

---

## `viral_summary.tsv`

Tab-separated, one row per detected virus.

| Column | Description |
|--------|-------------|
| `virus_name` | Human-readable virus name |
| `total_umi` | Total viral UMI across all cells |
| `infected_cells` | Number of cells with ≥ `--detection-threshold` viral UMI |
| `total_cells` | Total cells passing QC |
| `pct_infected` | `infected_cells / total_cells × 100` |
| `umi_per_10k` | `total_umi / total_umi_all × 10 000` |
| `cluster_pvalue` | Fisher's exact test p-value from `viral_neighbor_enrichment` |

---

## `per_cell_viral.tsv`

Tab-separated, one row per cell × detected virus combination.

| Column | Description |
|--------|-------------|
| `barcode` | Cell barcode |
| `virus_name` | Virus name |
| `viral_umi` | Viral UMI count for this cell |
| `total_umi` | Total UMI count for this cell |
| `viral_fraction` | `viral_umi / total_umi` |

---

## `report.html`

A self-contained HTML file with:

- Run metadata (date, parameters, sample names)
- QC summary table
- Per-virus detection table with normalised metrics
- Multimapping evidence table with unique and ambiguous viral support
- Embedded histogram plots (base64 PNG)
- Interpretation guidance
- `viral_neighbor_enrichment` p-values

Open in any modern browser — no internet connection required.

---

## `cell_type_enrichment.tsv`

Tab-separated, one row per detected virus and labeled cell type. Written only
when `--cell-types cell_types.csv` is supplied.

| Column | Description |
|--------|-------------|
| `virus` | Virus name |
| `cell_type` | Cell-type label from the CSV |
| `n_infected` | Infected labeled cells in this cell type |
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
| `unique_viral_umi` | UMI from ECs mapping only to this viral gene |
| `ambiguous_viral_umi` | Selected multimapper share assigned to this viral gene |
| `host_viral_ambiguous_umi` | Viral-compatible UMI from ECs containing host and viral genes |
| `corrected_viral_umi` | `unique_viral_umi + ambiguous_viral_umi` |
| `upper_bound_viral_umi` | Unique signal plus all viral-compatible ambiguous signal |
| `n_unique_viral_cells` | Cells with unique viral signal |
| `n_ambiguous_viral_cells` | Cells with viral-compatible ambiguous signal |
| `multimap_method` | `equal`, `host-conservative`, or `unique-weighted` |
| `call_confidence` | `strong`, `ambiguous`, `low_confidence`, or `not_detected` |

The default `multimap_method` is `host-conservative`. Confidence tiers
prioritize unambiguous viral signal. A `low_confidence` row is supported only
by host-virus ambiguous ECs and should be interpreted cautiously.

---

## AnnData files (`.h5ad`)

The AnnData objects can be loaded with [scanpy](https://scanpy.readthedocs.io/):

```python
import scanpy as sc

adata = sc.read_h5ad("output/sample/kb-python/counts_unfiltered/adata_multimap.h5ad")
print(adata)
# Layers: counts_original, counts_corrected
```

Key layers:

| Layer | Description |
|-------|-------------|
| `counts_original` | Raw kb count matrix (unique-mapping reads) |
| `counts_corrected` | Selected extra multi-mapped read share (additive correction) |
| `counts_multimap_equal` | Equal-split multimapper correction |
| `counts_multimap_host_conservative` | Correction excluding host-virus ambiguous mass from viral genes |
| `counts_multimap_unique_weighted` | Correction weighted by unique-gene evidence plus pseudocount |
| `counts_unique_viral` | Unambiguous viral signal used by `--multimap-primary-call unique-only` |
| `counts_host_viral_ambiguous` | Viral-compatible host-virus ambiguous signal |
| `counts_host_viral_selected` | Portion of the selected correction assigned to viral genes from host-virus ambiguous ECs |

`adata.X` = `counts_original + counts_corrected` (combined count matrix).

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
